# SETUP ####
library(tidyverse)
library(janitor)
library(dada2)
library(phyloseq)
source("./R/plot_bar2.R")
source("./R/bbdml_helper.R")
source("./R/palettes.R")

# clean up metadata ####
df <- 
  read_delim("./data/Sponge_Metadata.txt") %>% clean_names()
df <- 
  df[-1,]
df <- 
  df %>% mutate(pna = case_when(pna=="Y" ~ TRUE,
                              TRUE~FALSE),
              nano_drop = as.numeric(nano_drop))

# PARSE FILE PATHS ####

# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "./data/fastq" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "_R1_001.fastq.gz"))

# order files to match metadata
names(fns) <- fns %>% str_split("/") %>% map_chr(4) %>% str_split("_") %>% map_chr(1)
identical(df$sample_id,names(fns))
sample.names <- names(fns)

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
plotQualityProfile(fns[1:2])

# FILTER AND TRIM ####
filts <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
# fns <- unname(fns)

out <- filterAndTrim(fns, filts, # fnRs, filtRs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=4) # On Windows set multithread=FALSE


# sanity check  comparison of before and after filtration
plotQualityProfile(c(fns[1:2],filts[1:2]))

# LEARN ERROR RATES ####
# Since some samples may have had zero reads pass QC, reassign filts
filts <- sort(list.files(filtpath, full.names = TRUE))
errF <- learnErrors(filts, multithread=TRUE, MAX_CONSIST = 20)

# sanity check for error model
plotErrors(errF, nominalQ=TRUE)

# DEREPLICATION ####
derep <- derepFastq(filts, verbose=TRUE)


# Name the derep-class objects by the sample names
# If some samples were removed (no reads passed QC), reassign sample.names
sample.names <- unlist(map(strsplit(basename(filts), "_filt"), 1))
names(derep) <- sample.names

# SAMPLE INFERRENCE ####
dadaFs <- dada(derep, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")

# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(dadaFs)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[,1], sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track[,2])/track[,1]
head(track)
write.csv(track, file = "./output/read_counts_at_each_step.csv", row.names = TRUE)


# Save intermediate seqtable object
saveRDS(seqtab.nochim, "./output/seqtab.nochim.RDS")


# INSPECT METADATA ORDER ####
row.names(df) <- as.character(df$sample_id)
identical(row.names(df),row.names(seqtab.nochim))

# ASSIGN TAXONOMY ####
taxa <- assignTaxonomy(seqtab.nochim, "../Clayton_SRA/taxonomy/sh_general_release_dynamic_all_16.10.2022.fasta.gz", multithread=20)

# Save intermediate files
saveRDS(taxa, file = "./output/RDP_Taxonomy_from_dada2.RDS")

# export seqs for BLAST profiling
ASVs <- row.names(taxa)
names(ASVs) <- paste0("ASV_",seq_along(ASVs))
ShortRead::writeFasta(ASVs,"./output/ASV_Seqs.fasta")

# Run BLAST externally
# Import BLAST results (top hits)
blast <- read_delim("./output/UNITE_EUK_BLAST_OUT.txt",col_names = FALSE) %>% 
  select(X1,X2) %>% 
  mutate(ASV = X1,
         TopHit = X2 %>% str_split("\\|") %>% map_chr(5)) %>% 
  select(ASV, TopHit)
blast <- 
  blast %>% 
  separate(TopHit, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                   sep=";") %>% 
  select(-ASV)
row.names(blast) <- paste0("ASV_",1:nrow(blast))
tax2 <- tax_table(blast)


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(df)
row.names(met) <- row.names(df)


ps <- phyloseq(otu,met,tax)


# Find non-fungi
ps_nonfungi <- subset_taxa(ps, Kingdom != "k__Fungi")
ps2 <- ps %>% 
  subset_taxa(!is.na(Kingdom) & !is.na(Phylum))

ps2 %>% 
  merge_samples("pna") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Phylum") +
  labs(x="PNA") +
  theme_minimal() +
  scale_fill_manual(values=pal.discrete)

ggsave("./output/Phylum_Level_Taxonomic_Proportions.png",dpi=300)


ps_pnaTRUE <- 
  subset_samples(ps,pna==TRUE) %>% 
  subset_taxa(taxa_sums(ps_pnaTRUE) > 0)
ps_pnaFALSE <- 
  subset_samples(ps,pna==FALSE) %>% 
  subset_taxa(taxa_sums(ps_pnaFALSE) > 0)

PNA_Kingdoms <- full_join(
  ps_pnaTRUE@tax_table[,1] %>% table(useNA = "always") %>% as.data.frame() %>% mutate(PNA = TRUE,Proportion = Freq/sum(Freq)),
  ps_pnaFALSE@tax_table[,1] %>% table(useNA = "always") %>% as.data.frame() %>% mutate(PNA = FALSE,Proportion = Freq/sum(Freq))
)

names(PNA_Kingdoms) <- c("Kingdom","Freq","PNA","Proportion")
proportions <- 
  PNA_Kingdoms %>% 
  mutate(Kingdom = case_when(is.na(Kingdom) ~ "unidentified",
                             TRUE ~ Kingdom %>% str_remove("k__")))
  
kableExtra::kable(proportions) %>% 
  kableExtra::kable_classic(lightable_options = "hover")

proportions %>% 
  ggplot(aes(x=PNA,y=Proportion,fill=Kingdom)) +
  geom_col()


PNA_Kingdoms %>% 
  mutate(Kingdom = Kingdom %>% str_remove("k__"),
         Freq = Freq / ntaxa(ps)) %>% 
  ggplot(aes(x=PNA,y=Freq,fill=Kingdom)) +
  geom_col() 

# NMDS

# get rid of non-fungi
ps <- ps %>% subset_taxa(Kingdom == "k__Fungi")
ps <- ps %>% subset_taxa(taxa_sums(ps) > 0)

ord <- ordinate(ps,method = "NMDS")
p1 <- plot_ordination(ps,ord,color = "pna") + theme_minimal()
p2 <- plot_ordination(ps,ord,color = "genus") + theme_minimal() + labs(y="")

library(patchwork)
p1 + p2

# permanova
mat <- as.matrix(otu_table(ps))

vegan::adonis2(mat ~ ps@sam_data$genus * ps@sam_data$pna) %>% 
  broom::tidy() %>% 
  mutate(term = term %>% str_remove_all("ps@sam_data\\$")) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic(lightable_options = "hover")

# # REMOVE NON-FUNGI and empty samples/taxa ####
# ps <- subset_taxa(ps, Kingdom == "k__Fungi")
# ps <- subset_taxa(ps, taxa_sums(ps) > 0)
# ps <- subset_samples(ps, sample_sums(ps) > 0)
# 
# # Save DNA sequences apart from rownames (from subsetted ps object)
# seqs <- taxa_names(ps)
# seqs <- DNAStringSet(seqs)
# saveRDS(seqs,"./output/ASV_reference_sequences.RDS")
# 
# 
# pretty_names <- paste("FungalASV",1:length(taxa_names(ps)),":",
#                       tax_table(ps)[,2],
#                       tax_table(ps)[,3],
#                       tax_table(ps)[,4],
#                       tax_table(ps)[,5],
#                       tax_table(ps)[,6],
#                       tax_table(ps)[,7], sep="_") %>%
#   str_remove("k__") %>% str_remove("p__") %>% str_remove("c__") %>% str_remove("o__") %>% str_remove("f__") %>% str_remove("g__") %>% str_remove("s__") %>%
#   str_replace(pattern = "_:_",replacement = ": ")
# 
# df <- data.frame(TaxaName=pretty_names,Sequence=taxa_names(ps))
# saveRDS(df,"./output/SequenceNames_and_Taxonomy.RDS")
# 
# # Set Seawater as first level of Sponge_Species
# ps@sam_data$Sponge_Species <- factor(ps@sam_data$Sponge_Species, 
#                                      levels = c("Seawater","Chondrilla","Chondrosia","Crambe","Petrosia"))
# 
# 
# # Save RDS object for Phyloseq
# saveRDS(ps, file = "./output/clean_phyloseq_object.RDS")
# 
# 
