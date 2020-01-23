setwd("/Users/katyblumer/repos/SPQV")

# devtools::document('/Users/katyblumer/repos/SPQV')

library(plyr)
library(dplyr)
library(stringr)
# library(SPQV)


placement_type <- "extension"
trait <- "test_trait"



markers <- read.csv("example_data/TeoNAM_Marker_List.csv", stringsAsFactors = FALSE)
gene_list <- read.csv("example_data/Maize_Tiller_Gene_list.csv", stringsAsFactors = FALSE)
qtl <- read.csv("example_data/all_teonam_tb1_traits.csv", stringsAsFactors = FALSE)
chromosome_size <- read.csv("example_data/Maize_Chromosome_Size.csv", stringsAsFactors = FALSE)
wgd <- read.csv("example_data/Maize_WholeGenomeGeneDistribution.csv", stringsAsFactors = FALSE)


num_reps <- 100
num_genes <- 10
qtl_len <- as.integer(147483647)

# Sample genes
sampled_genes <- sample_n(gene_list, num_genes, replace=FALSE)
sampled_genes$ID <- sampled_genes$GeneID
sampled_genes$Base <- as.integer(sampled_genes$Locus)
sampled_genes$Trait <- trait

# Make QTL
valid_qtl_found <- FALSE
while (!valid_qtl_found) {
  start_marker <- sample_n(markers, size = 1)
  end_base <- start_marker[1, 'Base'] + qtl_len

  last_chr_marker <- chromosome_size[
    which(chromosome_size['Number']==start_marker[1, 'Chromosome']),
    'Last.Marker']
  if (end_base > last_chr_marker) {
    next
  } else {
    valid_qtl_found <- TRUE
  }
  sample_qtl <- as.data.frame(
    t(array(
      # Same array twice, since it behaves differently if only 1 qtl?
      c(
        start_marker[1, 'Chromosome'], start_marker[1, 'Base'], end_base,
        trait, "test_treatment", "test_method", "test_expt_type", qtl_len
      ),
      dim = c(length(colnames(qtl)), 1)
    )),
    stringsAsFactors = FALSE
  )
  colnames(sample_qtl) <- colnames(qtl)
}
sample_qtl$Chromosome <- as.integer(sample_qtl$Chromosome)
sample_qtl$LeftmostMarker <- as.integer(sample_qtl$Leftmost_Marker)
sample_qtl$RightmostMarker <- as.integer(sample_qtl$Rightmost_Marker)
sample_qtl$ExptType <- as.character(sample_qtl$Expt_Type)
sample_qtl$Length <- as.integer(sample_qtl$Length)

# Fix chromosome_size
chromosome_size$Chromosome <- chromosome_size$Number
chromosome_size$LeftmostMarker <- chromosome_size$First.Marker
chromosome_size$RightmostMarker <- chromosome_size$Last.Marker

# Fix wgd
wgd$GeneMiddle <- as.integer(wgd$GeneMiddle)


my_env <- new.env()

TEMP_output <- SPQValidate(
  qtl_list = sample_qtl,
  trait = trait,
  num_repetitions = num_reps,
  gene_list = sampled_genes,
  chromosome_size = chromosome_size,
  marker_list = markers,
  whole_genome_gene_dist = wgd,
  placement_type = placement_type,
  simulation_env = my_env
)
TEMP_output
