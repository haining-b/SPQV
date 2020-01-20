setwd("/Users/katyblumer/repos/SPQV")

devtools::document('/Users/katyblumer/repos/SPQV')

library(plyr)
library(dplyr)
library(stringr)
library(SPQV)


placement_type <- "extension"
trait <- "test_trait"



markers <- read.csv("example_data/TeoNAM_Marker_List.csv", stringsAsFactors = FALSE)
gene_list <- read.csv("example_data/Maize_Tiller_Gene_list.csv", stringsAsFactors = FALSE)
qtl <- read.csv("example_data/all_teonam_tb1_traits.csv", stringsAsFactors = FALSE)
chromosome_size <- read.csv("example_data/Maize_Chromosome_Size.csv", stringsAsFactors = FALSE)
wgd <- read.csv("example_data/Maize_WholeGenomeGeneDistribution.csv", stringsAsFactors = FALSE)


num_reps <- 5
num_genes <- 10
qtl_len <- 1e6

# Sample genes
sampled_genes <- sample_n(gene_list, num_genes, replace=FALSE)

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
      c(c(
        start_marker[1, 'Chromosome'], start_marker[1, 'Base'], end_base,
        "test_trait", "test_treatment", "test_method", "test_expt_type", qtl_len
      ), 
      c(
        start_marker[1, 'Chromosome'], start_marker[1, 'Base'], end_base,
        "test_trait", "test_treatment", "test_method", "test_expt_type", qtl_len
      ) ),
      dim = c(length(colnames(qtl)), 2)
    )),
    stringsAsFactors = FALSE
  )
  colnames(sample_qtl) <- colnames(qtl)
}
sample_qtl$Chromosome <- as.integer(sample_qtl$Chromosome)
sample_qtl$Leftmost_Marker <- as.integer(sample_qtl$Leftmost_Marker)
sample_qtl$Rightmost_Marker <- as.integer(sample_qtl$Rightmost_Marker)
sample_qtl$Length <- as.integer(sample_qtl$Length)


my_env <- new.env()

TEMP_output <- SPQValidate(
  QTL_with_Metadata = sample_qtl,
  Trait = trait,
  number_repetitions = num_reps,
  gene_list = sampled_genes,
  chromosome_size = chromosome_size,
  MarkerList = markers,
  WholeGenomeGeneDistribution = wgd,
  Placement_Type = placement_type,
  Simulation_Environment = my_env
)