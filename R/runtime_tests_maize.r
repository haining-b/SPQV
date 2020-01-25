setwd("/Users/katyblumer/repos/SPQV")

# devtools::document('/Users/katyblumer/repos/SPQV')
# devtools::install_github("olafmersmann/microbenchmarkCore")
# devtools::install_github("olafmersmann/microbenchmark")
# install.packages("Hmisc")

library(plyr)
library(dplyr)
library(microbenchmark)
library(reshape2)
library(stringr)
# library(SPQV)


placement_type <- "extension"
trait <- "test_trait"

markers <- read.csv("example_data/TeoNAM_Marker_List.csv", stringsAsFactors = FALSE)
gene_list <- read.csv("example_data/Maize_Tiller_Gene_list.csv", stringsAsFactors = FALSE)
qtl <- read.csv("example_data/all_teonam_tb1_traits.csv", stringsAsFactors = FALSE)
chromosome_size <- read.csv("example_data/Maize_Chromosome_Size.csv", stringsAsFactors = FALSE)
wgd <- read.csv("example_data/Maize_WholeGenomeGeneDistribution.csv", stringsAsFactors = FALSE)

# Fix gene_list
gene_list$ID <- gene_list$GeneID
gene_list$Base <- as.integer(gene_list$Locus)
gene_list$Trait <- trait

# Fix chromosome_size
chromosome_size$Chromosome <- chromosome_size$Number
chromosome_size$LeftmostMarker <- chromosome_size$First.Marker
chromosome_size$RightmostMarker <- chromosome_size$Last.Marker

# Fix wgd
wgd$GeneMiddle <- as.integer(wgd$GeneMiddle)


sampleQTL <- function(qtl_len, markers, chromosome_size) {
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

  return(sample_qtl)
}


my_env <- new.env()

num_reps <- 100
gene_nums <- c(5, 10, 45)  # only have 47, otherwise would be 50
marker_nums <- c(1e3, 3e3, 10e3)  # only have 13k, otherwise would be 2e3, 10e3, 100e3
qtl_lens <- c(1e6, 10e6, 100e6)  # max int is ~2e9



mbms <- c()

mbm_df <- data.frame(1:100)
temp_colname <- sprintf("%s_%s_%s", gene_nums[[1]], marker_nums[[1]], qtl_lens[[1]])
colnames(mbm_df) <- temp_colname

for (num_genes in gene_nums) {
  for (num_markers in marker_nums) {
    for (qtl_len in qtl_lens) {
      sampled_genes <- sample_n(gene_list, num_genes, replace=FALSE)
      sampled_markers <- sample_n(markers, num_markers, replace=FALSE)
      sample_qtl <- sampleQTL(qtl_len, markers, chromosome_size)

      mbm <- microbenchmark(spqv_and_sample = {
        sampled_genes <- sample_n(gene_list, num_genes, replace=FALSE)
        sampled_markers <- sample_n(markers, num_markers, replace=FALSE)
        sample_qtl <- sampleQTL(qtl_len, markers, chromosome_size)
        SPQValidate(
        qtl_list = sample_qtl,
        trait = trait,
        num_repetitions = num_reps,
        gene_list = sampled_genes,
        chromosome_size = chromosome_size,
        marker_list = markers,
        whole_genome_gene_dist = wgd,
        placement_type = placement_type,
        simulation_env = my_env,
        progress_bar = FALSE
        )
      },
        sample_only = {
          sampled_genes <- sample_n(gene_list, num_genes, replace=FALSE)
          sampled_markers <- sample_n(markers, num_markers, replace=FALSE)
          sample_qtl <- sampleQTL(qtl_len, markers, chromosome_size)
        }
      )
      mbms <- c(mbms, mbm)

      colname <- sprintf("%s_%s_%s", num_genes, num_markers, qtl_len)
      mbm_df[colname] <-  (
        mbm[which(mbm$expr=="spqv_and_sample"), 'time'] -
          median(mbm[which(mbm$expr=="sample_only"), 'time'])
      )

    }
  }
}

plot_df <- melt(mbm_df)
colnames(plot_df) <- c("Genes.Markers.QTL", "Milliseconds")
plot_df$Milliseconds <- plot_df$Milliseconds / 1e6  # Original is in nanoseconds
p <- ggplot(plot_df, aes(x=Genes.Markers.QTL, y=Milliseconds)) + geom_violin() + coord_flip()
p + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + ylim(-1, 50)
p + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + ylim(-1, 30)

medians <- data.frame(apply(mbm_df, MARGIN = 2, FUN = median))
colnames(medians) <- c('median')
print(medians)
