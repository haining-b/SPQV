## Results reported in Haining et al.
## Experiments run on version and system:
# platform       x86_64-apple-darwin15.6.0
# arch           x86_64
# os             darwin15.6.0
# system         x86_64, darwin15.6.0
# status
# major          3
# minor          6.2
# year           2019
# month          12
# day            12
# svn rev        77560
# language       R
# version.string R version 3.6.2 (2019-12-12)
# nickname       Dark and Stormy Night


setwd("/Users/katyblumer/repos/SPQV")

# # devtools::document('/Users/katyblumer/repos/SPQV')
# devtools::install_github("olafmersmann/microbenchmarkCore")
# devtools::install_github("olafmersmann/microbenchmark")
# install.packages("Hmisc")
# install.packages("ggplot2")

library(plyr)
library(dplyr)
library(microbenchmark)
library(reshape2)
library(stringr)
library("ggplot2")
# library(SPQV)


placement_type <- "extension"
trait <- "test_trait"

markers <- read.csv("example_data/TeoNAM_Marker_List.csv", stringsAsFactors = FALSE)
qtl <- read.csv("example_data/all_teonam_tb1_traits.csv", stringsAsFactors = FALSE)
chromosome_size <- read.csv("example_data/Maize_Chromosome_Size.csv", stringsAsFactors = FALSE)
wgd <- read.csv("example_data/Maize_WholeGenomeGeneDistribution.csv", stringsAsFactors = FALSE)

# Fix chromosome_size
chromosome_size$Chromosome <- chromosome_size$Number
chromosome_size$LeftmostMarker <- chromosome_size$First.Marker
chromosome_size$RightmostMarker <- chromosome_size$Last.Marker

# Fix wgd
wgd$GeneMiddle <- as.integer(wgd$GeneMiddle)

# Make wgd usable for gene_list
gene_list <- data.frame(wgd)
gene_list$ID <- "TestGeneId"
gene_list$Base <- gene_list$GeneMiddle
gene_list$Trait <- trait


sampleQTL <- function(qtl_len, num_qtl, markers, chromosome_size) {
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
        c(
          start_marker[1, 'Chromosome'], start_marker[1, 'Base'], end_base,
          trait, "test_treatment", "test_method", "test_expt_type", qtl_len
        ),
        dim = c(length(colnames(qtl)), num_qtl)
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


expTitle <- function(var_name, changed_val) {
  mult <- "  "
  if (changed_val / 1e6 >= 1) {
    mult <- "m"
    changed_val <- changed_val / 1e6
  } else if (changed_val / 1e3 >= 1) {
    mult <- "k"
    changed_val <- changed_val / 1e3
  }

  changed_val <- toString(changed_val)
  changed_val <- paste0(
        strrep(" ", 2 * (3 - nchar(changed_val))),
        changed_val
      )
  return(sprintf("%s: %s%s ", var_name, changed_val, mult))
}


my_env <- new.env()

num_benchmark_reps <- 100

# Make experiment parameter list
qtl_nums <- c(10, 3, 1)
bootstrap_nums <- c(100, 10, 1)
gene_nums <- c(50, 10, 5)
marker_nums <- c(10e3, 3e3, 1e3)  # only have 13k, otherwise would be 2e3, 10e3, 100e3
qtl_lens <- c(100e6, 10e6, 1e6)  # max int is ~2e9

vars <- list(qtl_nums, bootstrap_nums, gene_nums, marker_nums, qtl_lens)
var_titles <- c("QTL Count", "Bootstrap Reps", "Known Gene Count", 'Marker Count', 'QTL Length')

base_combo <- c()
middle_val_i <- 2
for (li in vars) {
  base_combo <- c(base_combo, li[middle_val_i])
}


mbms <- c()

mbm_df <- data.frame(1:num_benchmark_reps)
temp_colname <- expTitle(var_titles[[1]], vars[[1]][1])
colnames(mbm_df) <- temp_colname  # Column will be overwritten

for (li_i in 1:length(vars)) {
  var_li <- vars[[li_i]]
  var_title <- var_titles[[li_i]]
  for (val_i in 1:length(var_li)) {

    changed_val <- var_li[val_i]
    var_combo <- base_combo
    var_combo[li_i] <- changed_val

    colname <- expTitle(var_title, changed_val)

    num_qtl <- var_combo[1]
    num_bootstraps <- var_combo[2]
    num_genes <- var_combo[3]
    num_markers <- var_combo[4]
    qtl_len <- var_combo[5]

    mbm <- microbenchmark(
      spqv_and_sample = {
        sampled_genes <- sample_n(gene_list, num_genes, replace=FALSE)
        sampled_markers <- sample_n(markers, num_markers, replace=FALSE)
        sample_qtl <- sampleQTL(qtl_len, num_qtl, markers, chromosome_size)
        SPQValidate(
          qtl_list = sample_qtl,
          trait = trait,
          num_repetitions = num_bootstraps,
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
        sample_qtl <- sampleQTL(qtl_len, num_qtl, markers, chromosome_size)
      },
      times = num_benchmark_reps
    )

    mbms <- c(mbms, mbm)

    # Subtract median sample prep time from all runs
    mbm_df[colname] <-  (
      mbm[which(mbm$expr=="spqv_and_sample"), 'time'] -
        median(mbm[which(mbm$expr=="sample_only"), 'time'])
    )
  }
}

plot_df <- melt(mbm_df)
colnames(plot_df) <- c("Parameters", "Nanoseconds")
plot_df$Milliseconds <- (plot_df$Nanoseconds / 1e6)
p <- ggplot(plot_df, aes(x=Parameters, y=Milliseconds)) + geom_violin() + coord_flip()
p
p + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + ylim(-1, 50)
p + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  ylab("Milliseconds")+
  # ylim(-1, 1000)+
  xlab("Parameter Varied")+
  theme_classic()
#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

medians <- data.frame(apply(mbm_df, MARGIN = 2, FUN = median))
colnames(medians) <- c('median')
print(medians)
