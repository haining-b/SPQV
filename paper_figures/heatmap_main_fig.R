#Different types of bootstrapping for paper

setwd("/Users/katyblumer/repos/SPQV/SPQV")
devtools::document(".")

# Inputs ####
num_reps <- 10

trait <- "test_trait"

chromosome_size <- read.csv("example_data/Sviridis_ChromosomeSizes.csv", stringsAsFactors = FALSE)
marker_list <- read.csv("example_data/Sviridis_MarkerList.csv", stringsAsFactors = FALSE)

gene_list <- read.csv("example_data/FakeHBVGenes333.csv", stringsAsFactors = FALSE)
qtl_list <- read.csv("example_data/FakeQTL.csv", stringsAsFactors = FALSE)

wgd <- read.csv("example_data/allSiGeneswithends.csv", stringsAsFactors = FALSE)

# Regularize
chromosome_size <- dplyr::rename(chromosome_size, LeftmostMarker="First.Markers", RightmostMarker="Last.Markers", Chromosome="Number")

marker_list <- dplyr::rename(marker_list, ID="id")

gene_list <- dplyr::rename(gene_list, ID="GeneID")
gene_list$Trait <- trait
gene_list <- gene_list[, c("ID", "Trait", "Chromosome", "Base")]

qtl_list <- dplyr::rename(qtl_list,
              Chromosome="chromosome", LeftmostMarker="LCI_marker", RightmostMarker="RCI_pos",
              Trait="trait", Treatment="treatment", Method="type",
              ExptType="expt_type")
qtl_list$Length <-as.integer(0)
qtl_list$Trait <- trait
qtl_list <- qtl_list[, c("Chromosome", "LeftmostMarker", "RightmostMarker", "Trait", "Treatment",
               "Method", "ExptType", "Length" )]
qtl_list$Length <-qtl_list$RightmostMarker - qtl_list$LeftmostMarker

wgd$GeneMiddle <- as.integer(wgd$GeneStart + round((wgd$GeneEnd - wgd$GeneStart) /2, 0))

regenerate_results <- TRUE


# Run RWR and SPQV ####
if (regenerate_results) {
  sim_results_env <- new.env()

  RWR_results <- as.data.frame(matrix(nrow=16,ncol=nrow(qtl_list)))

  row_i <- 1
  for (markers_only in c(FALSE, TRUE)) {
    for (skip_hangovers in c(FALSE, TRUE)) {
      for (bidirectional in c(FALSE, TRUE)) {
        drop_tandem <- TRUE
        print(paste(
          format(Sys.time(), "%a %b %d %X %Y"),
          markers_only, skip_hangovers, bidirectional, drop_tandem
        ))
        RWR_results[row_i, ] <- RWR(
          markers_only = markers_only,
          skip_hangovers = skip_hangovers,
          bidirectional = bidirectional,
          drop_tandem = drop_tandem,
          qtl_list = qtl_list,
          gene_list = gene_list,
          marker_list = marker_list,
          chromosome_size = chromosome_size,
          n_reps = num_reps,
          intermediate_output_env = sim_results_env
        )
        row_i <- row_i + 1
      }
    }
  }

  write.csv(x = RWR_results, file = "../paper_figures/output_data/RWR_results.csv", row.names=FALSE)



  SPQV_results <- SPQValidate(
    qtl_list = qtl_list,
    trait = trait,
    num_repetitions = num_reps,
    placement_type = "extension",
    gene_list = gene_list,
    marker_list = marker_list,
    wgd,
    chromosome_size = chromosome_size,
    sim_results_env
  )

  write.csv(x = SPQV_results, file = "../paper_figures/output_data/SPQV_results.csv", row.names=FALSE)

}


###### Plots ----------- #####

# RWR results (and orig figure order):
# 1  no_m | ye_bb | unid | n_dup  1
# 2  no_m | ye_bb | bidi | n_dup  3
# 3  no_m | no_bb | unid | n_dup  2
# 4  no_m | no_bb | bidi | n_dup  4
# 5  ye_m | ye_bb | unid | n_dup  5
# 6  ye_m | ye_bb | bidi | n_dup  7
# 7  ye_m | no_bb | unid | n_dup  6
# 8  ye_m | no_bb | bidi | n_dup  8

orig_figure_order <- data.frame(matrix(nrow=8, ncol=4))
orig_i <- 1
for (markers_only in c(FALSE, TRUE)) {
  for (skip_hangovers in c(FALSE, TRUE)) {
    for (bidirectional in c(FALSE, TRUE)) {
      for (drop_tandem in c(FALSE, TRUE)) {
        orig_figure_order[orig_i, ] <- c(markers_only, skip_hangovers, bidirectional, drop_tandem)
        orig_i <- orig_i + 1
      }
    }
  }
}
new_figure_order <- numeric(16)
new_i <- 1
for (markers_only in c(FALSE, TRUE)) {  ## EDIT ORDER OF THESE
  for (skip_hangovers in c(FALSE, TRUE)) {
    for (drop_tandem in c(FALSE, TRUE)) {
      for (bidirectional in c(FALSE, TRUE)) {
      new_figure_order[new_i] <- which(
        orig_figure_order[1] == markers_only &
          orig_figure_order[2] == skip_hangovers &
          orig_figure_order[3] == bidirectional &
          orig_figure_order[4] == drop_tandem
      )
      new_i <- new_i + 1
      }
    }
  }
}

RWR_results <- read.csv("../paper_figures/output_data/RWR_results.csv", stringsAsFactors = FALSE)
SPQV_results <- read.csv("../paper_figures/output_data/SPQV_results.csv", stringsAsFactors = FALSE)
qtl_range_to_show <- 1:199
RWR_trunc <- RWR_results[new_figure_order, qtl_range_to_show]
SPQV_trunc <- SPQV_results[qtl_range_to_show,]

colnames(RWR_trunc) <- SPQV_trunc$QTL

spqv <- SPQV_trunc[, "Upper.95..CI"] # Special characters are removed when loading from CSV
if ((min(na.omit(spqv)) < 0)| min(na.omit(unlist(RWR_trunc))) < 0) {
  stop("RWR or SPQV upper confidence intervals contain negative values")
}

method_ratios <- RWR_trunc
pct_change_RS <- RWR_trunc
pct_change_SR <- RWR_trunc
pct_diff_avg <- RWR_trunc
pct_diff_max <- RWR_trunc
pct_diff_min <- RWR_trunc
subtract_diff <- RWR_trunc
ceiling_diff <- RWR_trunc
for (row_i in 1:nrow(method_ratios)) {
  rwr <- RWR_trunc[row_i, ]

  method_ratios[row_i, ] <-
    spqv / rwr

  pct_change_SR[row_i, ] <- 100 *
    (spqv - rwr) / rwr

  pct_change_RS[row_i, ] <- 100 *
    (spqv - rwr) / spqv

  pct_diff_max[row_i, ] <- 100 *
    (spqv - rwr) /
    max(max(spqv), max(rwr))  # spqv is list so do max first bc i trust that more than conversion

  # this is a bad idea, just doing it to see what happens
  pct_diff_min[row_i, ] <- 100 *
    (spqv - rwr) /
    min(min(spqv), min(rwr))

  # THIS IS THE ONE WE'RE USING
  pct_diff_avg[row_i, ] <- 100 *
    (spqv - rwr) /
    ((spqv + rwr)/2)

  subtract_diff[row_i, ] <-
    (spqv - rwr)

  ceiling_diff[row_i, ] <-
    (ceiling(spqv) - ceiling(unlist(rwr)))
}

# Repro plot from July ####
old_colors <- c(
  "#009292", "#3EACAD", "#7CC5C9", "#BADFE4", "#F8F8FF", "#F8F8FF", "#DFD5EF", "#C6B1E0", "#AD8ED0",
  "#946AC1", "#7B47B1", "#6223A2", "#490092")
old_breaks <- c(
  0.3744321, 0.5032979, 0.6321636, 0.7610294, 0.8898951, 1.0187609, 1.1476266, 1.2764924, 1.4053581,
  1.5342239, 1.6630896, 1.7919554, 1.9208211, 2.0496869)
showHeatmap(
  method_ratios,
  color_vals = old_colors, color_breaks = old_breaks
)
# July but with pct
old_breaks_pct <- (old_breaks-1)*100
showHeatmap(
  pct_change_SR,
  color_vals = old_colors, color_breaks = old_breaks_pct
)


# Ceiling plot ####
bin_width <- 1
max_val <- 5
color_breaks <- seq(bin_width, max_val, by = bin_width)
color_breaks <- sort(c(color_breaks, -color_breaks, -0.01, 0.01))
showHeatmap(
  ceiling_diff,
  color_breaks = color_breaks
)

showHeatmap(
  subtract_diff,
  num_colors = 41
)


# FINAL PLOT ####

# Each color bin represents a distance of bin_width from 0
bin_width <- 5
max_val <- 50
color_breaks <- seq(bin_width, max_val, by = bin_width)
color_breaks <- sort(c(color_breaks, -color_breaks))

final_plot <- showHeatmap(
  pct_diff_avg,
  color_breaks = color_breaks
)
