#Different types of bootstrapping for paper

library(bcaboot)
library(gplots)


# RWR functions ######


sampleUniformQTLLoc <- function(skip_hangovers, bidirectional, qtl_length, chromosome_size) {
  allowed_chr <- chromosome_size[which(chromosome_size$Length >= qtl_length), ]
  sampled_chr <- allowed_chr[sample(nrow(allowed_chr), 1, prob = allowed_chr$Length), ]
  chr_length <- sampled_chr$Length

  if (skip_hangovers) {
    start_base <- sample(chr_length - qtl_length, 1)
  } else {  # bounceback
    start_base <- sample(chr_length, 1)
    if (start_base + qtl_length > chr_length) {
      start_base <- chr_length - qtl_length
    }
  }

  qtl_range <- c(start_base, start_base + qtl_length)

  if (bidirectional) {
    # 50% chance to start range from opposite end
    if (sample(c(TRUE, FALSE), 1)) {
      qtl_range <- c(chr_length - (start_base + qtl_length) + 1, chr_length - start_base + 1)
    }
  }
  return(list(sampled_chr, qtl_range))
}


getAllowedMarkers <- function(
  skip_hangovers,
  qtl_length,
  marker_list
) {
  for (chr in unique(marker_list$Chromosome)) {
    chr_indices <- which(marker_list$Chromosome == chr)
    last_avail_base <- max(marker_list[chr_indices, "Base"]) - qtl_length

    # Either remove disallowed markers or do bounceback
    unavail_equiv_base <- if (skip_hangovers) 0 else last_avail_base
    marker_list[
      which(marker_list$Chromosome == chr & marker_list$Base > last_avail_base),
      "Base"
    ] <- unavail_equiv_base
  }

  marker_list <- marker_list[which(marker_list$Base > 0), ]
  return(marker_list)
}

#' TODO write docs
#'
#' @export
RWR <- function(
  markers_only,  # vs. random placement
  skip_hangovers,  # vs bounceback
  bidirectional,  # vs left-to-right only
  drop_tandem,
  qtl_list,
  gene_list,
  marker_list,
  chromosome_size,
  n_reps
) {

  qtl_list <- validateDf(qtl_list, list(
    c("Chromosome", "integer"),
    c("LeftmostMarker", "integer"),
    c("RightmostMarker", "integer"),
    c("Trait", "character"),
    c("Length", "integer")  # TODO is this superfluous? (or is Rightmost/leftmost? prob not bc need real gene count, but do we return that?)
  ))# need it later, doesn't need to be in the list at the start.

  gene_list <- validateDf(gene_list, list(  # TODO maybe rename as known_gene_list
    c("ID", "character"),
    c("Trait", "character"),
    c("Chromosome", "integer"),
    c("Base", "integer")  # TODO Maybe should be clearer that this is center
  ))

  marker_list <- validateDf(marker_list, list(
    c("ID", "character"),
    c("Chromosome", "integer"),
    c("Base", "integer")
  ))

  chromosome_size <- validateDf(chromosome_size, list(
    c("Chromosome", "integer"),
    c("LeftmostMarker", "integer"),  # TODO Actually we should be able to get this from marker_list - why here? #idk this is a holdover from the dawn of time
    c("RightmostMarker", "integer"),  # TODO superfluous? ditto w/comments below
    c("Length", "integer")
  ))

  dummy_trait <- "trait"

  gene_list$Trait <- dummy_trait
  if (drop_tandem) {
    gene_list <- FilterGeneList(trait = dummy_trait, gene_list = gene_list, marker_list = marker_list)
  }

  qtl_list$Trait <- dummy_trait
  qtl_list <- CountGenesFound(
    qtl_list,
    dummy_trait,
    gene_list,
    marker_list
  )

  qtl_cis <- rep(0, nrow(qtl_list))

  for (qtl_i in 1:nrow(qtl_list)) {
    qtl_length <- qtl_list[qtl_i, "Length"]
    obs_genes_found <- qtl_list[qtl_i, "NumGenes"]

    if (markers_only) {
      allowed_markers <- getAllowedMarkers(skip_hangovers, qtl_length, marker_list)
      allowed_markers <- cbind(allowed_markers, is_l_to_r=TRUE)

      if (bidirectional) {
        rl_markers <- cbind(marker_list, is_l_to_r=FALSE)

        chr_lengths <- merge(rl_markers, chromosome_size[, c("Chromosome", "Length")],
                            by="Chromosome")
        rl_markers$Base <- chr_lengths$Length - rl_markers$Base + 1  # reverse bases
        rl_markers <- getAllowedMarkers(skip_hangovers, qtl_length, rl_markers)
        chr_lengths <- merge(rl_markers, chromosome_size[, c("Chromosome", "Length")],
                             by="Chromosome")
        rl_markers$Base <- chr_lengths$Length - rl_markers$Base + 1  # reverse them back

        allowed_markers <- rbind(allowed_markers, rl_markers)
      }
    }

    exp_genes_found <- rep(0, n_reps)
    for (rep_i in 1:n_reps) {
      if (markers_only) {
        sampled_marker <- allowed_markers[sample(nrow(allowed_markers), 1), ]
        sampled_chr <- sampled_marker$Chromosome

        if (sampled_marker$is_l_to_r) {
          qtl_range <- c(sampled_marker$Base, sampled_marker$Base + qtl_length)
        } else {
          qtl_range <- c(sampled_marker$Base - qtl_length, sampled_marker$Base)
        }

      } else {  # !markers_only (all bases)
        chr_and_range <- sampleUniformQTLLoc(skip_hangovers, bidirectional, qtl_length, chromosome_size)
        sampled_chr <- chr_and_range[[1]]$Chromosome
        qtl_range <- chr_and_range[[2]]
      }

      genes_found <- nrow(gene_list[
          which((gene_list$Chromosome == sampled_chr) &
                gene_list$Base >= qtl_range[1] &
                gene_list$Base <= qtl_range[2])
        ,])
      exp_genes_found[rep_i] <- genes_found

    }  # endfor n_reps

    BCA_CIs <- as.data.frame(as.matrix(
      bcaboot::bcajack(exp_genes_found, n_reps, mean, verbose=F)$lims
    ))

    qtl_cis[qtl_i] <- BCA_CIs["0.95", "bca"]

  } # endfor qtl_list

  return(qtl_cis)
}

# Run RWR and SPQV ####

setwd("/Users/katyblumer/repos/SPQV/")

num_reps <- 1000

trait <- "test_trait"

chromosome_size <- read.csv("example_data/Sviridis_ChromosomeSizes.csv", stringsAsFactors = FALSE)
marker_list <- read.csv("example_data/Sviridis_MarkerList.csv", stringsAsFactors = FALSE)

gene_list <- read.csv("R/FakeHBVGenes333.csv", stringsAsFactors = FALSE)
qtl_list <- read.csv("R/FakeQTL.csv", stringsAsFactors = FALSE)

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

regenerate_results <- FALSE

if (regenerate_results) {
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
          n_reps = num_reps
        )
        row_i <- row_i + 1
      }
    }
  }

  write.csv(x = RWR_results, file = "example_data/RWR_results.csv", row.names=FALSE)



  SPQV_results <- SPQValidate(
    qtl_list = qtl_list,
    trait = trait,
    num_repetitions = num_reps,
    placement_type = "extension",
    gene_list = gene_list,
    marker_list = marker_list,
    wgd,
    chromosome_size = chromosome_size,
    new.env()
  )

  write.csv(x = SPQV_results, file = "example_data/SPQV_results.csv", row.names=FALSE)

}



###### Plot #####

showHeatmap <- function(data_df,
                        num_colors = 13, color_center = 0, color_range = NULL,
                        color_breaks = NULL, color_vals = NULL
                        ) {
  # Add row means
  row_means <- c()
  for(i in 1:nrow(data_df)){
    curr_row <- data_df[i,]
    curr_row <- curr_row[!is.na(curr_row)]
    curr_row <- abs(curr_row)
    row_means <- c(row_means, mean(curr_row))
  }
  # Repeat mean col to make it more visible
  data_df_w_mean <- cbind(data_df, row_means, row_means, row_means)
  colnames(data_df_w_mean)[(ncol(data_df)+1) : ncol(data_df_w_mean)] <- ''
  colnames(data_df_w_mean)[ncol(data_df)+2] <- 'Mean'
  data_df_w_mean <- sapply(data_df_w_mean, as.numeric)

  # Handle colors
  if (is.null(color_breaks)) {
    if (is.null(color_range)) {
      nums <- na.omit(as.numeric(unlist(data_df[1:nrow(data_df)-1,])))
      color_range <- c(min(nums), max(nums))
    }
    color_breaks <- seq(color_range[1],
                       color_range[2],
                       length.out=num_colors+1)
  }
  if (is.null(color_vals)) {
    gradient_g <- colorpanel(
      sum(color_breaks <= color_center), rgb(0,146,146,max=255), "ghostwhite" )
    gradient_p <- colorpanel(
      sum(color_breaks > color_center), "ghostwhite", rgb(73,0,146,max=255) )
    if (length(gradient_g) <= 1) {
      gradient_p <- gradient_p[2:length(gradient_p)]
    }
    color_vals <- unique(c(gradient_g,gradient_p))
  }

  # Handle colnames: Only show 1/10 of QTL lengths on plot, and round & pretty-print them
  col_plot_labels <- colnames(data_df_w_mean)
  for (i in 1:length(col_plot_labels)) {
    if (col_plot_labels[i] == "Mean" | col_plot_labels[i] == "") {
      next
    }
    if (i %% 10 == 1) {	 # always show first one
      qtl_len <- as.integer(col_plot_labels[i])
      qtl_len <- round(qtl_len, -3)
      if (qtl_len > 1e6) {
        qtl_len <- round(qtl_len, -6)
      }
      col_plot_labels[i] <- prettyNum(qtl_len, big.mark=",", scientific=FALSE)

    } else {
      col_plot_labels[i] <- ""
    }
  }

  # Plot
  par(mar=c(10, 8, 8, 3) + 0.1)
  heatmap.2(
    x = as.matrix(data_df_w_mean),
    col = color_vals, breaks = color_breaks,
    srtCol = 45, labCol = col_plot_labels,
    colsep = length(col_plot_labels) - 3, sepcolor = "black",sepwidth=0.5,
    main = 'Comparison of RWR methods to the SPQV',
    denscol = 'black',
    na.color = 'lightslategrey',
    dendrogram='none', trace="none", Colv=F, Rowv=F)
}

# RWR results (and orig figure order):
# 1  no_m | ye_bb | unid | n_dup  1
# 2  no_m | ye_bb | bidi | n_dup  3
# 3  no_m | no_bb | unid | n_dup  2
# 4  no_m | no_bb | bidi | n_dup  4
# 5  ye_m | ye_bb | unid | n_dup  5
# 6  ye_m | ye_bb | bidi | n_dup  7
# 7  ye_m | no_bb | unid | n_dup  6
# 8  ye_m | no_bb | bidi | n_dup  8

RWR_results <- read.csv("example_data/RWR_results.csv")
SPQV_results <- read.csv("example_data/SPQV_results.csv")
qtl_range_to_show <- 1:199
RWR_trunc <- RWR_results[c(1, 3, 2, 4, 5, 7, 6, 8), qtl_range_to_show]
SPQV_trunc <- SPQV_results[qtl_range_to_show,]

colnames(RWR_trunc) <- SPQV_trunc$QTL

method_ratios <- RWR_trunc
for (row_i in 1:nrow(method_ratios)) {
  method_ratios[row_i, ] <- RWR_trunc[row_i, ] / SPQV_trunc[, "Upper.95..CI"]  # Special characters are removed when loading from CSV
}

# Plot
old_colors <- c(
  "#009292", "#3EACAD", "#7CC5C9", "#BADFE4", "#F8F8FF", "#F8F8FF", "#DFD5EF", "#C6B1E0", "#AD8ED0",
  "#946AC1", "#7B47B1", "#6223A2", "#490092")
old_breaks <- c(
  0.3744321, 0.5032979, 0.6321636, 0.7610294, 0.8898951, 1.0187609, 1.1476266, 1.2764924, 1.4053581,
  1.5342239, 1.6630896, 1.7919554, 1.9208211, 2.0496869)
final_plot <- showHeatmap(
  method_ratios,
  color_vals = old_colors, color_breaks = old_breaks
)

