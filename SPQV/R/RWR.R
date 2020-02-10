#Different types of bootstrapping for paper

# # paperValidator<-Relevant_QTL[,c(1,2,3,5,4,6)]
# Relevant_QTL <- read.csv("R/FakeQTL.csv",stringsAsFactors = F)
# colnames(paperValidator)[1]<-'chromosome'
#
# BootStrapConfidenceIntervals<-as.data.frame(matrix(nrow=24,ncol=length(Relevant_QTL$ChrNum)))
# colnames(BootStrapConfidenceIntervals)<-paperValidator$Length
# library(bootstrap)
#
# target_gene_list_DUPS<-unique(read.csv("IonomicGeneList.csv",stringsAsFactors = F)[,c(1,3,4)])

library(bcaboot)
library(gplots)

#' TODO write docs
#'
#' @export



# RWR ######

#' TODO write docs
#'
#' @export
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

#' TODO write docs
#'
#' @export
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

# Run RWR ####

setwd("/Users/katyblumer/repos/SPQV/")

num_reps <- 100

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


BSCI <- as.data.frame(matrix(nrow=16,ncol=nrow(qtl_list)))

row_i <- 1
for (markers_only in c(FALSE, TRUE)) {
  for (skip_hangovers in c(FALSE, TRUE)) {
    for (bidirectional in c(FALSE, TRUE)) {
      for (drop_tandem in c(FALSE, TRUE)) {
        print(paste(markers_only, skip_hangovers, bidirectional, drop_tandem))
        BSCI[row_i, ] <- RWR(
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
        print(row_i)
        row_i <- row_i + 1
        print(BSCI[, 100:110])
      }
    }
  }
}

write.csv(x = BSCI, file = "example_data/RWR_results.csv", row.names=FALSE)


SPQV_results <- SPQValidate(
  qtl_list = qtl_list,
  trait = trait,
  num_repetitions = 10, # DO NOT SUBMIT
  placement_type = "extension",
  gene_list = gene_list,
  marker_list = marker_list,
  wgd,
  chromosome_size = chromosome_size,
  new.env()
)




###### Plot ? #####

showHeatmap <- function(ratios, num_colors= 13, color_center=0, color_range=NULL) {
  if (is.null(color_range)) {
    nums <- na.omit(as.numeric(unlist(ratios)))
    color_range <- c(min(nums), max(nums))
  }


  RowMeanPerc<-c()
  for(r in 1:nrow(ratios)){
    CurrentRow<-ratios[r,]
    CR<-CurrentRow[!is.na(CurrentRow)]
    CR<-abs(CR)
    RowMeanPerc<-c(RowMeanPerc,(mean(CR)))
  }
  ratios<-cbind(ratios,RowMeanPerc)

  color_breaks = seq(color_range[1],
                     color_range[2],
                     length.out=num_colors+1)
  gradient_g = colorpanel( sum( color_breaks<=color_center), rgb(0,146,146,max=255), "ghostwhite" )
  gradient_p = colorpanel( sum( color_breaks>color_center ), "ghostwhite", rgb(73,0,146,max=255) )
  print(c(gradient_g,gradient_p))
  if (length(gradient_g) <= 1) {
    gradient_p <- gradient_p[2:length(gradient_p)]
  }
  heatmap_colors = unique(c(gradient_g,gradient_p))
  print(color_breaks)
  print(heatmap_colors)

  par(mar=c(10, 8, 8, 3) + 0.1)
  heatmap.2(x=as.matrix(ratios),dendrogram='none',
            trace="none",Colv=F,Rowv=F,col=heatmap_colors,srtCol = 90,
            breaks = color_breaks,
            #cellnote = round(numHBVSinglesDifference), notecol = 'black',notecex=.7,
            main='Comparison of Bootstrapping methods to the Haining-Blumer Validator; Percent',
            denscol='black',
            na.color='lightslategrey')
}

#' 1  no_m | ye_bb | unid | y_dup
#' 2  no_m | ye_bb | unid | n_dup  1
#' 3  no_m | ye_bb | bidi | y_dup
#' 4  no_m | ye_bb | bidi | n_dup  3
#' 5  no_m | no_bb | unid | y_dup
#' 6  no_m | no_bb | unid | n_dup  2
#' 7  no_m | no_bb | bidi | y_dup
#' 8  no_m | no_bb | bidi | n_dup  4
#' 9  ye_m | ye_bb | unid | y_dup
#' 10 ye_m | ye_bb | unid | n_dup  5
#' 11 ye_m | ye_bb | bidi | y_dup
#' 12 ye_m | ye_bb | bidi | n_dup  7
#' 13 ye_m | no_bb | unid | y_dup
#' 14 ye_m | no_bb | unid | n_dup  6
#' 15 ye_m | no_bb | bidi | y_dup
#' 16 ye_m | no_bb | bidi | n_dup  8

BSCI <- read.csv("example_data/RWR_results.csv")
qtl_range_to_show <- 1:199
BSCI_trunc <- BSCI[c(2, 6, 4, 8, 10, 14, 12, 16), qtl_range_to_show]
SPQV_trunc <- SPQV_results[qtl_range_to_show, "Upper 95% CI"]

method_ratios <- BSCI_trunc
method_sub <- BSCI_trunc


for (row_i in 1:nrow(BSCI_trunc)) {
  method_ratios[row_i, ] <- BSCI_trunc[row_i, ] / SPQV_trunc
  method_sub[row_i, ] <- BSCI_trunc[row_i, ] - SPQV_trunc
}

showHeatmap(method_ratios, color_center=1, color_range = c(0.2, 2))
showHeatmap(log(method_ratios), color_range = c(-2, 2))

showHeatmap(method_sub)
showHeatmap(log(method_sub))


## Inspect individual dfs ####

showHeatmap(BSCI_trunc, color_center = 1, color_range = c(0, 58))
showHeatmap(BS_333, color_center = 1, color_range = c(0, 58))

showHeatmap(SPQV_results, color_center = 1, color_range = c(0, 58))
showHeatmap(HBV_333, color_center = 1, color_range = c(0, 58))


SPQV_ratio <- SPQV_results[qtl_range_to_show, 1:6]/HBV_333[, 1:6]
showHeatmap(
  as.data.frame(SPQV_ratio),
  color_center = 1,
  color_range = c(-0.1, 5)
  )

colnames(BS_ratio)
BS_ratio <- BSCI_trunc[, qtl_range_to_show] - BS_333[, 1:199]
showHeatmap(
  as.data.frame(BS_ratio),
  color_center = 0,
  color_range = c(-3, 3)
  )


colnames(HBV_333)
colnames(SPQV_results[1:199,])



# Try BBReso ####
SPQV_trunc_allcols <- SPQV_results[qtl_range_to_show, ]
PercentDiffMatrix<-as.data.frame(matrix(nrow=8,ncol=nrow(SPQV_trunc_allcols)))
colnames(PercentDiffMatrix)<-SPQV_trunc_allcols$QTL
for(shannon in 1:length(SPQV_trunc_allcols$QTL)){
  percentDiff<-BSCI_trunc[,shannon]/SPQV_trunc_allcols[shannon,"Upper 95% CI"] # 6
  PercentDiffMatrix[,shannon]<-percentDiff
}


HBVSinglesDifferenceMatrix<-as.data.frame(matrix(nrow=8,ncol=ncol(BSCI_trunc)))
# colnames(HBVSinglesDifferenceMatrix)<-Fakes$Length
for (exp_i in 1:ncol(BSCI_trunc)){
  DifferenceFromHBV=BSCI_trunc[,exp_i]-BSCI_trunc[exp_i,6]
  HBVSinglesDifferenceMatrix[,exp_i]<-DifferenceFromHBV
}
HBVSinglesDifferenceMatrix <- HBVSinglesDifferenceMatrix[, 1:ncol(HBVSinglesDifferenceMatrix)-1]


library(gplots)
breaks = seq(min(na.omit(as.numeric(unlist(PercentDiffMatrix)))),
             2.1, #max(abs(na.omit(as.numeric(unlist(PercentDiffMatrix))))),
             length.out=13)
gradient1 = colorpanel( sum( breaks<=1), rgb(0,146,146,max=255), "ghostwhite" )
gradient2 = colorpanel( sum( breaks>1 ), "ghostwhite", rgb(73,0,146,max=255) )
hm.colors = unique(c(gradient1,gradient2))
#Getting row means
RowMeanPerc<-c()
for(r in 1:8){
  CurrentRow<-PercentDiffMatrix[r,]
  CR<-CurrentRow[!is.na(CurrentRow)]
  CR<-abs(CR)
  RowMeanPerc<-c(RowMeanPerc,(mean(CR)))
}
PercentDiffMatrix<-cbind(PercentDiffMatrix,RowMeanPerc)
numPercentDiffMatrix<-sapply(PercentDiffMatrix,as.numeric)

old_breaks <- c(
  0.3744321, 0.5140367, 0.6536413, 0.7932458, 0.9328504, 1.0724549, 1.2120595, 1.3516641, 1.4912686,
  1.6308732, 1.7704777, 1.9100823, 2.0496869)
old_colors <- c(
  "#009292", "#3EACAD", "#7CC5C9", "#BADFE4", "#F8F8FF", "#F8F8FF", "#DFD5EF", "#C6B1E0", "#AD8ED0",
  "#946AC1", "#7B47B1", "#6223A2", "#490092")

colnames(numPercentDiffMatrix)[200]<-'Mean'
par(mar=c(10, 8, 8, 3) + 0.1)
heatmap.2(numPercentDiffMatrix,dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=old_colors,srtCol = 45,
          breaks=old_breaks,
          #cellnote = round(numHBVSinglesDifference), notecol = 'black',notecex=.7,
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator; Percent',
          denscol='black',
          na.color='lightslategrey')
