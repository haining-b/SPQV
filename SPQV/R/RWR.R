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

setwd("/Users/katyblumer/repos/SPQV/")

num_reps <- 100

trait <- "test_trait"

chromosome_size <- read.csv("example_data/Sviridis_ChromosomeSizes.csv", stringsAsFactors = FALSE)
marker_list <- read.csv("example_data/Sviridis_MarkerList.csv", stringsAsFactors = FALSE)

gene_list <- read.csv("example_data/forHBVPaper333UsedGenesWITHNAMES.csv", stringsAsFactors = FALSE)
qtl_list <- read.csv("example_data/FakeQTL.csv", stringsAsFactors = FALSE)

wgd <- read.csv("example_data/allSiGeneswithends.csv", stringsAsFactors = FALSE)

# Regularize
chromosome_size <- dplyr::rename(chromosome_size, LeftmostMarker="First.Markers", RightmostMarker="Last.Markers", Chromosome="Number")

marker_list <- dplyr::rename(marker_list, ID="id")

gene_list <- dplyr::rename(gene_list, ID="Gene.Name", Base="Midpoint")
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

write.csv(x = BSCI, file = "example_data/RWR_results.csv")

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

# ####
qtl_range_to_show <- 100:199
method_ratios <- BSCI[1:8, qtl_range_to_show+1]
SPQV_compare <- SPQV_results[qtl_range_to_show, "Upper 95% CI"]
for (row_i in 1:nrow(method_ratios)) {
  method_ratios[row_i, ] <- method_ratios[row_i, ] / SPQV_compare
}
method_ratios <- log(method_ratios)


###### Plot ? #####
#Getting row means
RowMeanPerc<-c()
for(r in 1:nrow(method_ratios)){
  CurrentRow<-method_ratios[r,]
  CR<-CurrentRow[!is.na(CurrentRow)]
  CR<-abs(CR)
  RowMeanPerc<-c(RowMeanPerc,(mean(CR)))
}
method_ratios<-cbind(method_ratios,RowMeanPerc)

nums <- na.omit(as.numeric(unlist(method_ratios)))

# num_colors = 5
# color_breaks_gr = seq(
#   min(nums[nums <= 0]),
#   max(nums[nums <= 0]),
#   length.out=num_colors+1)
# color_breaks_pu = seq(
#   min(nums[nums > 0]),
#   max(nums[nums > 0]),
#   length.out=num_colors+1)
# breaks <- c(color_breaks_gr[1:num_colors], 0, color_breaks_pu[2:num_colors+1])
# gradient_colors = colorpanel(length(breaks)-1, low=rgb(0,146,146,max=255), mid="ghostwhite", high=rgb(73,0,146,max=255) )

color_breaks = seq(-3,
             3,
             length.out=50)
gradient_g = colorpanel( sum( color_breaks<=0), rgb(0,146,146,max=255), "ghostwhite" )
gradient_p = colorpanel( sum( color_breaks>0 ), "ghostwhite", rgb(73,0,146,max=255) )
heatmap_colors = unique(c(gradient_g,gradient_p))


par(mar=c(10, 8, 8, 3) + 0.1)
heatmap.2(x=TEMP_method_ratios,dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=heatmap_colors,srtCol = 45,
          breaks = color_breaks,
          #cellnote = round(numHBVSinglesDifference), notecol = 'black',notecex=.7,
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator; Percent',
          denscol='black',
          na.color='lightslategrey')

#
#
#
#
# # Plot rounded #####
# RoundPercentDiffMatrix<-as.data.frame(matrix(nrow=8,ncol=length(HBV_333$QTL)))
# colnames(RoundPercentDiffMatrix)<-HBV_333$QTL
#
# for(shannon in 1: length(HBV_333$QTL)){
#   for(pants in 1:8){
#     BS<-round(as.numeric(BS_333[pants,shannon]))
#     HBV<-round(as.numeric(HBV_333[shannon,6]))
#
#     if(BS==1 &HBV==0){
#       print("inf!")
#       RoundPercentDiffMatrix[pants,shannon]<-NA
#     }
#     if(BS==0 &HBV==1){
#       print("WTF")
#       RoundPercentDiffMatrix[pants,shannon]<-NA
#     }
#     if(BS==0 & HBV==0){
#       print("0!")
#       RoundPercentDiffMatrix[pants,shannon]<-1
#     }
#     else{
#       RoundPercentDiffMatrix[pants,shannon]<-BS/HBV
#     }
#   }
#
# }
# RoundPercentDiffMatrix[ sapply(RoundPercentDiffMatrix,is.infinite)]<-NA
# color_breaks = seq(min(na.omit(as.numeric(unlist(RoundPercentDiffMatrix)))),
#              max(abs(na.omit(as.numeric(unlist(RoundPercentDiffMatrix))))),
#              length.out=10)
# gradient1 = colorpanel( sum( color_breaks<1), rgb(0,146,146,max=255), "ghostwhite" )
# gradient2 = colorpanel( sum( color_breaks>1 ), "ghostwhite", rgb(73,0,146,max=255) )
# hm.colors = c(gradient1,gradient2)
#
# rRowMeanPerc<-c()
# for(r in 1:8){
#   CurrentRow<-RoundPercentDiffMatrix[r,]
#   CR<-CurrentRow[!is.na(CurrentRow)]
#   CR<-abs(CR)
#   rRowMeanPerc<-c(rRowMeanPerc,(mean(CR)))
# }
# RoundPercentDiffMatrix<-cbind(RoundPercentDiffMatrix,rRowMeanPerc)
# numRoundPercentDiffMatrix<-sapply(RoundPercentDiffMatrix,as.numeric)
#
# heatmap.2 (numRoundPercentDiffMatrix,dendrogram='none',
#            trace="none",Colv=F,Rowv=F,col=hm.colors,srtCol = 45,
#            #cellnote = round(numHBVSinglesDifference), notecol = 'black',notecex=.7,
#            main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
#            denscol='black',
#            na.color='lightslategrey')
#
#
#
#
#
#
#
#
#
#
# #w/ rounded nums ####
#
# color_breaks = seq(min(na.omit(as.numeric(unlist(ceiling(HBVSinglesDifferenceMatrix))))),
#              max(abs(na.omit(as.numeric(unlist(ceiling(HBVSinglesDifferenceMatrix)))))),
#              length.out=9)
# gradient1 = colorpanel( sum( color_breaks[-1]<=0)+1, "deeppink4", "floralwhite" )
# gradient2 = colorpanel( sum( color_breaks[-1]>0 )-1, "floralwhite", "aquamarine4" )
# hm.colors = c(gradient1,gradient2)
# formatC(numb, format = "e", digits = 2)
# colnames(numHBVSinglesDifference)<-formatC(as.numeric(colnames(numHBVSinglesDifference)), format = "e", digits = 2)
# colnames(numHBVSinglesDifference)[52]<-"Mean"
# cf<-as.data.frame(numHBVSinglesDifference)
# celifloor<-numHBVSinglesDifference
# for(x in 1:length(cf)){
#   for (y in 1:length(cf$Mean)){
#     if(cf[y,x]>0){
#       celifloor[y,x]<-ceiling(cf[y,x])
#     }
#     if(cf[y,x]<0){
#       celifloor[y,x]<-floor(cf[y,x])
#     }
#   }
# }
#
# color_breaks = seq(min(na.omit(as.numeric(unlist(celifloor)))),
#              max(abs(na.omit(as.numeric(unlist(celifloor))))),
#              length.out=9)
# gradient1 = colorpanel( sum( color_breaks[-1]<=0), "deeppink4", "floralwhite" )
# gradient2 = colorpanel( sum( color_breaks[-1]>0 ), "floralwhite", "aquamarine4" )
# gradient2<-gradient2[2:4]
# hm.colors = c(gradient1,gradient2)
# heatmap.2(celifloor,dendrogram='none',#ceiling(numHBVSinglesDifference),dendrogram='none',
#           trace="none",Colv=F,Rowv=F,col=hm.colors,srtCol = 45,
#           #cellnote = celifloor, notecol = 'black',notecex=.7,
#           sepwidth=c(0.01,0.01),
#           sepcolor="white",
#           colsep=1:ncol(numHBVSinglesDifference),
#           rowsep=1:nrow(numHBVSinglesDifference),
#           main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
#           na.color='lightslategrey')
#
# #o....kay? absolute values? #####
#
# AbsDiffMatrix<-abs(numHBVSinglesDifference)
# color_breaks = seq(min(na.omit(as.numeric(unlist(AbsDiffMatrix)))),
#              max(abs(na.omit(as.numeric(unlist(AbsDiffMatrix))))),
#              length.out=101)
# gradient2 = colorpanel( sum( color_breaks[-1]>0 ), "white", "aquamarine4" )
# abs.hm.colors = gradient2
# heatmap.2(AbsDiffMatrix,dendrogram='none',
#           trace="none",Colv=F,Rowv=F,col=abs.hm.colors,srtCol = 45,
#           sepwidth=c(0.01,0.01),
#           cellnote = round(AbsDiffMatrix,1), notecol = 'black',notecex=.7,
#           sepcolor="white",
#           colsep=1:ncol(numHBVSinglesDifference),
#           rowsep=1:nrow(numHBVSinglesDifference),
#           main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
#           na.color='lightslategrey')
#
# color_breaks = seq(min(na.omit(as.numeric(unlist(round(AbsDiffMatrix))))),
#              max(abs(na.omit(as.numeric(unlist(round(AbsDiffMatrix)))))),
#              length.out=101)
# gradient2 = colorpanel( sum( color_breaks[-1]>0 ), "floralwhite", "aquamarine4" )
# abs.hm.colors = gradient2
#
#
# heatmap.2(round(AbsDiffMatrix),dendrogram='none',
#           trace="none",Colv=F,Rowv=F,col=abs.hm.colors,srtCol = 45,
#           sepwidth=c(0.01,0.01),
#           cellnote = round(AbsDiffMatrix), notecol = 'black',notecex=.7,
#           sepcolor="white",
#           colsep=1:ncol(numHBVSinglesDifference),
#           rowsep=1:nrow(numHBVSinglesDifference),
#           main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
#           na.color='lightslategrey')
# #
