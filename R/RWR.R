#Different types of bootstrapping for paper

# paperValidator<-Relevant_QTL[,c(1,2,3,5,4,6)]
Relevant_QTL <- read.csv("R/FakeQTL.csv",stringsAsFactors = F)
colnames(paperValidator)[1]<-'chromosome'

BootStrapConfidenceIntervals<-as.data.frame(matrix(nrow=24,ncol=length(Relevant_QTL$ChrNum)))
colnames(BootStrapConfidenceIntervals)<-paperValidator$Length
library(bootstrap)
ds<-c()

target_gene_list_DUPS<-unique(read.csv("IonomicGeneList.csv",stringsAsFactors = F)[,c(1,3,4)])



complex_Fun<-function(x){
  for(i in length(x)){
    present<-x[i]
    ales<-x[-i]
    alesDot<-sum(ales)/(length(x)-1)

    dSubI<-present-alesDot
    ds<-c(ds,dSubI)

  }
  return(ds)
}


# RWR ######

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
    marker_list[chr_indices & marker_list$Base > last_avail_base, "Base"] <- unavail_equiv_base
  }

  marker_list <- marker_list[which(marker_list$Base > 0)]
  return(marker_list)
}

RWR <- function(
  markers_only,  # vs. random placement
  skip_hangovers,  # vs bounceback
  bidirectional,  # vs left-to-right only
  drop_tandem,
  qtl_list,
  gene_list,
  marker_list,
  chromosome_size
) {

  gene_list$Trait <- "trait"
  if (drop_tandem) {
    gene_list <- FilterGeneList("trait", "extension", gene_list, marker_list)
  }

  qtl_cis <- rep(0, nrow(qtl_list))

  for (qtl_i in nrow(qtl_list)) {
    qtl_length <- qtl_list[qtl_i, "Length"]
    obs_genes_found <- qtl_list[qtl_i, "NumGenes"]  # TODO have to calc this?
    exp_genes_found <- rep(0, n_reps)

    if (markers_only) {
      allowed_markers <- getAllowedMarkers(skip_hangovers, qtl_length, marker_list)
      allowed_markers <- cbind(allowed_markers, is_l_to_r=TRUE)

      if (bidirectional) {
        rl_markers <- cbind(marker_list, is_l_to_r=FALSE)

        chr_lengths <- merge(rl_markers, chromosome_size[, c("Number", "Length")],
                            by.x="Chromosome", by.y="Number")
        rl_markers$Base <- chr_lengths$Length - rl_markers$Base + 1  # reverse bases
        rl_markers <- getAllowedMarkers(skip_hangovers, qtl_length, rl_markers)
        rl_markers$Base <- chr_lengths$Length - rl_markers$Base + 1  # reverse them back

        allowed_markers <- rbind(allowed_markers, rl_markers)
      }
    }

    for (rep_i in n_reps) {
      if (markers_only) {
        sampled_marker <- allowed_markers[sample(nrow(allowed_markers, 1)), ]
        sampled_chr <- sampled_marker$Chromosome

        if (sampled_marker$is_r_to_l) {
          qtl_range <- c(sampled_marker$Base, sampled_marker$Base + qtl_length)
        } else {
          qtl_range <- c(sampled_marker$Base - qtl_length, sampled_marker$Base)
        }

      } else {  # !markers_only (all bases)
        chr_and_range <- sampleUniformQTLLoc(skip_hangovers, bidirectional, qtl_length, chromosome_size)
        sampled_chr <- chr_and_range[[1]]
        qtl_range <- chr_and_range[[2]]
      }

      genes_found = nrow(target_gene_list[
          which(target_gene_list$Chromosome == sampled_chr) &
          which(target_gene_list$Base >= qtl_range[1]) &
          which(target_gene_list$Base <= qtl_range[2])
        ,])
      exp_genes_found[rep_i] <- genes_found

    }  # endfor n_reps

    # no touchy the math
    BS_dist <- exp_genes_found
    GHat <- ecdf(BS_dist)
    z0 <- qnorm(GHat(obs_genes_found))

    dSubI <- bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
    accelConstant <- (1/6) * sum(dSubI$jack.values ^ 3) / ((sum(dSubI$jack.values ^ 2)) ^ (2/3))

    bcaNum <- z0 + 1.96
    bcaDen <- 1 - accelConstant*(z0+1.96)
    if (z0 == Inf) {
      toQuantile <- NA
    } else {
      allIn <- z0 + bcaNum / bcaDen
      toQuantile <- pnorm(allIn)
    }

    qtl_cis[qtl_i] <- quantile(BS_dist, toQuantile)
  } # endfor qtl_list
  return(qtl_cis)
}

# ?? #####
ceiling(BootStrapConfidenceIntervals)
colnames(BootStrapConfidenceIntervals)<-Relevant_QTL$Length

BS_333<-read.csv('333 Genes Bootstrap BCa CIs.csv',stringsAsFactors = F)
HBV_333<-read.csv("HBVOutputfor333_200QTL_6_28.csv",stringsAsFactors = F)
PercentDiffMatrix<-as.data.frame(matrix(nrow=8,ncol=length(HBV_333$QTL)))
colnames(PercentDiffMatrix)<-HBV_333$QTL
for(shannon in 1: length(HBV_333$QTL)){
  percentDiff<-BS_333[,shannon]/HBV_333[shannon,6]
  PercentDiffMatrix[,shannon]<-percentDiff
}


library(gplots)
breaks = seq(min(na.omit(as.numeric(unlist(PercentDiffMatrix)))),
             max(abs(na.omit(as.numeric(unlist(PercentDiffMatrix))))),
             length.out=13) 
gradient1 = colorpanel( sum( breaks<=1), rgb(0,146,146,max=255), "ghostwhite" )
gradient2 = colorpanel( sum( breaks>1 ), "ghostwhite", rgb(73,0,146,max=255) )
hm.colors = c(gradient1,gradient2)
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

colnames(numPercentDiffMatrix)[200]<-'Mean'
par(mar=c(10, 8, 8, 3) + 0.1)
heatmap.2(numPercentDiffMatrix,dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=hm.colors,srtCol = 45,
          #cellnote = round(numHBVSinglesDifference), notecol = 'black',notecex=.7,
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator; Percent', 
          denscol='black',
          na.color='lightslategrey')






RoundPercentDiffMatrix<-as.data.frame(matrix(nrow=8,ncol=length(HBV_333$QTL)))
colnames(RoundPercentDiffMatrix)<-HBV_333$QTL

for(shannon in 1: length(HBV_333$QTL)){
  for(pants in 1:8){
    BS<-round(as.numeric(BS_333[pants,shannon]))
    HBV<-round(as.numeric(HBV_333[shannon,6]))
    
    if(BS==1 &HBV==0){
      print("inf!")
      RoundPercentDiffMatrix[pants,shannon]<-NA
    }
    if(BS==0 &HBV==1){
      print("WTF")
      RoundPercentDiffMatrix[pants,shannon]<-NA
    }
    if(BS==0 & HBV==0){
      print("0!")
      RoundPercentDiffMatrix[pants,shannon]<-1
    }
    else{
      RoundPercentDiffMatrix[pants,shannon]<-BS/HBV
    }
  }
  
}
RoundPercentDiffMatrix[ sapply(RoundPercentDiffMatrix,is.infinite)]<-NA
breaks = seq(min(na.omit(as.numeric(unlist(RoundPercentDiffMatrix)))),
             max(abs(na.omit(as.numeric(unlist(RoundPercentDiffMatrix))))),
             length.out=10)
gradient1 = colorpanel( sum( breaks<1), rgb(0,146,146,max=255), "ghostwhite" )
gradient2 = colorpanel( sum( breaks>1 ), "ghostwhite", rgb(73,0,146,max=255) )
hm.colors = c(gradient1,gradient2)

rRowMeanPerc<-c()
for(r in 1:8){
  CurrentRow<-RoundPercentDiffMatrix[r,]
  CR<-CurrentRow[!is.na(CurrentRow)]
  CR<-abs(CR)
  rRowMeanPerc<-c(rRowMeanPerc,(mean(CR)))
}
RoundPercentDiffMatrix<-cbind(RoundPercentDiffMatrix,rRowMeanPerc)
numRoundPercentDiffMatrix<-sapply(RoundPercentDiffMatrix,as.numeric)

heatmap.2 (numRoundPercentDiffMatrix,dendrogram='none',
           trace="none",Colv=F,Rowv=F,col=hm.colors,srtCol = 45,
           #cellnote = round(numHBVSinglesDifference), notecol = 'black',notecex=.7,
           main='Comparison of Bootstrapping methods to the Haining-Blumer Validator', 
           denscol='black',
           na.color='lightslategrey')










#w/ rounded nums

breaks = seq(min(na.omit(as.numeric(unlist(ceiling(HBVSinglesDifferenceMatrix))))),
             max(abs(na.omit(as.numeric(unlist(ceiling(HBVSinglesDifferenceMatrix)))))),
             length.out=9)
gradient1 = colorpanel( sum( breaks[-1]<=0)+1, "deeppink4", "floralwhite" )
gradient2 = colorpanel( sum( breaks[-1]>0 )-1, "floralwhite", "aquamarine4" )
hm.colors = c(gradient1,gradient2)
formatC(numb, format = "e", digits = 2)
colnames(numHBVSinglesDifference)<-formatC(as.numeric(colnames(numHBVSinglesDifference)), format = "e", digits = 2)
colnames(numHBVSinglesDifference)[52]<-"Mean"
cf<-as.data.frame(numHBVSinglesDifference)
celifloor<-numHBVSinglesDifference
for(x in 1:length(cf)){
  for (y in 1:length(cf$Mean)){
    if(cf[y,x]>0){
      celifloor[y,x]<-ceiling(cf[y,x])
    }
    if(cf[y,x]<0){
      celifloor[y,x]<-floor(cf[y,x])
    }
  }
}

breaks = seq(min(na.omit(as.numeric(unlist(celifloor)))),
             max(abs(na.omit(as.numeric(unlist(celifloor))))),
             length.out=9)
gradient1 = colorpanel( sum( breaks[-1]<=0), "deeppink4", "floralwhite" )
gradient2 = colorpanel( sum( breaks[-1]>0 ), "floralwhite", "aquamarine4" )
gradient2<-gradient2[2:4]
hm.colors = c(gradient1,gradient2)
heatmap.2(celifloor,dendrogram='none',#ceiling(numHBVSinglesDifference),dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=hm.colors,srtCol = 45,
          #cellnote = celifloor, notecol = 'black',notecex=.7,
          sepwidth=c(0.01,0.01),
          sepcolor="white",
          colsep=1:ncol(numHBVSinglesDifference),
          rowsep=1:nrow(numHBVSinglesDifference),
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
          na.color='lightslategrey')

#o....kay? absolute values?

AbsDiffMatrix<-abs(numHBVSinglesDifference)
breaks = seq(min(na.omit(as.numeric(unlist(AbsDiffMatrix)))),
             max(abs(na.omit(as.numeric(unlist(AbsDiffMatrix))))),
             length.out=101)
gradient2 = colorpanel( sum( breaks[-1]>0 ), "white", "aquamarine4" )
abs.hm.colors = gradient2
heatmap.2(AbsDiffMatrix,dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=abs.hm.colors,srtCol = 45,
          sepwidth=c(0.01,0.01),
          cellnote = round(AbsDiffMatrix,1), notecol = 'black',notecex=.7,
          sepcolor="white",
          colsep=1:ncol(numHBVSinglesDifference),
          rowsep=1:nrow(numHBVSinglesDifference),
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
          na.color='lightslategrey')

breaks = seq(min(na.omit(as.numeric(unlist(round(AbsDiffMatrix))))),
             max(abs(na.omit(as.numeric(unlist(round(AbsDiffMatrix)))))),
             length.out=101)
gradient2 = colorpanel( sum( breaks[-1]>0 ), "floralwhite", "aquamarine4" )
abs.hm.colors = gradient2


heatmap.2(round(AbsDiffMatrix),dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=abs.hm.colors,srtCol = 45,
          sepwidth=c(0.01,0.01),
          cellnote = round(AbsDiffMatrix), notecol = 'black',notecex=.7,
          sepcolor="white",
          colsep=1:ncol(numHBVSinglesDifference),
          rowsep=1:nrow(numHBVSinglesDifference),
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
          na.color='lightslategrey')
#






faketargets<-WholeGenomeGeneDistribution[sample(1:length(WholeGenomeGeneDistribution$GeneStart),500),]

interpolatedgenes<-as.data.frame(matrix(ncol=4,nrow=0))
colnames(interpolatedgenes)<-colnames(MarkerList)
for(i in 1:9){
  MarkersOnChromosome<-eval(as.name(paste0('MarkersOnChromosome',i)))
  MarkersOnChromosome$type<-'marker'
  targetsonChr<-faketargets[which(faketargets$Chromosome==i),]
  targetsonChr$type<-'gene'
  targetsonChr$GeneID<-'genenomen'
  targetsonChr<-targetsonChr[,c(4,1,2,3)]
  colnames(targetsonChr)<-colnames(MarkersOnChromosome)
  snoop<-rbind(MarkersOnChromosome,targetsonChr)
  snoop<-snoop[order(snoop$Base),]
  interpolatedgenes<-rbind(interpolatedgenes,snoop)
}

todrop<-c()
for (j in 2:length(interpolatedgenes$id)){
  type1<-interpolatedgenes[j-1,4]
  type2<-interpolatedgenes[j,4]
  if(type1 == 'gene' & type2=='gene'){
    todrop<-c(todrop,j)
  }
}

interpolatedgenes<-interpolatedgenes[-todrop,]
length(todrop)
#using 365 genes

head(WholeGenomeGeneDistribution)
withEnds<-read.csv('allSiGeneswithends.csv',stringsAsFactors = F)
GenesForTheInternet<-as.data.frame(matrix(ncol=3))
colnames(GenesForTheInternet)<-colnames(withEnds)
for (i in 1:9){
  GenesOnChr<-withEnds[which(withEnds$Chromosome==i),]
  GenesOnChr<-GenesOnChr[order(GenesOnChr$GeneStart),]
  threesonChr<-threes[which(threes$Chromosome==i),]
  for (gene in threesonChr$Base){
    p<-GenesOnChr[max(which(GenesOnChr$GeneStart<as.numeric(as.character(gene)))),]
    GenesForTheInternet<-rbind(GenesForTheInternet,p)
  }
  
}

write.csv(GenesForTheInternet,'forHBVPaper333UsedGenes.csv',row.names = F)
head(GenesForTheInternet)

listBoi<-apply(GenesForTheInternet,1,paste,":",sep='',collapse="")
whet<-gsub('.{1}$', '', listBoi)
write.csv(whet,file='startandEnds333geneshbv.csv',row.names = F)
apply(df, 1, paste, collapse="")
length(listBoi)
lilListBoi<-c()
for(i in 1:333){
  
  lilListBoi<-c(lilListBoi,substr(listBoi[i],3,nchar(listBoi[i])))
}
head(lilListBoi)
tail(lilListBoi)
threesGeneDeets<-c()
for(g in 1:333){
  threesGeneDeets<-c(threesGeneDeets, geneDeets[grep(lilListBoi[g],geneDeets)])
}
head(threesGeneDeets)
geneNames<-c()
for(u in 1:333){
  entryLvl<-unlist(strsplit(threesGeneDeets[u]," "))
  futuristic<-entryLvl[grep('gene:',entryLvl)]
  futuristic<-substr(futuristic,6,nchar(futuristic))
  geneNames<-c(geneNames,futuristic)
}
#GenesForTheInternet<-GenesForTheInternet[-1,]
WITHNAMES<-cbind(geneNames,GenesForTheInternet)
write.csv(WITHNAMES,'forHBVPaper333UsedGenesWITHNAMES.csv',row.names = F)
