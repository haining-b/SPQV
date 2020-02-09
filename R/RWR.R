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


#Singles Heat Map

HBVSinglesDifferenceMatrix<-as.data.frame(matrix(nrow=8,ncol=length(BootStrapConfidenceIntervals)))
colnames(HBVSinglesDifferenceMatrix)<-Relevant_QTL$Length

for (Experiment in 1:length(BootStrapConfidenceIntervals)){
  DifferenceFromHBV=ceiling(BootStrapConfidenceIntervals[1:8,Experiment])-ceiling(EndGameSingles[Experiment,4])
  HBVSinglesDifferenceMatrix[,Experiment]<-DifferenceFromHBV
}

HBVSinglesDifferenceMatrix<-HBVSinglesDifferenceMatrix[,order(Relevant_QTL$Length) ]

numHBVSinglesDifference<-sapply(HBVSinglesDifferenceMatrix,as.numeric)

AbsDiffMatrix<-abs(numHBVSinglesDifference)

SinglesAbsDiffMatrix<-as.data.frame(as.matrix(AbsDiffMatrix))
SinglesAbsDiffMatrix[is.na(SinglesAbsDiffMatrix)]<-10

heatmap.2(as.matrix(SinglesAbsDiffMatrix),dendrogram='none',trace="none",Colv=F,Rowv=F,col=coul,main='withOUT dups',na.color='lightblue')
#w/o the big bois
heatmap.2(as.matrix(SinglesAbsDiffMatrix[,2:17]),dendrogram='none',trace="none",Colv=F,Rowv=F,col=coul,main='withOUT dups',na.color='lightblue')

#Doubles Heat Map
HBVDifferenceMatrixDups<-as.data.frame(matrix(nrow=8,ncol=length(BootStrapConfidenceIntervals)))
colnames(HBVDifferenceMatrixDups)<-Relevant_QTL$Length

for (Experiment in 1:length(BootStrapConfidenceIntervals)){
  DifferenceFromHBV=ceiling(BootStrapConfidenceIntervals[1:8,Experiment])-ceiling(EndGame[Experiment,4])
  HBVDifferenceMatrixDups[,Experiment]<-DifferenceFromHBV
}

HBVDifferenceMatrixDups<-HBVDifferenceMatrixDups[,order(Relevant_QTL$Length)]

numHBVDifferenceDups<-sapply(HBVDifferenceMatrixDups,as.numeric)
AbsDiffMatrixDups<-abs(numHBVDifferenceDups)


heatmap.2(as.matrix(DupsAbsDiffMatrix),dendrogram='none',trace="none",Colv=F,Rowv=F,col=coul,main='with dups',na.color='lightblue')
heatmap.2(as.matrix(DupsAbsDiffMatrix[,2:17]),dendrogram='none',trace="none",Colv=F,Rowv=F,col=coul,main='withOUT dups',na.color='lightblue')


#Singles and Doubles Heat Map

bothDandS<-cbind(SinglesAbsDiffMatrix,DupsAbsDiffMatrix)
heatmap.2(as.matrix(bothDandS)[,c(2:17,23:38)],dendrogram='none',trace="none",Colv=F,Rowv=F,col=coul,main='with dups',na.color='lightblue')


#What if we only look at the ones we know are realistic

heatmap.2(as.matrix(bothDandS[5:8,]),dendrogram='none',trace="none",Colv=F,Rowv=F,col=coul,main='with dups',na.color='lightblue')

#what happens if I average things
for(r in 1:length(SinglesAbsDiffMatrix$`89563`)){
  CurrentRow<-SinglesAbsDiffMatrix[r,]
  CR<-CurrentRow[!is.na(CurrentRow)]
  print(mean(CR))
}
for(r in 1:length(SinglesAbsDiffMatrix$`89563`)){
  CurrentRow<-DupsAbsDiffMatrix[r,]
  CR<-CurrentRow[!is.na(CurrentRow)]
  print(mean(CR))
}

#W/ buckets

Tento100<-c()
for(r in 1:8){
  CurrentRow<-SinglesAbsDiffMatrix[r,2:4]
  CR<-CurrentRow[!is.na(CurrentRow)]
  Tento100<-c(Tento100,(mean(CR)))
}

Tento100D<-c()
for(r in 1:8){
  CurrentRow<-DupsAbsDiffMatrix[r,2:4]
  CR<-CurrentRow[!is.na(CurrentRow)]
  Tento100D<-c(Tento100D,(mean(CR)))
}

MilliontoTenmil<-c()
for(r in 1:8){
  CurrentRow<-SinglesAbsDiffMatrix[r,5:13]
  CR<-CurrentRow[!is.na(CurrentRow)]
  MilliontoTenmil<-c(MilliontoTenmil,(mean(CR)))
}

MilliontoTenmilD<-c()
for(r in 1:8){
  CurrentRow<-DupsAbsDiffMatrix[r,5:13]
  CR<-CurrentRow[!is.na(CurrentRow)]
  MilliontoTenmilD<-c(MilliontoTenmilD,(mean(CR)))
}

most<-c()
for(r in 1:8){
  CurrentRow<-SinglesAbsDiffMatrix[r,14:17]
  CR<-CurrentRow[!is.na(CurrentRow)]
  most<-c(most,(mean(CR)))
}

mostD<-c()
for(r in 1:8){
  CurrentRow<-DupsAbsDiffMatrix[r,14:17]
  CR<-CurrentRow[!is.na(CurrentRow)]
  mostD<-c(mostD,(mean(CR)))
}

Bucketed_SinglesAbsDiffMatrix<-cbind(Tento100,MilliontoTenmil,most)
heatmap.2(Bucketed_SinglesAbsDiffMatrix,dendrogram='none',trace="none",Colv=F,Rowv=F,col=coul,main='with dups',na.color='lightblue')

Bucketed_DupsAbsDiffMatrix<-cbind(Tento100D,MilliontoTenmilD,mostD)
heatmap.2(Bucketed_DupsAbsDiffMatrix,dendrogram='none',trace="none",Colv=F,Rowv=F,col=coul,main='with dups',na.color='lightblue')
fintimes<-cbind(Bucketed_SinglesAbsDiffMatrix,Bucketed_DupsAbsDiffMatrix)
heatmap.2(fintimes,dendrogram='none',trace="none",Colv=F,Rowv=F,col=coul,main='with dups',na.color='lightblue')











fakebelow10K<-abs(sing.Plor_rchrom[,1])
fakeTento100<-rowMeans(abs(sing.Plor_rchrom[,2:4]))
fakeMilliontoTenmil<-rowMeans(abs(sing.Plor_rchrom[,5:13]))
fakemost<-rowMeans(abs(sing.Plor_rchrom[,14:17]))
faketooBig<-rowMeans(abs(sing.Plor_rchrom[,18:21]))

fakeBucketed_sing.Plor<-cbind(fakebelow10K,fakeTento100,fakeMilliontoTenmil,fakemost,faketooBig)
levelplot(t(fakeBucketed_sing.Plor),
          col.regions = colorRampPalette(c('white','grey','black'),1),
          scale=list(x=list(rot=45)))

fakeDbelow10K<-abs(dupPlor_rchrom[,1])
fakeDTento100<-rowMeans(abs(dupPlor_rchrom[,2:4]))
fakeDMilliontoTenmil<-rowMeans(abs(dupPlor_rchrom[,5:13]))
fakeDmost<-rowMeans(abs(dupPlor_rchrom[,14:17]))
fakeDtooBig<-rowMeans(abs(dupPlor_rchrom[,18:21]))
fakeBucketed_dupPlor<-cbind(fakeDbelow10K,fakeDTento100,fakeDMilliontoTenmil,fakeDmost,fakeDtooBig)
levelplot(t(fakeBucketed_dupPlor),
          col.regions = colorRampPalette(c('white','grey','black'),1),
          scale=list(x=list(rot=45)))


Dbelow10K<-abs(dupPlor[,1])
DTento100<-rowMeans(abs(dupPlor[,2:4]))
DMilliontoTenmil<-rowMeans(abs(dupPlor[,5:13]))
Dmost<-rowMeans(abs(dupPlor[,14:17]))
DtooBig<-rowMeans(abs(dupPlor[,18:21]))
Bucketed_dupPlor<-cbind(Dbelow10K,DTento100,DMilliontoTenmil,Dmost,DtooBig)
levelplot(t(Bucketed_dupPlor),
          col.regions = colorRampPalette(c('white','grey','black'),1),
          scale=list(x=list(rot=45)))

Buckets_real_chrom<-cbind(Bucketed_sing.Plor,Bucketed_dupPlor)

Buckets_fake_chrom<-cbind(fakeBucketed_sing.Plor,fakeBucketed_dupPlor)

levelplot(t(Buckets_real_chrom),
          col.regions = colorRampPalette(c('white','grey','black'),1),
          scale=list(x=list(rot=45)))
levelplot(t(Buckets_fake_chrom),
          col.regions = colorRampPalette(c('white','grey','black'),1),
          scale=list(x=list(rot=45)))
levelplot(t(cbind(Buckets_real_chrom,Buckets_fake_chrom)),
          col.regions = colorRampPalette(c('white','grey','black'),1),
          scale=list(x=list(rot=45)))

#randomizing QTL loci
AllTog<-allTOGETHER[c(1,3,5,7,9,11),c(1:21,43:63)]







ACTUAL_CHR_RELE_QTL<-Relevant_QTL
ACTUAL_CHR_RELE_QTL$ChrNum
Relevant_QTL$ChrNum
for(friendo in 1:length(Relevant_QTL$ChrNum)){
  EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))
  while(length(EveryOptionalMarker$id)<1){
    ChosenChromosome<-sample(1:9,1)
    MarkersOnChromosome <-
      eval(as.name(paste0(
        'MarkersOnChromosome', ChosenChromosome
      )))
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    L_MarkerList <-
      MarkersOnChromosome[which(
        MarkersOnChromosome$Base < (upperLimitChosenChromosome - QTLength) &
          MarkersOnChromosome$Base >
          lowerLimitChosenChromosome
      ), ]
    R_MarkerList <-
      MarkersOnChromosome[which(
        MarkersOnChromosome$Base > (lowerLimitChosenChromosome + QTLength) &
          MarkersOnChromosome$Base <
          upperLimitChosenChromosome
      ), ]
    EveryOptionalMarker <- unique(rbind(L_MarkerList, R_MarkerList))
  }
  Relevant_QTL[friendo,1]<-ChosenChromosome
}


#ok let's get real

howmanyMarkers<-172
num1<-56764
numlast<-41169562

avgqtl<-as.integer(median(sort(Relevant_QTL$Length)[1:17]))
n_bins<-ceiling(chromosomeSize[1,3]/10000)

BinnedMaerker<-c()
for (mark in 1:length(MarkersOnChromosome1$Base)){
  BONFIRE<-max(skipforCounting[skipforCounting[,1]<MarkersOnChromosome1[mark,3],1])
  BinnedMaerker<-c(BinnedMaerker,BONFIRE)
}
CBM<-count(BinnedMaerker)
###1. bounceback
# ?? #####
startpoints<-sample(1:numlast,1000,replace=T)
for (p in 1:length(startpoints)){
  if (startpoints[p] >(chromosomeSize[1,2]-avgqtl)){
    startpoints[p]<-chromosomeSize[1,2]-avgqtl
  }
}
ends<-startpoints+avgqtl
endIT<-cbind(startpoints,ends)
forCounting<-as.data.frame(matrix(nrow=n_bins,ncol = 2))
forCounting[,2]<-0
forCounting[,1]<-seq(from=1,to=41169562,by=10000)


for (i in 1:length(endIT[,1])){
  starting<-endIT[i,1]
  ending<-endIT[i,2]
  stpt<-max(forCounting[which(forCounting[,1]<starting),1])
  if(ending>41160001){
    endpt=max(forCounting[,1])
  }else{ endpt<-min(forCounting[which(forCounting[,1]>ending),1])}
  yeses<-which(forCounting[,1]==stpt)
  endes<-which(forCounting[,1]==endpt)
  forCounting[c(yeses:endes),2]<-forCounting[c(yeses:endes),2]+1

}
head(forCounting)
playingGames<-forCounting[,2]/30
head(MarkersOnChromosome1)
ggplot()+geom_line(aes(x=forCounting[,1],forCounting[,2]))+
  geom_histogram(aes(x=Base, y=(..count..)*25),color="cornflowerblue",fill="cornflowerblue",alpha=.5, data=MarkersOnChromosome1,bins=100)+
  theme_classic()+
  xlab("Position")+
  ylab("Coverage")




#2. skip the end #####
skipstarts<-sample(1:(chromosomeSize[1,2]-avgqtl),1000,replace=T)
skipends<-skipstarts+avgqtl
skipendIT<-cbind(skipstarts,skipends)
skipforCounting<-as.data.frame(matrix(nrow=n_bins,ncol = 2))
skipforCounting[,2]<-0
skipforCounting[,1]<-seq(from=1,to=41169562,by=10000)


for (i in 1:length(skipendIT[,1])){
  starting<-skipendIT[i,1]
  ending<-skipendIT[i,2]
  stpt<-max(skipforCounting[which(skipforCounting[,1]<starting),1])
  if(ending>41160001){
    endpt=max(skipforCounting[,1])
  }else{ endpt<-min(skipforCounting[which(skipforCounting[,1]>ending),1])}
  yeses<-which(skipforCounting[,1]==stpt)
  endes<-which(skipforCounting[,1]==endpt)
  skipforCounting[c(yeses:endes),2]<-skipforCounting[c(yeses:endes),2]+1

}

ggplot()+geom_line(aes(x=skipforCounting[,1],skipforCounting[,2]))+
  geom_histogram(aes(x=Base, y=(..count..)*25),color="cornflowerblue",fill="cornflowerblue",alpha=.5, data=MarkersOnChromosome1,bins=100)+
  theme_classic()+
  xlab("Position")+
  ylab("Coverage")




#3. unidir.markers ######

wellp<-MarkersOnChromosome1[which(MarkersOnChromosome1[,3]<=(chromosomeSize[1,3]-avgqtl)),3]
unimarkers<-sample(wellp, 1000,replace=T)
uniends<-unimarkers+avgqtl
uni<-cbind(unimarkers,uniends)
uniforCounting<-as.data.frame(matrix(nrow=n_bins,ncol = 2))
uniforCounting[,2]<-0
uniforCounting[,1]<-seq(from=1,to=41169562,by=10000)
where<-rep(MarkersOnChromosome1$Base,5)


for (i in 1:length(uni[,1])){
  starting<-uni[i,1]
  ending<-uni[i,2]
  stpt<-max(uniforCounting[which(uniforCounting[,1]<uni[i,1]),1])
  endpt<-min(uniforCounting[which(uniforCounting[,1]>uni[i,2]),1])
  yeses<-which(uniforCounting$V1==stpt)
  endes<-which(uniforCounting$V1==endpt)
  uniforCounting[c(yeses:endes),2]<-uniforCounting[c(yeses:endes),2]+1

}

ggplot()+geom_line(aes(x=uniforCounting$V1,uniforCounting$V2))+
  #geom_point(aes(x=CBM$x,y=CBM$freq))+
  geom_histogram(aes(x=Base, y=(..count..)*25),color="cornflowerblue",fill="cornflowerblue",alpha=.5, data=MarkersOnChromosome1,bins=100)+
  theme_classic()+xlab("Position")+
  ylab("Coverage")

#4. bidir markers ######

bimarkers<-sample(MarkersOnChromosome1[,3], 1000,replace=T)

CATS<-c()
bimarkers<-as.integer(bimarkers)
avgqtl<-as.integer(avgqtl)
for(f in 1:length(bimarkers)){
  if(avgqtl<=bimarkers[f] & bimarkers[f] <=(chromosomeSize[1,3]-avgqtl)){
    direc<-sample(c((bimarkers[f] +avgqtl),(bimarkers[f] -avgqtl)),1)


  }
  if(avgqtl>bimarkers[f] ){
    direc<-(bimarkers[f] +avgqtl)

  }
  if((chromosomeSize[1,3]-avgqtl)<bimarkers[f] ){
    direc<-(bimarkers[f]-avgqtl)
  }
  CATS<-c(CATS,direc)
}

biBoys<-cbind(CATS,bimarkers)

biBoys<-t(apply(biBoys, 1, sort))
biBoys<-apply(biBoys,2,sort)
hist(count(biBoys)[,3])
biforCounting<-as.data.frame(matrix(nrow=n_bins,ncol = 2))
biforCounting[,2]<-0
biforCounting[,1]<-seq(from=1,to=41169562,by=10000)
for (i in 1:length(biBoys[,1])){
  starting<-biBoys[i,1]
  ending<-biBoys[i,2]
  stpt<-max(biforCounting[which(biforCounting[,1]<starting),1])
  if(ending>41160001){
    endpt=max(biforCounting[,1])
  }else{ endpt<-min(biforCounting[which(biforCounting[,1]>ending),1])}

  yeses<-which(biforCounting[,1]==stpt)

  endes<-which(biforCounting[,1]==endpt)
  biforCounting[c(yeses:endes),2]<-biforCounting[c(yeses:endes),2]+1

}

ggplot()+geom_line(aes(x=biforCounting$V1,biforCounting$V2))+
  geom_histogram(aes(x=Base, y=(..count..)*25),color="cornflowerblue",fill="cornflowerblue",alpha=.5, data=MarkersOnChromosome1,bins=100)+
  theme_classic()+
  xlab("Position")+
  ylab("Coverage")


#5. wraparound #####
#:(
toWrapM<-sample(MarkersOnChromosome1[,3], 1000,replace=T)
extras<-as.data.frame(matrix(ncol=2))
CATS<-c()
for(f in 1:length(toWrapM)){
  direc<-sample(c((toWrapM[f] +avgqtl),(toWrapM[f] -avgqtl)),1)
  extrabits<-c()
  if(0>direc){
    newdir<-chromosomeSize[1,2]
    newNub<-newdir+direc
    direc=0
    extrabits<- c(newdir,newNub)
  }
  if(direc>chromosomeSize[1,2]){
    newdir<-0
    newNub<-direc-chromosomeSize[1,2]
    direc<-chromosomeSize[1,2]
    extrabits<- c(newdir,newNub)
}
  CATS<-c(CATS,direc)
  if(length(extrabits)>0){
    extras[length(extras$V1)+1,]<-extrabits
  }
}

Wrapped<-cbind(toWrapM,CATS)
colnames(Wrapped)<-c('start','end')
colnames(extras)<-c('start','end')
Wrapp<-rbind(Wrapped,extras)
Wrapp<-na.omit(Wrapp)
Wrap<-t(apply(Wrapp, 1, sort))
wrapforCounting<-as.data.frame(matrix(nrow=n_bins,ncol = 2))
wrapforCounting[,2]<-0
wrapforCounting[,1]<-seq(from=1,to=41169562,by=10000)
for (i in 1:length(Wrap[,1])){
  starting<-Wrap[i,1]
  ending<-Wrap[i,2]
  if(starting<=0){
    stpt=min(wrapforCounting[,1])
  }else{stpt<-max(wrapforCounting[which(wrapforCounting[,1]<starting),1])}
  if(ending>41160001){
    endpt=max(wrapforCounting[,1])
  }else{ endpt<-min(wrapforCounting[which(wrapforCounting[,1]>ending),1])}

  yeses<-which(wrapforCounting[,1]==stpt)

  endes<-which(wrapforCounting[,1]==endpt)
  wrapforCounting[c(yeses:endes),2]<-wrapforCounting[c(yeses:endes),2]+1

}
# plot ####
ggplot()+geom_line(aes(x=wrapforCounting$V1,wrapforCounting$V2))+
  geom_histogram(aes(x=Base, y=(..count..)*25),color="cornflowerblue",fill="cornflowerblue",alpha=.5, data=MarkersOnChromosome1,bins=100)+
  theme_classic()+
  xlab("Position")+
  ylab("Coverage")







to_string <- as_labeller(c(`1` = "Chromosome 1", `2` = "Chromosome 2",`3` = "Chromosome 3", `4` = "Chromosome 4",`5` = "Chromosome 5", `6` = "Chromosome 6",`7` = "Chromosome 7", `8` = "Chromosome 8",`9` = "Chromosome 9"))

ggplot(data=MarkerList)+
  facet_wrap(~Chromosome,ncol=3,scales = 'free_x',labeller=to_string)+
  geom_histogram(aes(x=Base, y=(..count..)),color="#00ae5f",fill="#00ae5f",stat='bin',alpha=.5,bins=100)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Position")+
  ylab("Number of Markers")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
