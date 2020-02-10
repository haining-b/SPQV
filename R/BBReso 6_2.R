#BadBootstrapping reorder


#making fake QTL for exprimentation

#The scaffolds [id,Chromosome, Base]
MarkerList<-read.csv("example_data/Sviridis_MarkerList.csv",stringsAsFactors = F)
#chromosome parameters [Numba, Longness, Last Markers, First Markers]
chromosomeSize<-read.csv("example_data/Sviridis_ChromosomeSizes.csv",stringsAsFactors = F)
number_chromosomes<-length(chromosomeSize$Number)

# #### Regularize
marker_list <- dplyr::rename(MarkerList, ID="id")
marker_list$Base<-as.integer(marker_list$Base)

#Splitting Markers



sectioned_markers <- SectionMarkers(marker_list = marker_list, num_chromosomes = number_chromosomes)

FakeLengths<-exp(seq(log(5000),log(55365370-90984),length.out=200))
n_reps=1
FakeMarkersL<-c()
FakeMarkersR<-c()
FakeChromosomes<-c()
L_MarkerList<-data.frame()
R_MarkerList<-data.frame()

for (QTL in 1:length(FakeLengths)){
  QTLength<-FakeLengths[QTL]
  for(BB2Aii in 1:n_reps){
    L_MarkerList$id<-c()
    R_MarkerList$id<-c()
    EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))

    while(length(EveryOptionalMarker$id)<1){
      ChosenChromosome<-sample(1:9,1)
      MarkersOnChromosome <-
        sectioned_markers[[ChosenChromosome]]
      upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
      lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
      L_MarkerList <-
        MarkersOnChromosome[which(
          as.numeric(as.character(MarkersOnChromosome$Base)) <= as.numeric(as.character(upperLimitChosenChromosome - QTLength))+1 &
            MarkersOnChromosome$Base >= lowerLimitChosenChromosome
        ), ]
      R_MarkerList <-
        MarkersOnChromosome[which(
          MarkersOnChromosome$Base >= (lowerLimitChosenChromosome + QTLength)-1 &
            MarkersOnChromosome$Base <=
            upperLimitChosenChromosome
        ), ]
      EveryOptionalMarker <- rbind( R_MarkerList,L_MarkerList)
    }
    if(length(EveryOptionalMarker$Base)==1){
      fakeQTLspot=EveryOptionalMarker$Base
    } else{
      fakeQTLspot<-sample(EveryOptionalMarker$Base,1)}

    L_List<-as.numeric(length(which(L_MarkerList$Base==fakeQTLspot)))
    R_List<-as.numeric(length(which(R_MarkerList$Base==fakeQTLspot)))
    if (L_List>0 & R_List>0){
      print("Freeeeeeee")
      secondEndpoint<-sample(c(fakeQTLspot+QTLength,fakeQTLspot-QTLength),1)}
    if (L_List>0 & R_List==0){
      print("L")
      secondEndpoint<-fakeQTLspot+QTLength}
    if (L_List==0 & R_List>0){
      print('r')
      secondEndpoint<-fakeQTLspot-QTLength}
    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    FakeMarkersL<-c(FakeMarkersL, edges[1])
    FakeMarkersR<-c(FakeMarkersR, edges[2])
    FakeChromosomes<-c(FakeChromosomes, ChosenChromosome)

  }

}
FakeTreatments<-sample(as.character(QTLData$treatment),length(FakeLengths),replace=T)
FakeTraits<-c()
FakeTypes<-c()
FakeExptTypes<-c()

for (i in FakeTreatments){
  if (i=="sparse"|i=='dense'){
    etype= 'raw'
    yr=sample(c('DN13','DN14'),1)
    exptType="DN"
  }
  if(i=='wet'|i=='dry'){
    etype='raw'
    yr='DR14'
    exptType="DR"
  }
  if(i=='diff'){
    etype='diff'
    yr=sample(c('DN13','DN14','DR14'),1)
    exptType=substr(yr,1,2)
  }
  trait<-paste0(sample(element_short_names[2:length(element_short_names)],1),'_',yr)
  FakeTraits<-c(FakeTraits,trait)
  FakeTypes<-c(FakeTypes,etype)
  FakeExptTypes<-c(FakeExptTypes,exptType)

}
Fakes<-as.data.frame(cbind(FakeChromosomes,round(FakeMarkersL),round(FakeMarkersR),
                           FakeTraits,FakeTreatments,FakeExptTypes,FakeTypes))
colnames(Fakes)<-colnames(QTLData)


Fakes<-read.csv("example_data/FakeQTL.csv",stringsAsFactors = F)


#Fake Gene List

WholeGenomeGeneList<-read.csv("example_data/allSiGeneswithends.csv",stringsAsFactors = F)
head(WholeGenomeGeneList)
WGGqL<-(WholeGenomeGeneList$GeneStart+WholeGenomeGeneList$GeneEnd)/2
WGGL<-as.data.frame(cbind(WholeGenomeGeneList$Chromosome,WGGqL))
colnames(WGGL)<-c("Chromosome","GeneMiddle")
WholeGenomeGeneDistribution<-WGGL

FakeGenes<-sample(WholeGenomeGeneDistribution$GeneMiddle,400, replace=F)
FakeG<-WholeGenomeGeneDistribution[which(WholeGenomeGeneDistribution$GeneMiddle%in% FakeGenes),]
head(FakeG)
FakeG$id<-'genenomen'
FakeG<-FakeG[,c(3,1,2)]
colnames(FakeG)<-colnames(MarkerList)
head(MarkerList)
interleaved<-rbind(unique(FakeG),MarkerList)
interleavec1<-interleaved[order(interleaved$Chromosome,interleaved$Base),]
GeneSpots<-grep("genenomen",interleavec1$id)

head(GeneSpots)

for(i in 2:length(GeneSpots)){
  GeneSpots<-grep("genenomen",interleavec1$id)
  if (GeneSpots[i]-1==GeneSpots[i-1]){
    toPop<-sample(c(i,i-1),1)
    print(i)
    interleavec1<-interleavec1[-GeneSpots[toPop],]
  }
}
FakeGenes<-interleavec1[GeneSpots,]
FakeGenes$Base<-ceiling(FakeGenes$Base)
head(FakeGenes)
colnames(FakeGenes)[1]<-"GeneID"
Fakes$trait<-'PC'
system.time()
#write.csv(FakeGenes,file="FakeHBVGenes333.csv",row.names = F)
FakeGenes<-read.csv("R/FakeHBVGenes333.csv",stringsAsFactors = F)
# Fakes<-read.csv()
library(plyr)
starttime<-Sys.time()

FakeGeneOutput2<-HBVFunLoop(AllQTL = hmmmm[24:25,],TraitofInterest =  "PC",
                            number_repetitions =  1000, gene_list = FakeGenes,
                            chromosome_size =  chromosomeSize,MarkerList = MarkerList,
                            WholeGenomeGeneDistribution = WholeGenomeGeneDistribution,GeneCount='y')
endtime<-Sys.time()
system("say Just finished!")
starttime-endtime

nextone<-rbind(nextone,FakeGeneOutput2)


#hmm let's do this right
nextone<-FakeGeneOutput2[1,]

FakesCounted<-GeneCounter(Fakes, FakeGenes, "PC")
for(i in 26:40){
  print(i)
  rcureenet<-c(i*5-4,i*5)
  forRun<-frrn[rcureenet[1]:rcureenet[2],]
  FakeGeneOutput2<-HBVFunLoop(forRun,"PC",1000,FakeGenes,
                              chromosomeSize,MarkerList,
                              WholeGenomeGeneDistribution,GeneCount='n')
  nextone<-rbind(nextone,FakeGeneOutput2)
  write.csv(nextone,file="example_data/HBVOutputforNowON.csv",row.names = F)


}
unique(nextone[-1,])
FakeGOutput<-HBVFunLoop(frrn,"PC",10,FakeGenes,
                        chromosomeSize,MarkerList,
                        WholeGenomeGeneDistribution,GeneCount='n')


frrn<-GeneCounter(Fakes,FakeGenes,'PC')

forRun=Fakes[131:135,]


HBVFunLoop(Fakes,"PC",1000,FakeGenes,
           chromosomeSize,MarkerList,
           WholeGenomeGeneDistribution,GeneCount='y')

ionGenes<-read.csv('example_data/DirectRepeats5_28_19.csv',stringsAsFactors = F)

sadbrianna11<-HBVFunLoop(Saveme[191:200,],'PC',1000,ionGenes,
           chromosomeSize,MarkerList,
           WholeGenomeGeneDistribution,'n')
head(Fakes)

TheGaMe<-rbind(sadbrianna,sadbrianna2,sadbrianna3,
               sadbrianna4,sadbrianna5,sadbrianna6,
               sadbrianna7,sadbrianna8,sadbrianna9,
               sadbrianna10,sadbrianna11)
#write.csv(TheGaMe,file='IonHBVVALUES.csv',row.names = F)
read.csv('example_data/IonHBVVALUES.csv')
Saveme<-GeneCounter(Fakes,ionGenes,'PC')

length(Fakes$chromosome)
# write.csv(Fakes,"R/FakeQTL.csv",row.names = F)
Fakes<-Fakes[201:400,]
paper
colnames(BootStrapConfidenceIntervals)<-paperValidator$Length
head(FakeSings)
paperValidator<-Fakes
paperValidator$Length<-as.numeric(as.character(Fakes$RCI_pos))-as.numeric(as.character(Fakes$LCI_marker))
paperValidator<-paperValidator[,c(1:3,8,4,6)]
paperValidator<-assay1[,c(1,2,3,5,4,6)]
colnames(paperValidator)[1]<-'chromosome'
huh<-merge(paperValidator,frrn,by="Length")
head(huh)
huh<-huh[,c(2,3,4,1,5,11)]
paperValidator<-huh
BootStrapConfidenceIntervals<-as.data.frame(matrix(nrow=8,ncol=length(paperValidator$chromosome)))
target_gene_list<-FakeGenes
target_gene_list<-interpolatedgenes[,1:3]
target_gene_list<-ionGenes
n_reps=1000
#1. 3,1,4,2,7,5,8,6
############ # 1 # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##B. Bounceback,
###i. R-> L  mapping,cut off whatever's beyond the end of the chromosome
####a. No Dups
###(BB1Bi)
n_hit<-target_gene_list
n_start<-MarkerList
n_hit$hit<-0
n_start$hit<-0
d=1
BB1BiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB1Bi in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
    lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
      lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    }
    fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
    endSpot<-fakeQTLspot+QTLength

    if(endSpot>upperLimitChosenChromosome){
      endSpot=upperLimitChosenChromosome
      fakeQTLspot=endSpot-QTLength
    }

    disorderededges<-c(fakeQTLspot,endSpot)
    edges<-disorderededges[order(disorderededges)]
    n_start[which(n_start$Base==fakeQTLspot),4]<- n_start[which(n_start$Base==fakeQTLspot),4]+1

    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(as.numeric(as.character(target_gene_list[poss,2]))==ChosenChromosome){
        locus<-as.numeric(as.character(target_gene_list[poss,3]))
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1
          n_hit[poss,4]<-n_hit[poss,4]+1
        }
      }

    }
    BB1BiBootstrapper[QTL,BB1Bi]<-GenesWithin
  }
  what<-as.data.frame(as.matrix(bcaboot::bcajack(unlist(BB1BiBootstrapper[QTL,]),1000,mean,verbose=F)$lims))
  BootStrapConfidenceIntervals[d,QTL]<-what[8,1]
}

ggplot(n_hit,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[1])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#####
############ # 2 # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##A. No bounceback,
###i. R-> L mapping only
####a. No Dups
###(BB1Ai)

n_hit$hit<-0
n_start$hit<-0
d=2
BB1AiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB1Ai in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    endSpot=10000000000000000
    upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
    lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
      lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    }
    while(endSpot>upperLimitChosenChromosome){

      fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
      endSpot<-fakeQTLspot+QTLength
    }

    edges<-c(fakeQTLspot,endSpot)
    n_start[which(n_start$Base==fakeQTLspot),4]<- n_start[which(n_start$Base==fakeQTLspot),4]+1
    GenesWithin=0
    #for(poss in 1:length(known_gene_list$GeneID)){
    #for DUPLICATES
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-as.numeric(as.character(target_gene_list[poss,3]))
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1
          n_hit[poss,4]<-n_hit[poss,4]+1
        }
      }

    }
    BB1AiBootstrapper[QTL,BB1Ai]<-GenesWithin
  }

  what<-as.data.frame(as.matrix(bcaboot::bcajack(unlist(BB1AiBootstrapper[QTL,]),1000,mean,verbose=F)$lims))
  BootStrapConfidenceIntervals[d,QTL]<-what[8,1]


}
ggplot(n_hit,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[2])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#####
############ # 3 # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##B. Bounceback,
###ii. R-> L, L->R
####a. No Dups
###(BB1Bii)
d=3
n_hit$hit<-0
n_start$hit<-0

paperValidator<-paperValidator[-200,]

BB1BiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  for(BB1Bii in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
    lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
      lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    }
    fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
    endSpot<-sample(c(fakeQTLspot+QTLength,fakeQTLspot-QTLength),1)

    if(endSpot>upperLimitChosenChromosome){
      endSpot=upperLimitChosenChromosome
      fakeQTLspot=endSpot-QTLength
    }
    if(endSpot<lowerLimitChosenChromosome){
      endSpot=lowerLimitChosenChromosome
      fakeQTLspot=endSpot+QTLength
    }


    disorderededges<-c(fakeQTLspot,endSpot)
    edges<-disorderededges[order(disorderededges)]
    n_start[which(n_start$Base==fakeQTLspot),4]<- n_start[which(n_start$Base==fakeQTLspot),4]+1

    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-as.numeric(as.character(target_gene_list[poss,3]))
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1
          n_hit[poss,4]<-n_hit[poss,4]+1
        }
      }

    }
    BB1BiiBootstrapper[QTL,BB1Bii]<-GenesWithin
  }
  what<-as.data.frame(as.matrix(bcaboot::bcajack(unlist(BB1BiiBootstrapper[QTL,]),1000,mean,verbose=F)$lims))
  BootStrapConfidenceIntervals[d,QTL]<-what[8,1]
}

ggplot(n_hit,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[3])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#####
############ # 4 # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##A. No bounceback,
###ii. R-> L OR L-> R mapping,cut off whatever's beyond the end of the chromosome
####a. No Dups
###(BB1Aii)
n_hit$hit<-0
n_start$hit<-0
d=4
BB1AiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]

  for(BB1Aii in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    endSpot=10000000000000000
    upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
    lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
      lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    }
    while(endSpot>upperLimitChosenChromosome| endSpot<lowerLimitChosenChromosome){
      fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
      endSpot<-sample(c(fakeQTLspot+QTLength,fakeQTLspot-QTLength),1)
    }

    disorderededges<-c(fakeQTLspot,endSpot)
    edges<-disorderededges[order(disorderededges)]
    n_start[which(n_start$Base==fakeQTLspot),4]<- n_start[which(n_start$Base==fakeQTLspot),4]+1

    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-as.numeric(as.character(target_gene_list[poss,3]))
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1
          n_hit[poss,4]<-n_hit[poss,4]+1
        }
      }

    }
    BB1AiiBootstrapper[QTL,BB1Aii]<-GenesWithin
  }
  what<-as.data.frame(as.matrix(bcaboot::bcajack(unlist(BB1AiiBootstrapper[QTL,]),1000,mean,verbose=F)$lims))
  BootStrapConfidenceIntervals[d,QTL]<-what[8,1]
}
ggplot(n_hit,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[4])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#####
############ # 5 # ###############
#Bad Bootstrapping 2: Marker Placement
##B.Bounceback,
###i. L->R mapping
####a. No Dups
###(BB2Bi)
d=5
n_hit$hit<-0
n_start$hit<-0
BB2BiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  print(QTL)
  obs_value= paperValidator[QTL,6]
  for(BB2Bi in 1:n_reps){

    EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))


    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
    lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
      lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    }


    MarkersOnChromosome <-
      eval(as.name(paste0(
        'MarkersOnChromosome', ChosenChromosome
      )))

    fakeQTLspot<-sample(MarkersOnChromosome$Base,1)
    secondEndpoint<-fakeQTLspot+QTLength
    if(secondEndpoint>=upperLimitChosenChromosome){
      fakeQTLspot<-max(MarkersOnChromosome$Base)
      secondEndpoint<-fakeQTLspot-QTLength
    }

    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    n_start[which(n_start$Base==fakeQTLspot),4]<- n_start[which(n_start$Base==fakeQTLspot),4]+1

    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-as.numeric(as.character(target_gene_list[poss,3]))
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1
          n_hit[poss,4]<-n_hit[poss,4]+1

        }
      }

    }
    BB2BiBootstrapper[QTL,BB2Bi]<-GenesWithin
  }
  what<-as.data.frame(as.matrix(bcaboot::bcajack(unlist(BB2BiBootstrapper[QTL,]),1000,mean,verbose=F)$lims))
  BootStrapConfidenceIntervals[d,QTL]<-what[8,1]
}
ggplot(n_hit,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[5])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(n_start,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[5])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#####
############ # 6 # ###############
#Bad Bootstrapping 2: Marker Placement
##A. No bounceback,
###i. R-> L mapping only , cut off whatever's beyond the end of the chromosome
####a. No ds
###(BB2Ai)
d=6

n_hit$hit<-0
n_start$hit<-0
BB2AiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  print(QTL)
  QTLength<-paperValidator[QTL,4]
  for(BB2Ai in 1:n_reps){
    L_MarkerList$id<-c()
    while(length(L_MarkerList$id)<1){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
      lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
      while(upperLimitChosenChromosome-lowerLimitChosenChromosome< QTLength){
        ChosenChromosome<-sample(1:9,1)
        upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
        lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
      }
      EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))
      MarkersOnChromosome <-
        eval(as.name(paste0(
          'MarkersOnChromosome', ChosenChromosome
        )))

      L_MarkerList <-
        MarkersOnChromosome[which(
          MarkersOnChromosome$Base <= as.numeric(as.character(chromosomeSize[ChosenChromosome, 4])) - QTLength+1 &
            MarkersOnChromosome$Base >= lowerLimitChosenChromosome), ]
    }

    EveryOptionalMarker<-L_MarkerList

    fakeQTLspot<-sample(EveryOptionalMarker$Base,1)

    secondEndpoint<-fakeQTLspot+QTLength
    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    n_start[which(n_start$Base==fakeQTLspot),4]<- n_start[which(n_start$Base==fakeQTLspot),4]+1
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-as.numeric(as.character(target_gene_list[poss,3]))
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1
          n_hit[poss,4]<-n_hit[poss,4]+1
        }
      }

    }
    BB2AiBootstrapper[QTL,BB2Ai]<-GenesWithin
  }
  what<-as.data.frame(as.matrix(bcaboot::bcajack(unlist(BB2AiBootstrapper[QTL,]),1000,mean,verbose=F)$lims))
  BootStrapConfidenceIntervals[d,QTL]<-what[8,1]
}

ggplot(n_hit,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[6])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(n_start,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[6])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#####
############ # 7 # ###############
#Bad Bootstrapping 2: Marker Placement
##B.Bounceback,
###ii. bidirectional mapping
####a. No Dups
###(BB2Bii)
d=7
n_start$hit<-0
n_hit$hit<-0
BB2BiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  print(QTL)
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Bii in 1:n_reps){

    EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))


    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
    lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
      lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
    }


    MarkersOnChromosome <-
      eval(as.name(paste0(
        'MarkersOnChromosome', ChosenChromosome
      )))

    fakeQTLspot<-sample(MarkersOnChromosome$Base,1)
    secondEndpoint<-sample(c(fakeQTLspot+QTLength,fakeQTLspot-QTLength),1)
    if(secondEndpoint>upperLimitChosenChromosome){
      fakeQTLspot<-max(MarkersOnChromosome$Base)
      secondEndpoint<-fakeQTLspot-QTLength
    }
    if(secondEndpoint<lowerLimitChosenChromosome){
      fakeQTLspot<-min(MarkersOnChromosome$Base)
      secondEndpoint<-fakeQTLspot+QTLength
    }

    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    n_start[which(n_start$Base==fakeQTLspot),4]<- n_start[which(n_start$Base==fakeQTLspot),4]+1

    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-as.numeric(as.character(target_gene_list[poss,3]))
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1
          n_hit[poss,4]<-n_hit[poss,4]+1
        }
      }

    }
    BB2BiiBootstrapper[QTL,BB2Bii]<-GenesWithin
  }

  what<-as.data.frame(as.matrix(bcaboot::bcajack(unlist(BB2BiiBootstrapper[QTL,]),1000,mean,verbose=F)$lims))
  BootStrapConfidenceIntervals[d,QTL]<-what[8,1]
}
ggplot(n_hit,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[7])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(n_start,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[7])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#####
############ # 8 # ###############
#Bad Bootstrapping 2: Marker Placement
##A. No bounceback,
###i. R-> L, L->R mapping
####a. No Dups
###(BB2Aii)

d=8
BB2AiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))

n_hit$hit<-0
n_start$hit<-0
for (QTL in 1:length(paperValidator$chromosome)){
  print(QTL)
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Aii in 1:n_reps){
    L_MarkerList$id<-c()
    R_MarkerList$id<-c()
    EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))

    while(length(EveryOptionalMarker$id)<1){
      ChosenChromosome<-sample(1:9,1)
      MarkersOnChromosome <-
        eval(as.name(paste0(
          'MarkersOnChromosome', ChosenChromosome
        )))
      upperLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))
      lowerLimitChosenChromosome <- as.numeric(as.character(chromosomeSize[ChosenChromosome, 3]))
      L_MarkerList <-
        MarkersOnChromosome[which(
          MarkersOnChromosome$Base <= (as.numeric(as.character(chromosomeSize[ChosenChromosome, 4]))- QTLength+1) &
            MarkersOnChromosome$Base >
            lowerLimitChosenChromosome
        ), ]
      R_MarkerList <-
        MarkersOnChromosome[which(
          MarkersOnChromosome$Base >= (lowerLimitChosenChromosome + QTLength-1) &
            MarkersOnChromosome$Base <=
            upperLimitChosenChromosome
        ), ]
      EveryOptionalMarker <- rbind( R_MarkerList,L_MarkerList)
    }
    if(length(EveryOptionalMarker$Base)==1){
      fakeQTLspot=EveryOptionalMarker$Base
    } else{
      fakeQTLspot<-sample(EveryOptionalMarker$Base,1)}

    L_List<-as.numeric(length(which(L_MarkerList$Base==fakeQTLspot)))
    R_List<-as.numeric(length(which(R_MarkerList$Base==fakeQTLspot)))
    if (L_List>0 & R_List>0){
      secondEndpoint<-sample(c(fakeQTLspot+QTLength,fakeQTLspot-QTLength),1)}
    if (L_List>0 & R_List==0){
      secondEndpoint<-fakeQTLspot+QTLength}
    if (L_List==0 & R_List>0){
      secondEndpoint<-fakeQTLspot-QTLength}
    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    n_start[which(n_start$Base==fakeQTLspot),4]<- n_start[which(n_start$Base==fakeQTLspot),4]+1
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-as.numeric(as.character(target_gene_list[poss,3]))
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1
          n_hit[poss,4]<-n_hit[poss,4]+1
        }
      }

    }
    BB2AiiBootstrapper[QTL,BB2Aii]<-GenesWithin
  }
  what<-as.data.frame(as.matrix(bcaboot::bcajack(unlist(BB2AiiBootstrapper[QTL,]),1000,mean,verbose=F)$lims))
  BootStrapConfidenceIntervals[d,QTL]<-what[8,1]
}

ggplot(n_hit,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[8])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(n_start,aes(x=Base,y=hit))+
  geom_point(col=rainbowCols[8])+
  facet_rep_wrap(~Chromosome,ncol=3,scales='free_y', repeat.tick.labels = 'left',labeller=to_string)+
  theme_minimal()+
  theme(axis.line = element_line())+
  scale_x_continuous()+
  xlab("Position")+
  ylab("Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#####

# Plot 1? #######

#for 71 genes, 51 QTL
#don't be an idiot don't run this code jesus
#seventyoneBS<-BootStrapConfidenceIntervals
#for 365Genes, 51 QTL
#threesixtyfiveBS<-BootStrapConfidenceIntervals
#write.csv(BootStrapConfidenceIntervals,file = "71 Genes Bootstrap BCa CIs.csv",row.names = F)

head(nextone)
forPlay<-nextone[order(nextone$QTL),]
forPlay<-TheGaMe[order(TheGaMe$QTL),]
nextone<-nextone[-1,]
HBVSinglesDifferenceMatrix<-as.data.frame(matrix(nrow=8,ncol=length(BootStrapConfidenceIntervals)))
colnames(HBVSinglesDifferenceMatrix)<-Fakes$Length
for (Experiment in 1:length(BootStrapConfidenceIntervals)){
  DifferenceFromHBV=BootStrapConfidenceIntervals[,Experiment]-forPlay[Experiment,6]
  HBVSinglesDifferenceMatrix[,Experiment]<-DifferenceFromHBV
}


#Color palette
library(gplots)
breaks = seq(min(na.omit(as.numeric(unlist(HBVSinglesDifferenceMatrix[,c(115:199,201)])))),
             max(abs(na.omit(as.numeric(unlist(HBVSinglesDifferenceMatrix[,c(115:199,201)]))))),
             length.out=51)
gradient1 = colorpanel( sum( breaks[-1]<=0), rgb(0,146,146,max=255), "ghostwhite" )
gradient2 = colorpanel( sum( breaks[-1]>0 ), "ghostwhite", rgb(73,0,146,max=255) )
hm.colors = c(gradient1,gradient2)

#Getting row means
RowMeanHBV<-c()
for(r in 1:8){
  CurrentRow<-HBVSinglesDifferenceMatrix[r,115:199]
  CR<-CurrentRow[!is.na(CurrentRow)]
  CR<-abs(CR)
  RowMeanHBV<-c(RowMeanHBV,(mean(CR)))
}
HBVSinglesDifferenceMatrix<-cbind(HBVSinglesDifferenceMatrix,RowMeanHBV)
numHBVSinglesDifference<-sapply(HBVSinglesDifferenceMatrix,as.numeric)
colnames(numHBVSinglesDifference)[52]<-'Mean'
 colnames(numHBVSinglesDifference)[202]<-'Mean'
 par(mar=c(10, 8, 8, 3) + 0.1)
heatmap.2(numHBVSinglesDifference[,c(115:199,202)],dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=hm.colors,srtCol = 45,
          #cellnote = round(numHBVSinglesDifference), notecol = 'black',notecex=.7,
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
          denscol='black',
          na.color='lightslategrey')

#Percent Diff - best ##############

BS_333<-read.csv('R/333 Genes Bootstrap BCa CIs.csv',stringsAsFactors = F)
HBV_333<-read.csv("R/HBVOutputfor333_200QTL_6_28.csv",stringsAsFactors = F)

# qtl_range_to_show <- 1:199
# BS_333 <- BSCI[1:8, qtl_range_to_show+1]
# HBV_333 <- SPQV_results[qtl_range_to_show, ]

PercentDiffMatrix<-as.data.frame(matrix(nrow=8,ncol=nrow(HBV_333)))
colnames(PercentDiffMatrix)<-HBV_333$QTL
for(shannon in 1: length(HBV_333$QTL)){
  percentDiff<-BS_333[,shannon]/HBV_333[shannon,6] # "Upper 95% CI"
  PercentDiffMatrix[,shannon]<-percentDiff
}


HBVSinglesDifferenceMatrix<-as.data.frame(matrix(nrow=8,ncol=ncol(BS_333)))
# colnames(HBVSinglesDifferenceMatrix)<-Fakes$Length
for (exp_i in 1:ncol(BS_333)){
  DifferenceFromHBV=BS_333[,exp_i]-HBV_333[exp_i,6]
  HBVSinglesDifferenceMatrix[,exp_i]<-DifferenceFromHBV
}
HBVSinglesDifferenceMatrix <- HBVSinglesDifferenceMatrix[, 1:ncol(HBVSinglesDifferenceMatrix)-1]


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
old_hmap <- heatmap.2(numPercentDiffMatrix,dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=hm.colors,srtCol = 45,
          #cellnote = round(numHBVSinglesDifference), notecol = 'black',notecex=.7,
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator; Percent',
          denscol='black',
          na.color='lightslategrey')




# RoundPercentDiff?? - gives almost all 0 #############

RoundPercentDiffMatrix<-as.data.frame(matrix(nrow=8,ncol=length(HBV_333$QTL)))
colnames(RoundPercentDiffMatrix)<-HBV_333$QTL

for(shannon in 1: length(HBV_333$QTL)){
  for(pants in 1:8){
    BS<-round(as.numeric(BS_333[pants,shannon]))
    HBV<-round(as.numeric(HBV_333[shannon,6]))

    # if (is)

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

heatmap.2(numRoundPercentDiffMatrix,dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=hm.colors,srtCol = 45,
          #cellnote = round(numHBVSinglesDifference), notecol = 'black',notecex=.7,
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
          denscol='black',
          na.color='lightslategrey')










#w/ rounded nums ##############

breaks = seq(min(na.omit(as.numeric(unlist(ceiling(HBVSinglesDifferenceMatrix))))),
             max(abs(na.omit(as.numeric(unlist(ceiling(HBVSinglesDifferenceMatrix)))))),
             length.out=9)
gradient1 = colorpanel( sum( breaks[-1]<=0)+1, "deeppink4", "floralwhite" )
gradient2 = colorpanel( sum( breaks[-1]>0 )-1, "floralwhite", "aquamarine4" )
hm.colors = c(gradient1,gradient2)
# formatC(numb, format = "e", digits = 2)
numHBVSinglesDifference <- HBVSinglesDifferenceMatrix
# colnames(numHBVSinglesDifference)<-formatC(as.numeric(colnames(numHBVSinglesDifference)), format = "e", digits = 2)
# colnames(numHBVSinglesDifference)[52]<-"Mean"
cf<-as.data.frame(numHBVSinglesDifference)
celifloor<-numHBVSinglesDifference
for(x in 1:ncol(cf)){
  for (y in 1:nrow(cf)){
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
heatmap.2(x=as.matrix(celifloor),dendrogram='none',#ceiling(numHBVSinglesDifference),dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=hm.colors,srtCol = 45,
          #cellnote = celifloor, notecol = 'black',notecex=.7,
          sepwidth=c(0.01,0.01),
          sepcolor="white",
          colsep=1:ncol(numHBVSinglesDifference),
          rowsep=1:nrow(numHBVSinglesDifference),
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
          na.color='lightslategrey')

#o....kay? absolute values? ##############

AbsDiffMatrix<-abs(numHBVSinglesDifference)
breaks = seq(min(na.omit(as.numeric(unlist(AbsDiffMatrix)))),
             max(abs(na.omit(as.numeric(unlist(AbsDiffMatrix))))),
             length.out=101)
gradient2 = colorpanel( sum( breaks[-1]>0 ), "white", "aquamarine4" )
abs.hm.colors = gradient2
heatmap.2(as.matrix(AbsDiffMatrix),dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=abs.hm.colors,srtCol = 45,
          sepwidth=c(0.01,0.01),
          # cellnote = round(AbsDiffMatrix,1), notecol = 'black',notecex=.7,
          # sepcolor="white",
          colsep=1:ncol(numHBVSinglesDifference),
          rowsep=1:nrow(numHBVSinglesDifference),
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
          na.color='lightslategrey')

##?? abs val? ############
breaks = seq(min(na.omit(as.numeric(unlist(round(AbsDiffMatrix))))),
             max(abs(na.omit(as.numeric(unlist(round(AbsDiffMatrix)))))),
             length.out=101)
gradient2 = colorpanel( sum( breaks[-1]>0 ), "floralwhite", "aquamarine4" )
abs.hm.colors = gradient2

heatmap.2(as.matrix(round(AbsDiffMatrix)),dendrogram='none',
          trace="none",Colv=F,Rowv=F,col=abs.hm.colors,srtCol = 45,
          sepwidth=c(0.01,0.01),
          cellnote = round(AbsDiffMatrix), notecol = 'black',notecex=.7,
          sepcolor="white",
          colsep=1:ncol(numHBVSinglesDifference),
          rowsep=1:nrow(numHBVSinglesDifference),
          main='Comparison of Bootstrapping methods to the Haining-Blumer Validator',
          na.color='lightslategrey')







#Getting rid of dups?? ###########

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
withEnds<-read.csv('example_data/allSiGeneswithends.csv',stringsAsFactors = F)
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
