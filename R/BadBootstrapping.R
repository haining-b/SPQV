#Different types of bootstrapping for paper

paperValidator<-Relevant_QTL[,c(1,2,3,5,4,6)]
colnames(paperValidator)[1]<-'chromosome'

BootStrapConfidenceIntervals<-as.data.frame(matrix(nrow=24,ncol=length(Relevant_QTL$ChrNum)))
colnames(BootStrapConfidenceIntervals)<-paperValidator$Length
library(bootstrap)
ds<-c()
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
target_gene_list_DUPS<-unique(read.csv("IonomicGeneList.csv",stringsAsFactors = F)[,c(1,3,4)])


############ # 1  | rand    | no_bb | RL+cutoff   | no_dup | BB1Ai  # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##A. No bounceback,
###i. R-> L mapping only,cut off whatever's beyond the end of the chromosome
####a. No Dups
###(BB1Ai)

#paperValidator<-toValidate[c(4,5,14),]


d=1
BB1AiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB1Ai in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    endSpot=10000000000000000
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    }
    while(endSpot>upperLimitChosenChromosome){

      fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
      endSpot<-fakeQTLspot+QTLength
    }

    edges<-c(fakeQTLspot,endSpot)
    GenesWithin=0
    #for(poss in 1:length(known_gene_list$GeneID)){
    #for DUPLICATES
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB1AiBootstrapper[QTL,BB1Ai]<-GenesWithin
  }

  BS_dist<-unlist(BB1AiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught <-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum <-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)


}


############ # 2  | rand    | no_bb | both+cutoff | no_dup | BB1Aii # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##A. No bounceback,
###ii. R-> L OR L-> R mapping,cut off whatever's beyond the end of the chromosome
####a. No Dups
###(BB1Aii)

d=d+1
BB1AiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]

  for(BB1Aii in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    endSpot=10000000000000000
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    }
    while(endSpot>upperLimitChosenChromosome| endSpot<lowerLimitChosenChromosome){
      fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
      endSpot<-sample(c(fakeQTLspot+QTLength,fakeQTLspot-QTLength),1)
    }

    disorderededges<-c(fakeQTLspot,endSpot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB1AiiBootstrapper[QTL,BB1Aii]<-GenesWithin
  }
  BS_dist<-unlist(BB1AiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}


############ # 3  | rand    | bb    | RL+cutoff   | no_dup | BB1Bi  # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##B. Bounceback,
###i. R-> L  mapping,cut off whatever's beyond the end of the chromosome
####a. No Dups
###(BB1Bi)
d=d+1
BB1BiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB1Bi in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    }
    fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
    endSpot<-fakeQTLspot+QTLength

    if(endSpot>upperLimitChosenChromosome){
      endSpot=upperLimitChosenChromosome
      fakeQTLspot=endSpot-QTLength
    }

    disorderededges<-c(fakeQTLspot,endSpot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB1BiBootstrapper[QTL,BB1Bi]<-GenesWithin
  }
  BS_dist<-unlist(BB1BiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}

############ # 4  | rand    | bb    | both        | no_dup | BB1Bii # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##B. Bounceback,
###ii. R-> L, L->R
####a. No Dups
###(BB1Bii)
d=d+1
BB1BiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB1Bii in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
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
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB1BiiBootstrapper[QTL,BB1Bii]<-GenesWithin
  }
  BS_dist<-unlist(BB1BiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}

############ # 5  | mrk     | no_bb | RL+cutoff   | no_dup | BB2Ai  # ###############
#Bad Bootstrapping 2: Marker Placement
##A. No bounceback,
###i. R-> L mapping only , cut off whatever's beyond the end of the chromosome
####a. No Dups
###(BB2Ai)
d=d+1
BB2AiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Ai in 1:n_reps){
    L_MarkerList$id<-c()
    while(length(L_MarkerList$id)<1){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
      while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    }
      EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))
      MarkersOnChromosome <-
        eval(as.name(paste0(
          'MarkersOnChromosome', ChosenChromosome
        )))

      L_MarkerList <-
        MarkersOnChromosome[which(
          MarkersOnChromosome$Base < (upperLimitChosenChromosome - QTLength) &
            MarkersOnChromosome$Base >
            lowerLimitChosenChromosome), ]
    }

    EveryOptionalMarker<-L_MarkerList

    fakeQTLspot<-sample(EveryOptionalMarker$Base,1)

    secondEndpoint<-fakeQTLspot+QTLength
    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB2AiBootstrapper[QTL,BB2Ai]<-GenesWithin
  }
  BS_dist<-unlist(BB2AiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}



############ # 6  | mrk     | no_bb | both        | no_dup | BB2Aii # ###############
#Bad Bootstrapping 2: Marker Placement
##A. No bounceback,
###i. R-> L, L->R mapping
####a. No Dups
###(BB2Aii)
d=d+1
BB2AiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Aii in 1:n_reps){
    L_MarkerList$id<-c()
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
    if(length(EveryOptionalMarker$Base)==1){
      fakeQTLspot=EveryOptionalMarker$Base}
    else{
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
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB2AiiBootstrapper[QTL,BB2Aii]<-GenesWithin
  }
  BS_dist<-unlist(BB2AiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}



############ # 7  | mrk     | bb    | LR          | no_dup | BB2Bi  # ###############
#Bad Bootstrapping 2: Marker Placement
##B.Bounceback,
###i. L->R mapping
####a. No Dups
###(BB2Bi)
d=d+1
BB2BiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Bi in 1:n_reps){

    EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))


    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
      while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
        ChosenChromosome<-sample(1:9,1)
        upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
        lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
      }


      MarkersOnChromosome <-
        eval(as.name(paste0(
          'MarkersOnChromosome', ChosenChromosome
        )))

      fakeQTLspot<-sample(MarkersOnChromosome$Base,1)
      secondEndpoint<-fakeQTLspot+QTLength
    if(secondEndpoint>upperLimitChosenChromosome){
      fakeQTLspot<-max(MarkersOnChromosome$Base)
      secondEndpoint<-fakeQTLspot-QTLength
    }

    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB2BiBootstrapper[QTL,BB2Bi]<-GenesWithin
  }
  BS_dist<-unlist(BB2BiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}



############ # 8  | mrk     | bb    | bidi        | no_dup | BB2Bii # ###############
#Bad Bootstrapping 2: Marker Placement
##B.Bounceback,
###ii. bidirectional mapping
####a. No Dups
###(BB2Bii)
d=d+1
BB2BiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Bii in 1:n_reps){

    EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))


    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
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
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB2BiiBootstrapper[QTL,BB2Bii]<-GenesWithin
  }
  BS_dist<-unlist(BB2BiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}



############ # 9  | qtl_mrk | no_bb | LR          | no_dup | BB3Ai  # ###############

#Bad Bootstrapping 3: QTL specific Marker Placement
##A.No bounceback,
###i. L->R mapping
####a. No Dups
###(BB3Ai)
d=d+1
BB3AiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  ChosenChromosome<-paperValidator[QTL,1]
  EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))
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

  EveryOptionalMarker<-L_MarkerList
  for(catastrophe in 1:n_reps){

    fakeQTLspot<-sample(EveryOptionalMarker$Base,1)
    secondEndpoint<-fakeQTLspot+QTLength

    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){

        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB3AiBootstrapper[QTL,catastrophe]<-GenesWithin
  }
  BS_dist<-unlist(BB3AiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)

}


############ # 10 | qtl_mrk | no_bb | both        | no_dup | BB3Aii # ###############
#Bad Bootstrapping 3: QTL specific Marker Placement
##A.No bounceback,
###ii. L->R, R->L mapping
####a. No Dups
###(BB3Aii)
d=d+1
BB3AiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  ChosenChromosome<-paperValidator[QTL,1]
  EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))

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

  for(catastrophe in 1:n_reps){

    if(length(EveryOptionalMarker$Base)==1){
      fakeQTLspot=EveryOptionalMarker$Base
    } else{
      fakeQTLspot<-sample(EveryOptionalMarker$Base,1)
    }
    L_List<-as.numeric(length(which(L_MarkerList$Base==fakeQTLspot)))
    R_List<-as.numeric(length(which(R_MarkerList$Base==fakeQTLspot)))
    if (L_List>0 & R_List>0){
      secondEndpoint<-sample(c(fakeQTLspot+QTLength,fakeQTLspot-QTLength),1)
    }
    if (L_List>0 & R_List==0){
      secondEndpoint<-fakeQTLspot+QTLength
    }
    if (L_List==0 & R_List>0){
      secondEndpoint<-fakeQTLspot-QTLength
    }
    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){

        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB3AiiBootstrapper[QTL,catastrophe]<-GenesWithin
  }
  BS_dist<-unlist(BB3AiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}


############ # 11 | qtl_mrk | bb    | LR          | no_dup | BB3Bi  # ###############
#Bad Bootstrapping 3: QTL specific Marker Placement
##B.Bounceback,
###i. L->R mapping
####a. No Dups
###(BB3Bi)
d=d+1
BB3BiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  ChosenChromosome<-as.numeric(paperValidator[QTL,1])
  obs_value= paperValidator[QTL,6]

  MarkersOnChromosome <-
    eval(as.name(paste0(
      'MarkersOnChromosome', ChosenChromosome
    )))
  upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
  lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]


  for(catastrophe in 1:n_reps){
    fakeQTLspot<-sample(MarkersOnChromosome$Base,1)
    secondEndpoint<-fakeQTLspot+QTLength
    if(secondEndpoint>upperLimitChosenChromosome){
      fakeQTLspot<-max(MarkersOnChromosome$Base)
      secondEndpoint<-fakeQTLspot-QTLength
    }
    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){

        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB3BiBootstrapper[QTL,catastrophe]<-GenesWithin
  }
  BS_dist<-unlist(BB3BiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}


############ # 12 | qtl_mrk | bb    | both        | no_dup | BB3Bii # ###############
#Bad Bootstrapping 3: QTL specific Marker Placement
##B.Bounceback
###ii. L->R, R->L mapping
####a. No Dups
###(BB3Bii)
d=d+1
BB3BiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  ChosenChromosome<-paperValidator[QTL,1]


  MarkersOnChromosome <-
    eval(as.name(paste0(
      'MarkersOnChromosome', ChosenChromosome
    )))
  upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
  lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]


  for(catastrophe in 1:n_reps){
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
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){

        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB3BiiBootstrapper[QTL,catastrophe]<-GenesWithin
  }
  BS_dist<-unlist(BB3BiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}







#########
target_gene_list<-target_gene_list_DUPS
############ # 13 | rand    | no_bb | RL+cutoff   | no_dup | BB1Ai  (assumed bc identical to #1) # ###############

d=13
BB1AiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]

  for(BB1Ai in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    endSpot=10000000000000000
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    }
    while(endSpot>upperLimitChosenChromosome){

      fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
      endSpot<-fakeQTLspot+QTLength
    }

    edges<-c(fakeQTLspot,endSpot)
    GenesWithin=0
    #for(poss in 1:length(known_gene_list$GeneID)){
    #for DUPLICATES
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB1AiBootstrapper[QTL,BB1Ai]<-GenesWithin
  }
  BS_dist<-unlist(BB1AiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}


############ # 14 | rand    | no_bb | both+cutoff | no_dup | BB1Aii # ###############

#Bad Bootstrapping 1: Random Placement (no markers involved)
##A. No bounceback,
###ii. R-> L OR L-> R mapping,cut off whatever's beyond the end of the chromosome
####a. No Dups
###(BB1Aii)

d=d+1
BB1AiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]

  for(BB1Aii in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    endSpot=10000000000000000
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    }
    while(endSpot>upperLimitChosenChromosome| endSpot<lowerLimitChosenChromosome){
      fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
      endSpot<-sample(c(fakeQTLspot+QTLength,fakeQTLspot-QTLength),1)
    }

    disorderededges<-c(fakeQTLspot,endSpot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB1AiiBootstrapper[QTL,BB1Aii]<-GenesWithin
  }
  BS_dist<-unlist(BB1AiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}


############ # 15 | rand    | bb    | RL+cutoff   | no_dup | BB1Bi  # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##B. Bounceback,
###i. R-> L  mapping,cut off whatever's beyond the end of the chromosome
####a. No Dups
###(BB1Bi)
d=d+1
BB1BiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB1Bi in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    }
    fakeQTLspot<-sample(lowerLimitChosenChromosome:upperLimitChosenChromosome,1)
    endSpot<-fakeQTLspot+QTLength

    if(endSpot>upperLimitChosenChromosome){
      endSpot=upperLimitChosenChromosome
      fakeQTLspot=endSpot-QTLength
    }

    disorderededges<-c(fakeQTLspot,endSpot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB1BiBootstrapper[QTL,BB1Bi]<-GenesWithin
  }
  BS_dist<-unlist(BB1BiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}


############ # 16 | rand    | bb    | both        | no_dup | BB1Bii # ###############
#Bad Bootstrapping 1: Random Placement (no markers involved)
##B. Bounceback,
###ii. R-> L, L->R
####a. No Dups
###(BB1Bii)
d=d+1
BB1BiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB1Bii in 1:n_reps ){
    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
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
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB1BiiBootstrapper[QTL,BB1Bii]<-GenesWithin
  }
  BS_dist<-unlist(BB1BiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}



############ # 17 | mrk     | no_bb | RL+cutoff   | no_dup | BB2Ai  # ###############
#Bad Bootstrapping 2: Marker Placement
##A. No bounceback,
###i. R-> L mapping only , cut off whatever's beyond the end of the chromosome
####a. No Dups
###(BB2Ai)
d=d+1
BB2AiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Ai in 1:n_reps){
    L_MarkerList$id<-c()
    while(length(L_MarkerList$id)<1){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
      while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
        ChosenChromosome<-sample(1:9,1)
        upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
        lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
      }
      EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))
      MarkersOnChromosome <-
        eval(as.name(paste0(
          'MarkersOnChromosome', ChosenChromosome
        )))

      L_MarkerList <-
        MarkersOnChromosome[which(
          MarkersOnChromosome$Base < (upperLimitChosenChromosome - QTLength) &
            MarkersOnChromosome$Base >
            lowerLimitChosenChromosome
        ), ]
    }

    EveryOptionalMarker<-L_MarkerList

    fakeQTLspot<-sample(EveryOptionalMarker$Base,1)

    secondEndpoint<-fakeQTLspot+QTLength
    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB2AiBootstrapper[QTL,BB2Ai]<-GenesWithin
  }
  BS_dist<-unlist(BB2AiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}



############ # 18 | mrk     | no_bb | both+cutoff | no_dup | BB2Aii # ###############
#Bad Bootstrapping 2: Marker Placement
##A. No bounceback,
###i. R-> L, L->R mapping
####a. No Dups
###(BB2Aii)
d=d+1
BB2AiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Aii in 1:n_reps){
    L_MarkerList$id<-c()
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
    if(length(EveryOptionalMarker$Base)==1){
      fakeQTLspot=EveryOptionalMarker$Base}
    else{
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
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB2AiiBootstrapper[QTL,BB2Aii]<-GenesWithin
  }
  BS_dist<-unlist(BB2AiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}


############ # 19 | mrk     | bb    | LR          | no_dup | BB2Bi  # ###############

#Bad Bootstrapping 2: Marker Placement
##B.Bounceback,
###i. L->R mapping
####a. No Dups
###(BB2Bi)
d=d+1
BB2BiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Bi in 1:n_reps){

    EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))


    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    }


    MarkersOnChromosome <-
      eval(as.name(paste0(
        'MarkersOnChromosome', ChosenChromosome
      )))

    fakeQTLspot<-sample(MarkersOnChromosome$Base,1)
    secondEndpoint<-fakeQTLspot+QTLength
    if(secondEndpoint>upperLimitChosenChromosome){
      fakeQTLspot<-max(MarkersOnChromosome$Base)
      secondEndpoint<-fakeQTLspot-QTLength
    }

    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB2BiBootstrapper[QTL,BB2Bi]<-GenesWithin
  }
  BS_dist<-unlist(BB2BiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}


############ # 20 | mrk     | bb    | bidi        | no_dup | BB2Bii # ###############
#Bad Bootstrapping 2: Marker Placement
##B.Bounceback,
###ii. bidirectional mapping
####a. No Dups
###(BB2Bii)
d=d+1
BB2BiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  for(BB2Bii in 1:n_reps){

    EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))


    ChosenChromosome<-sample(1:9,1)
    upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
    lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
    while(upperLimitChosenChromosome-lowerLimitChosenChromosome<QTLength){
      ChosenChromosome<-sample(1:9,1)
      upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
      lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]
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
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){
        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB2BiiBootstrapper[QTL,BB2Bii]<-GenesWithin
  }

  BS_dist<-unlist(BB2BiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
  }



####
############ # 21 | qtl_mrk | no_bb | LR          | no_dup | BB3Ai  # ###############

#Bad Bootstrapping 3: QTL specific Marker Placement
##A.No bounceback,
###i. L->R mapping
####a. No Dups
###(BB3Ai)
d=d+1
BB3AiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  ChosenChromosome<-paperValidator[QTL,1]
  EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))
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

  EveryOptionalMarker<-L_MarkerList
  for(catastrophe in 1:n_reps){

    fakeQTLspot<-sample(EveryOptionalMarker$Base,1)
    secondEndpoint<-fakeQTLspot+QTLength

    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){

        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB3AiBootstrapper[QTL,catastrophe]<-GenesWithin
  }
  BS_dist<-unlist(BB3AiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}

############ # 22 | qtl_mrk | no_bb | both        | no_dup | BB3Aii # ###############

#Bad Bootstrapping 3: QTL specific Marker Placement
##A.No bounceback,
###ii. L->R, R->L mapping
####a. No Dups
###(BB3Aii)
d=d+1
BB3AiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  obs_value= paperValidator[QTL,6]
  ChosenChromosome<-paperValidator[QTL,1]
  EveryOptionalMarker<-as.data.frame(matrix(nrow=0,ncol=3))

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

  for(catastrophe in 1:n_reps){

    if(length(EveryOptionalMarker$Base)==1){
      fakeQTLspot=EveryOptionalMarker$Base
    } else{
      fakeQTLspot<-sample(EveryOptionalMarker$Base,1)
    }
    L_List<-as.numeric(length(which(L_MarkerList$Base==fakeQTLspot)))
    R_List<-as.numeric(length(which(R_MarkerList$Base==fakeQTLspot)))
    if (L_List>0 & R_List>0){
      secondEndpoint<-sample(c(fakeQTLspot+QTLength,fakeQTLspot-QTLength),1)
    }
    if (L_List>0 & R_List==0){
      secondEndpoint<-fakeQTLspot+QTLength
    }
    if (L_List==0 & R_List>0){
      secondEndpoint<-fakeQTLspot-QTLength
    }
    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){

        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB3AiiBootstrapper[QTL,catastrophe]<-GenesWithin
  }
  BS_dist<-unlist(BB3AiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}



############ # 23 | qtl_mrk | bb    | LR          | no_dup | BB3bi  # ###############
#Bad Bootstrapping 3: QTL specific Marker Placement
##B.Bounceback,
###i. L->R mapping
####a. No Dups
###(BB3Bi)
d=d+1
BB3BiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  ChosenChromosome<-paperValidator[QTL,1]
  obs_value= paperValidator[QTL,6]

  MarkersOnChromosome <-
    eval(as.name(paste0(
      'MarkersOnChromosome', ChosenChromosome
    )))
  upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
  lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]


  for(catastrophe in 1:n_reps){
    fakeQTLspot<-sample(MarkersOnChromosome$Base,1)
    secondEndpoint<-fakeQTLspot+QTLength
    if(secondEndpoint>upperLimitChosenChromosome){
      fakeQTLspot<-max(MarkersOnChromosome$Base)
      secondEndpoint<-fakeQTLspot-QTLength
    }
    disorderededges<-c(secondEndpoint,fakeQTLspot)
    edges<-disorderededges[order(disorderededges)]
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){

        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB3BiBootstrapper[QTL,catastrophe]<-GenesWithin
  }
  BS_dist<-unlist(BB3BiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)

}


############ # 24 | qtl_mrk | bb    | both        | no_dup | BB3Bii # ###############
#Bad Bootstrapping 3: QTL specific Marker Placement
##B.Bounceback
###ii. L->R, R->L mapping
####a. No Dups
###(BB3Bii)
d=d+1
BB3BiiBootstrapper<-as.data.frame(matrix(data=0,nrow=length(paperValidator$chromosome),ncol=n_reps))
for (QTL in 1:length(paperValidator$chromosome)){
  QTLength<-paperValidator[QTL,4]
  ChosenChromosome<-paperValidator[QTL,1]
  obs_value= paperValidator[QTL,6]

  MarkersOnChromosome <-
    eval(as.name(paste0(
      'MarkersOnChromosome', ChosenChromosome
    )))
  upperLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 4]
  lowerLimitChosenChromosome <- chromosomeSize[ChosenChromosome, 3]


  for(catastrophe in 1:n_reps){
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
    GenesWithin=0
    for(poss in 1:length(target_gene_list$GeneID)){
      if(target_gene_list[poss,2]==ChosenChromosome){

        locus<-target_gene_list[poss,3]
        if(locus>=edges[1]& locus<=edges[2]){
          GenesWithin=GenesWithin+1

        }
      }

    }
    BB3BiiBootstrapper[QTL,catastrophe]<-GenesWithin
  }
  BS_dist<-unlist(BB3BiiBootstrapper[QTL,])
  GHat<-ecdf(BS_dist)
  zNaught<-qnorm(GHat(obs_value))

  dSubI<-bootstrap::jackknife(as.vector(BS_dist),complex_Fun)
  accelConstant<-(1/6)*sum(dSubI$jack.values^3)/((sum(dSubI$jack.values^2))^(2/3))

  bcaNum<-zNaught+1.96
  bcaDen<-1- accelConstant*(zNaught+1.96)
  if (zNaught==Inf){
    toQuantile<-NA
  }else{
    allIn<-zNaught+bcaNum/bcaDen
    toQuantile<-pnorm(allIn)
  }

  BootStrapConfidenceIntervals[d,QTL]<-quantile(BS_dist,toQuantile)
}

# ?? #####
ceiling(BootStrapConfidenceIntervals)
colnames(BootStrapConfidenceIntervals)<-Relevant_QTL$Length

RWR <- function(
  markers_only, # bool (vs. random placement)
  no_bounceback, # bool
  cutoff_beyond_end, # bool (was this a paste error - isnt it bounceback OR cutoff OR skip?)
  bidirectional, # bool
  drop_duplicates # bool
) {
  temp -> "TEMP"
}

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
