
#' Create chromosome specific matrices of QTL mapping marker locations
#'
#' This function produces N matrices containing the locations of the
#' markers used in a QTL mapping experiment, where N refers to the number
#' of chromosomes in the organism of interest.
#'
#' @param number_chromosomes Numeric value refering the the number
#' of chromosomes in the organism of interest
#'
#' @param MarkerList A 3 column matrix containing the marker loci
#' (ID, Chromosome, and Base)).
#'
#' @return N matrices of marker loci
#' @return Sectioned_List, a list of the names of the N matrices of marker loci
#' @export
Marker_Sectioner<-function(number_chromosomes,MarkerList){
  colnames(MarkerList)<-c("ID","Chromosome","Base")
  Sectioned_List<-c()
  for(chromosome in 1:number_chromosomes){
    assign(paste("MarkersOnChromosome", chromosome, sep = ""), MarkerList[which(MarkerList$Chromosome==chromosome),],envir = .GlobalEnv)
    Sectioned_List<-c(Sectioned_List,paste("MarkersOnChromosome", chromosome, sep = ""))
  }
  return(Sectioned_List)

}


#' Create known gene lists appropriate for specific traits and
#' your style of QTL mapping
#'
#' This function takes a gene list and refines it based on the specific trait you
#' would like to examine. Will remove duplicate rows. Additionally, the function
#' will act to reduce the final list in the case of QTL mapping that involved extension of the 95% confidence
#' intervals to the nearest marker, so as to avoid the artificial inflation of
#' identified genes that linkage causes. This reduction is performed such that a single
#' random gene from any set of genes that lie between two consecutive markers remains.
#'
#' @param MarkerList A 3 column matrix containing the marker loci
#' (ID, Chromosome, and Base)).
#' @param gene_list A 4 column matrix containing a list of known genes.
#' Columns should include "GeneID', "Trait", "Chromosome", and "Start_Site".
#' @param Placement_Type Either 'extension' or 'centered'. This depends on the
#' mapping style of the original QTL experiment.
#' @param Trait The trait for which you would like to produce a gene list.
#'
#' @return A 4 column matrix containing "GeneID', "Trait", "Chromosome",
#' and "Start_Site" for each remaining gene in the trait specified.
#' @export

Gene_List_Groomer<-function(MarkerList,gene_list,Placement_Type,Trait){

  temp_gene_list<-unique(gene_list)

  colnames(temp_gene_list)<-c('GeneID',"Trait","Chromosome","Start_Site")

  temp_gene_list<-temp_gene_list[grep(Trait,temp_gene_list$Trait),]


  if(length(temp_gene_list$Trait)==0){
    return("No genes for this trait.")
  }
  if(length(temp_gene_list$Trait)==1){
    return(temp_gene_list)
  }
  if(Placement_Type=='centered'){
    Gene_List<-temp_gene_list
  }
  if(Placement_Type=='extension'){

    colnames(MarkerList)<-c('GeneID',"Chromosome","Start_Site")
    MarkerList$Trait<-"Marker"
    LociList<-rbind(MarkerList,temp_gene_list)
    LociList<-LociList[order(LociList$Chromosome,LociList$Start_Site),]
    GeneSpots<-grep(Trait,LociList$Trait)
    Tandem_Arrays<-c()
    for (i in 2:length(GeneSpots)){
      if(GeneSpots[i]==GeneSpots[i-1]+1){
        Tandem_Arrays[i]<-1
        Tandem_Arrays[i-1]<-1
      }else{
        Tandem_Arrays[i]<-0
      }
    }
    if(is.na(Tandem_Arrays[1])){
      Tandem_Arrays[1]<-0
    }
    SequentialIncreases<-rle(Tandem_Arrays)
    vals<-SequentialIncreases$values
    longs<-SequentialIncreases$lengths
    presenceAbsence<-c()
    for(seqi in 1:length(vals)){
      if(vals[seqi]==0){
        presenceAbsence<-c(presenceAbsence,rep(1,longs[seqi]))
      }
      else{
        posSelected<-sample(1:longs[seqi],1)
        lowSubstituteZeroes<-rep(0,(posSelected-1))
        highSubstituteZeroes<-rep(0,longs[seqi]-posSelected)
        AllSubs<-c(lowSubstituteZeroes,1,highSubstituteZeroes)
        presenceAbsence<-c(presenceAbsence,AllSubs)
      }
    }

    Gene_List<-unique(temp_gene_list[grep(1,presenceAbsence),])

  }
  return(Gene_List)
}

#' Determine the probability of placing a particular QTL on any marker in the genome
#'
#' This function determines the likelihood of placing a QTL on every individual
#' marker in the genome. It takes into account the sizes of your identified QTL,
#' the lengths of each chromosome, the locations of markers, and the type of mapping
#' that originally produced the QTL.
#'
#' @param QTLofInterest A matrix containing (at minimum) the chromosome on which the QTL
#' was originally found ('Chromosome'), the Leftmost Confidence interval's location,
#' and the Rightmost Conficence interval's location.
#' @param chromosome_size A matrix containing the number of each chromosome,
#'  their total lengths, and the first and last markers used for the original
#'  mapping experiment.
#' @param Sectioned_List The list of the names of the N matrices of marker loci
#'  resulting from use of the Marker_Sectioner function.
#' @param Placement_Type Either 'extension' or 'centered'. This depends on the
#' mapping style of the original QTL experiment.
#' @param MarkerList A 3 column matrix containing the marker loci
#' (ID, Chromosome, and Base)).
#' @return A list of the probability of placing each QTL on any individual marker
#' in the genome
#' @export
QTL_Placement_Probabilities<-function(QTLofInterest,chromosome_size,Sectioned_List,Placement_Type){
  number_chromosomes<-as.numeric(length(chromosome_size$First.Markers))
  colnames(QTLofInterest)[1:3]<-c('Chromosome',"Leftmost_Marker",
                             "Rightmost_Marker")
  if(length(QTLofInterest$Length)==0){
    QTLofInterest$Length<-as.numeric(QTLofInterest$Rightmost_Marker)-as.numeric(QTLofInterest$Leftmost_Marker)

  }
  PossiblePositions<-as.data.frame(matrix(nrow=length(QTLofInterest$Chromosome),ncol=1+number_chromosomes))
  colnames(PossiblePositions)<-c("Length_QTL",1:number_chromosomes)
  PossiblePositions$Length_QTL<-QTLofInterest$Length

  if(Placement_Type=='centered'){
    for (QTL in 1:length(QTLofInterest$Length)){
      half_QTL_Length<-as.numeric(as.character(QTLofInterest$Length[QTL]))/2

      for (Chromosome in  1:number_chromosomes){
        FirstMarker<-as.numeric(as.character(chromosome_size$First.Markers[Chromosome]))
        LastMarker<-as.numeric(as.character(chromosome_size$Last.Markers[Chromosome]))

        firstAvailableMarker<- as.numeric(FirstMarker+ half_QTL_Length)
        lastAvailableMarker <-as.numeric(LastMarker - half_QTL_Length)

        Relevant_Markers<-as.data.frame(as.matrix(eval(as.name(Sectioned_List[Chromosome]))))


        LR_Remaining_Markers<-Relevant_Markers[which(as.numeric(as.character(Relevant_Markers$Base))<lastAvailableMarker),]
        RL_Remaining_Markers<-Relevant_Markers[which(as.numeric(as.character(Relevant_Markers$Base))>firstAvailableMarker),]

        Remaining_Markers<-rbind(LR_Remaining_Markers,RL_Remaining_Markers)
        Length_Remaining_Markers<-length(Remaining_Markers$ID)

        PossiblePositions[QTL,Chromosome+1]<-Length_Remaining_Markers

        Probabilities<-1/rowSums(PossiblePositions[,2:10])
      }
    }
  }
  if(Placement_Type=='extension'){
    for (QTL in 1:length(QTLofInterest$Length)){
      QTL_Length<-as.numeric(as.character(QTLofInterest$Length[QTL]))
      for (Chromosome in  1:number_chromosomes){
        FirstMarker<-as.numeric(as.character(chromosome_size$First.Markers[Chromosome]))
        LastMarker<-as.numeric(as.character(chromosome_size$Last.Markers[Chromosome]))

        firstAvailableMarker<- as.numeric(FirstMarker+ QTL_Length)
        lastAvailableMarker <-as.numeric(LastMarker - QTL_Length)

        Relevant_Markers<-as.data.frame(as.matrix(eval(as.name(Sectioned_List[Chromosome]))))


        LR_Remaining_Markers<-Relevant_Markers[which(as.numeric(as.character(Relevant_Markers$Base))<lastAvailableMarker),]
        RL_Remaining_Markers<-Relevant_Markers[which(as.numeric(as.character(Relevant_Markers$Base))>firstAvailableMarker),]

        Remaining_Markers<-rbind(LR_Remaining_Markers,RL_Remaining_Markers)
        Length_Remaining_Markers<-length(Remaining_Markers$ID)

        PossiblePositions[QTL,Chromosome+1]<-Length_Remaining_Markers

        Probabilities<-1/rowSums(PossiblePositions[,2:10])
      }
    }
  }
  return(Probabilities)
}




#' Count the number of known genes under a particular QTL
#'
#' This function determines the number of genes found by a QTL in the original
#' mapping experiment. Relies on the Gene_List_Groomer function to ensure
#' the appropriate genes are being used.
#'
#' @param QTL_with_Metadata A 7 column matrix containing the chromosome
#' on which the QTL was originally found ('Chromosome'), the Leftmost
#'  Confidence interval's location,the Rightmost Conficence interval's location,
#'  the trait for which the QTL were discovered, the treatment that was applied,
#'  the mapping method, the overall type of experiment, and the lengths of the QTL.
#'
#' @param gene_list A 4 column matrix containing a list of known genes.
#' Columns should include "GeneID', "Trait", "Chromosome", and "Start_Site".
#' @param Placement_Type Either 'extension' or 'centered'. This depends on the
#' mapping style of the original QTL experiment.
#' @param Trait The trait for which you would like to produce a gene list.
#' @param MarkerList A 3 column matrix containing the marker loci
#' (ID, Chromosome, and Base)).
#'
#' @return A matrix containing the number of genes identified within a QTL
#' for a specific trait, as well as the number of QTL that were identified
#' for that trait
#' @export
GeneCounter<-function(QTL_with_Metadata,gene_list,Trait,Placement_Type, MarkerList){

  colnames(QTL_with_Metadata)<-c('Chromosome',"Leftmost_Marker",
                       "Rightmost_Marker","Trait","Treatment","QTL_Type")
  QTL_with_Metadata$N_Genes<-c(0)
  temp_gene_list<-Gene_List_Groomer(MarkerList,gene_list,Placement_Type,Trait)
  if(length(temp_gene_list)< 2 ){

    return(temp_gene_list)

  }


  if(as.numeric(as.character(length(QTL_with_Metadata[,1])))>1){
    TraitsandTreatments<-subset(QTL_with_Metadata, select= c(Trait, Treatment))
    sepTreatments_QTL_Count<-plyr::ddply(TraitsandTreatments,.(TraitsandTreatments$Trait,TraitsandTreatments$Treatment),nrow)
    colnames(sepTreatments_QTL_Count)<-c("Trait",'Treatment','Frequency')



    QTLData<-QTL_with_Metadata[grep(Trait, TraitsandTreatments$Trait),]
    if(length(QTLData$Chromosome)==0){
      return("No QTL found for this trait.")
    }

    identified_genes<-data.frame(matrix(ncol=length(QTLData)))
    colnames(identified_genes)<-colnames(QTLData)


    for (i in 1:length(QTLData$Chromosome)){
      QTLchromosome<-as.numeric(as.character(QTLData$Chromosome[i]))
      QTLLCI<-as.numeric(as.character(QTLData$Leftmost_Marker[i]))
      QTLRCI<-as.numeric(as.character(QTLData$Rightmost_Marker[i]))

      for (Gene in 1: length(temp_gene_list$GeneID)){
        if(as.numeric(as.character(temp_gene_list$Chromosome[Gene]))==as.numeric(as.character(QTLchromosome))){
          GeneStartSite<-as.numeric(temp_gene_list$Start_Site[Gene])
          if (GeneStartSite>=QTLLCI & GeneStartSite<=QTLRCI){
              QTLData$N_Genes[i]<-paste0(QTLData$N_Genes[i]," and ",temp_gene_list$GeneID[Gene])
              identified_genes<-rbind(identified_genes,QTLData[i,])

          }
        }
      }
    }
    identified_genes<-identified_genes[-1,]
    identified_genes$Length<-as.numeric(as.character(identified_genes$Rightmost_Marker))-as.numeric(as.character(identified_genes$Leftmost_Marker))


    identified_genes$N_Genes<-stringr::str_count(identified_genes$N_Genes,'and')
    identified_genes2<-identified_genes
    if(length(identified_genes[,1])>1){
      for(possible_Duplicates in 1:(length(identified_genes[,1])-1)){
        if(all(identified_genes[possible_Duplicates,c(1:6)]==identified_genes[possible_Duplicates+1,c(1:6)])){
          identified_genes2[possible_Duplicates,c(1:7)]<-0
        }
      }
      if(length(which(identified_genes2$N_Genes==0))>0){
        CountedIdentifiedGenes<-identified_genes2[-which(identified_genes2$N_Genes==0),]
      } else{
        CountedIdentifiedGenes<-identified_genes2
      }
      CountedIdentifiedGenes$N_QTL<-c(0)
      count<-c()

      for(identified in 1:length(CountedIdentifiedGenes$Chromosome)){
        sep_Treat_QTL_4_u<-sepTreatments_QTL_Count$Frequency[which(sepTreatments_QTL_Count$Trait==CountedIdentifiedGenes$Trait[identified] &
                                                           sepTreatments_QTL_Count$Treatment==CountedIdentifiedGenes$Treatment[identified])]
        count<-c(count,sep_Treat_QTL_4_u)
      }

      CountedIdentifiedGenes$N_QTL<-count

      CountedIdentifiedGenesOutput<-CountedIdentifiedGenes
      toComp<-QTL_with_Metadata
      colnames(toComp)<-colnames(CountedIdentifiedGenes)[1:7]
      toComp$Chromosome<-as.double(toComp$Chromosome)
      toComp$Leftmost_Marker<-as.double(toComp$Leftmost_Marker)
      toComp$Rightmost_Marker<-as.double(toComp$Rightmost_Marker)
      toComp$Treatment<-as.character(toComp$Treatment)


      toComp$N_Genes<-0
      NoGenes<-setdiff(QTLData[1:7],CountedIdentifiedGenes[1:7])
      if(length(NoGenes$Chromosome)>0){
        NoGenes<-NoGenes[which(NoGenes$N_Genes==0),]
        NoGenes$Length<-as.numeric(as.character(NoGenes$Rightmost_Marker))-as.numeric(as.character(NoGenes$Leftmost_Marker))
        NoGenes$Last<-"NA"
        NoGenesOutput<-NoGenes
        colnames(NoGenesOutput)<-c("Chromosome", "Leftmost_Marker",
                                   "Rightmost_Marker", "Trait",
                                   "Treatment","QTL_Type","N_Genes",
                                   "Length", "Number_Trait_QTL")
        colnames(CountedIdentifiedGenesOutput)<-c("Chromosome", "Leftmost_Marker",
                                                  "Rightmost_Marker", "Trait", "Treatment", "QTL_Type",
                                                  "N_Genes", "Length", "Number_Trait_QTL")


        CountedIdentifiedGenesOutput2<-rbind(CountedIdentifiedGenesOutput,NoGenesOutput)
      }else{
        CountedIdentifiedGenesOutput2<-CountedIdentifiedGenesOutput
      }

      return(CountedIdentifiedGenesOutput2)
    }
    if(length(identified_genes[,1])==1){
      count<-c()

      for(identified in 1:length(identified_genes$Chromosome)){
        sep_Treat_QTL_4_u<-sepTreatments_QTL_Count[which(sepTreatments_QTL_Count$Trait==identified_genes$Trait[identified] &
                                                           sepTreatments_QTL_Count$Treatment==identified_genes$Treatment[identified]),3]
        count<-c(count,sep_Treat_QTL_4_u)
      }

      identified_genes$QTL_Number<-count
      forReturn<-identified_genes
      colnames(forReturn)<-c("Chromosome", "Leftmost_Marker",
                             "Rightmost_Marker", "Trait", "Treatment", "QTL_Type",
                             "N_Genes", "Length","Number_Trait_QTL")
      return(forReturn)
    } else{
      return("No identified genes for this trait")
    }
  }
  else{
    i=1
    QTL_with_Metadata$Length<-as.numeric(as.character(QTL_with_Metadata$Rightmost_Marker))-as.numeric(as.character(QTL_with_Metadata$Leftmost_Marker))
    QTL_with_Metadata$N_Genes<-as.character(QTL_with_Metadata$N_Genes)
    QTLchromosome<-as.numeric(as.character(QTL_with_Metadata$Chromosome[i]))
    QTLLCI<-as.numeric(as.character(QTL_with_Metadata$Leftmost_Marker[i]))
    QTLRCI<-as.numeric(as.character(QTL_with_Metadata$Rightmost_Marker[i]))

    identified_genes<-data.frame(matrix(ncol=length(QTL_with_Metadata)))
    colnames(identified_genes)<-colnames(QTL_with_Metadata)

    for (Gene in 1: length(gene_list$GeneID)){
      if(as.numeric(as.character(temp_gene_list$Chromosome[Gene]))==as.numeric(as.character(QTLchromosome))){
        GeneStartSite<-as.numeric(temp_gene_list$Start_Site[Gene])
        if (GeneStartSite>=QTLLCI & GeneStartSite<=QTLRCI){
          QTL_with_Metadata$N_Genes[i]<-paste0(QTL_with_Metadata$N_Genes[i]," and ",temp_gene_list$GeneID[Gene])
          identified_genes<-rbind(identified_genes,QTL_with_Metadata[i,])
        }
      }
    }
    identified_genes<-identified_genes[-1,]
    identified_genes$N_Genes<-stringr::str_count(identified_genes$N_Genes,'and')
    identified_genes2<-identified_genes
    if(length(as.numeric(as.character(identified_genes[,1])))>1){
      for(possible_Duplicates in 1:(length(identified_genes[,1])-1)){
        if(all(identified_genes[possible_Duplicates,c(1:6)]==identified_genes[possible_Duplicates+1,c(1:6)])){
          identified_genes2[possible_Duplicates,c(1:7)]<-0
        }
      }

    }
    if(length(which(identified_genes2$N_Genes==0))>0){
      CountedIdentifiedGenes<-identified_genes2[-which(identified_genes2$N_Genes==0),]
    } else{
      CountedIdentifiedGenes<-identified_genes2
    }
    CountedIdentifiedGenes<-CountedIdentifiedGenes[,c(1:4,8,7)]

    CountedIdentifiedGenes$N_QTL<-c(1)

    colnames(CountedIdentifiedGenes)<-c("Chromosome", "Leftmost_Marker",
                                        "Rightmost_Marker", "Trait",'Treatment',
                                        "Length","QTL_Type","N_Genes", "Number_Trait_QTL")
    return(CountedIdentifiedGenes)
 } }



#' Validate QTL using a list of previously known genes
#'
#' This function determines the number of genes expected to be identified by a QTL
#' of a particular length, given the distribution of genes and mapping markers in
#' the genome. Relies on the Gene_List_Groomer, GeneCounter, MarkerSectioner,
#' and QTL_Placement_Probabilities functions.
#'
#' @param QTL_with_Metadata A 7 column matrix containing the chromosome
#' on which the QTL was originally found ('Chromosome'), the Leftmost
#'  Confidence interval's location,the Rightmost Conficence interval's location,
#'  the trait for which the QTL were discovered, the treatment that was applied,
#'  the mapping method, the overall type of experiment, and the lengths of the QTL.
#' @param number_repetitions An integer defining the number of simulations to run.
#' @param gene_list A 4 column matrix containing a list of known genes.
#' Columns should include "GeneID', "Trait", "Chromosome", and "Start_Site".
#' @param Placement_Type Either 'extension' or 'centered'. This depends on the
#' mapping style of the original QTL experiment.
#' @param Trait The trait for which you would like to analyze your QTL.
#' @param MarkerList A 3 column matrix containing the marker loci
#' (ID, Chromosome, and Base)).
#' @param chromosome_size A matrix containing the number of each chromosome,
#'  their total lengths, and the first and last markers used for the original
#'  mapping experiment
#' @param WholeGenomeGeneDistribution A 4 column matrix of all genes in the
#'  genome. Include Chromosome, GeneStart, GeneEnd, and GeneMiddle. The SPQValidate
#'  function uses GeneMiddle to simulate gene placement.
#' @return A matrix containing the lengths of QTL (as an identifier), the true observed
#' value for identified gene number, the mean value resulting from the simulations,
#' and the upper and lower 95% confidence intervals for the simulation distributions.
#' @export
SPQValidate<-function(QTL_with_Metadata,
                     Trait,
                     number_repetitions,
                     gene_list,
                     chromosome_size,
                     MarkerList,
                     WholeGenomeGeneDistribution,
                     Placement_Type){
  TrueGeneList<-gene_list
  QTLofInterest<-GeneCounter(QTL_with_Metadata,gene_list,Trait,Placement_Type,MarkerList)

  if(length(QTLofInterest)<=1){
    return(QTLofInterest)
  }
  gene_list<-Gene_List_Groomer(MarkerList,gene_list,Placement_Type,Trait)
  gene_list<-subset(gene_list,select=c(GeneID,Chromosome,Start_Site))
  gene_list$EGN<-c(0)
  QTLofInterest<-QTLofInterest[grep(Trait,QTLofInterest$Trait),]
  num_QTL<-length(QTLofInterest$Chromosome)
  if(length(QTLofInterest$Length)==0){
    QTLofInterest$Length<-as.numeric(QTLofInterest$Rightmost_Marker)-as.numeric(QTLofInterest$Leftmost_Marker)

  }
  num_Genes<-length(gene_list$GeneID)
  number_chromosomes<-length(chromosome_size[,1])
  toOutput<-as.data.frame(matrix(nrow= num_QTL,ncol=number_repetitions,data=0))
  pb   <- txtProgressBar(0, 1, style=3)


  Sectioned_List<-Marker_Sectioner(number_chromosomes,MarkerList)
  Probabilities<-QTL_Placement_Probabilities(QTLofInterest,chromosome_size,Sectioned_List,Placement_Type)

  if(Placement_Type=='centered'){
    for(QTL in 1:num_QTL){
      half_QTL_Length<-as.numeric(QTLofInterest$Length[QTL])/2
      for(SimRun in 1:number_repetitions){
        for(Gene in 1:num_Genes){
          SelectedLocus<-as.numeric(as.character(sample(1:length(WholeGenomeGeneDistribution$GeneMiddle),1)))
          gene_list$Start_Site[Gene]<-as.numeric(as.character(WholeGenomeGeneDistribution$GeneMiddle[SelectedLocus]))
          gene_list$Chromosome[Gene]<-as.numeric(as.character(WholeGenomeGeneDistribution$Chromosome[SelectedLocus]))
          BP_locus<-gene_list$Start_Site[Gene]
          Chr_number<-gene_list$Chromosome[Gene]

          MarkersOnChromosome<-as.data.frame(as.matrix(eval(as.name(Sectioned_List[Chr_number]))))

          chromosomeMarkerPositions<-as.integer(as.character(MarkersOnChromosome$Base))

          knownLikelihood<-as.numeric(Probabilities[QTL])


          RangeL<- BP_locus-half_QTL_Length
          RangeR<- BP_locus+half_QTL_Length
          MarkersL<-MarkersOnChromosome[which(chromosomeMarkerPositions + half_QTL_Length <= as.numeric(chromosome_size$Last.Markers[Chr_number]) &
                                                chromosomeMarkerPositions <= BP_locus &
                                                chromosomeMarkerPositions >= RangeL),3]
          LengthMarkersL<-as.numeric(length(MarkersL))

          MarkersR<-MarkersOnChromosome[which(chromosomeMarkerPositions - half_QTL_Length >= as.numeric(chromosome_size$First.Markers[Chr_number]) &
                                                chromosomeMarkerPositions > as.numeric(BP_locus) &  # it's not >= so we don't double-count locus
                                                chromosomeMarkerPositions <= RangeR ),3]
          LengthMarkersR<-as.numeric(length(MarkersR))

          AvailableMarkers <- LengthMarkersL + LengthMarkersR
          Number_Expected<- AvailableMarkers*knownLikelihood
          if(length(Number_Expected)==0){
            Number_Expected<-0
          }
          gene_list$EGN[Gene]<-gene_list$EGN[Gene]+ Number_Expected
        }
        SimEGN<-sum(gene_list$EGN)
        toOutput[QTL,SimRun]<-SimEGN
        gene_list$EGN<-0
      }
      setTxtProgressBar(pb, QTL/length(QTLofInterest$Chromosome))
    }
  }
  if(Placement_Type=='extension'){
    for(QTL in 1:num_QTL){
      QTL_Length<-as.numeric(QTLofInterest$Length[QTL])
      for(SimRun in 1:number_repetitions){
        for(Gene in 1:num_Genes){
          SelectedLocus<-as.numeric(as.character(sample(1:length(WholeGenomeGeneDistribution$GeneMiddle),1)))
          gene_list$Start_Site[Gene]<-as.numeric(as.character(WholeGenomeGeneDistribution$GeneMiddle[SelectedLocus]))
          gene_list$Chromosome[Gene]<-as.numeric(as.character(WholeGenomeGeneDistribution$Chromosome[SelectedLocus]))
          BP_locus<-gene_list$Start_Site[Gene]
          Chr_number<-gene_list$Chromosome[Gene]

          MarkersOnChromosome<-as.data.frame(as.matrix(eval(as.name(Sectioned_List[Chr_number]))))

          chromosomeMarkerPositions<-as.integer(as.character(MarkersOnChromosome$Base))

          knownLikelihood<-as.numeric(Probabilities[QTL])


          RangeL<- BP_locus-QTL_Length
          RangeR<- BP_locus+QTL_Length
          MarkersL<-MarkersOnChromosome[which(chromosomeMarkerPositions + QTL_Length <= as.numeric(chromosome_size$Last.Markers[Chr_number]) &
                                                chromosomeMarkerPositions <= BP_locus &
                                                chromosomeMarkerPositions >= RangeL),3]
          LengthMarkersL<-as.numeric(length(MarkersL))

          MarkersR<-MarkersOnChromosome[which(chromosomeMarkerPositions - QTL_Length >= as.numeric(chromosome_size$First.Markers[Chr_number]) &
                                                chromosomeMarkerPositions > as.numeric(BP_locus) &  # it's not >= so we don't double-count locus
                                                chromosomeMarkerPositions <= RangeR ),3]
          LengthMarkersR<-as.numeric(length(MarkersR))

          AvailableMarkers <- LengthMarkersL + LengthMarkersR
          Number_Expected<- AvailableMarkers*knownLikelihood
          if(length(Number_Expected)==0){
            Number_Expected<-0
          }
          gene_list$EGN[Gene]<-gene_list$EGN[Gene]+ Number_Expected
        }
        SimEGN<-sum(gene_list$EGN)
        toOutput[QTL,SimRun]<-SimEGN
        gene_list$EGN<-0
      }
      setTxtProgressBar(pb, QTL/length(QTLofInterest$Chromosome))
    }
  }



  SDs<-apply(toOutput,1,sd)
  upper<-rowMeans(toOutput)+1.96*SDs
  lower<-rowMeans(toOutput)-1.96*SDs
  CIs<-as.data.frame(as.matrix(cbind(QTLofInterest$Length,QTLofInterest$N_Genes,rowMeans(toOutput),SDs,lower,upper)))
  colnames(CIs)<-c("QTL","Observed Value","Mean","SEM","Lower 95% CI","Upper 95% CI")
  return(CIs)


}



