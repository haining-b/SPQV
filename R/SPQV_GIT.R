
#' Create chromosome specific matrices of QTL mapping marker locations
#'
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

  ############################
  ######Checking inputs#######
  ############################

  colnames(MarkerList)<-c("ID","Chromosome","Base")

  #tests Chromosome column before allowing a full run
  if(typeof(MarkerList$Chromosome[1])!='integer' &
     typeof(MarkerList$Chromosome[1])!='numeric'){
    warning("Chromosome column contains neither integer nor numeric types.")

  }
  if(mean(nchar(as.character(MarkerList$Chromosome)))>4){
    warning("Chromosome number is unusually large. Column may represent base instead.")
  }

  #tests base column before allowing a full run
  if(typeof(MarkerList$Base[1])!='integer' &
     typeof(MarkerList$Base[1])!='numeric'){
    warning("Base column contains neither integer nor numeric types.")

  }
  if(mean(nchar(as.character(MarkerList$Base)))<3){
    warning("Base number is unusually small. Column may represent chromosome instead.")
  }

  #Testing number_chromosomes
  Apparent_N_Chr<-as.numeric(MarkerList$Chromosome[order(MarkerList$Chromosome)[length(MarkerList$Chromosome)]])
  number_chromosomes<-as.numeric(number_chromosomes)
  if(Apparent_N_Chr != number_chromosomes){
    warning('Input value for number_chromosomes not congruent with MarkerList.')
  }
  ##############
  ##############
  ##############


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
#' Columns should include "GeneID', "Trait", "Chromosome", and "Locus".
#' @param Placement_Type Either 'extension' or 'centered'. This depends on the
#' mapping style of the original QTL experiment.
#' @param Trait The trait for which you would like to produce a gene list.
#'
#' @return A 4 column matrix containing "GeneID', "Trait", "Chromosome",
#' and "Locus" for each remaining gene in the trait specified.
#' @export

Gene_List_Groomer<-function(MarkerList,gene_list,Placement_Type,Trait){

  ############################
  ######Checking inputs#######
  ############################

  #Check Placement_Type
  #####
  if(typeof(Placement_Type) != 'character'){
    stop("Placement_Type must be either 'extension' or 'centered'.")
  }
  Placement_Type<-tolower(Placement_Type)
  if(Placement_Type!='extension' & Placement_Type!='centered'){
    stop("Placement_Type must be either 'extension' or 'centered'.")
  }
  temp_gene_list<-unique(gene_list)
  if(length(temp_gene_list)<4){
    stop("Gene list does not contain all required columns.")
  }

  colnames(temp_gene_list)<-c('GeneID',"Trait","Chromosome","Locus")
  #####

  #gene_list testing
  #####
  #tests Chromosome column before allowing a full run
  if(typeof(temp_gene_list$Chromosome[1])=='double'){
    temp_gene_list$Chromosome<-as.integer(temp_gene_list$Chromosome)
  }
  if(typeof(temp_gene_list$Chromosome[1])!='integer' &
     typeof(temp_gene_list$Chromosome[1])!='numeric'){
    warning("gene_list Chromosome column contains neither integer nor numeric types.")

  }
  if(mean(nchar(temp_gene_list$Chromosome))>4){
    warning("Chromosome number is unusually large. Column may represent base instead.")
  }

  #tests Locus column before allowing a full run
  if(typeof(temp_gene_list$Locus[1])=='double'){
    temp_gene_list$Locus<-as.integer(temp_gene_list$Locus)
  }
  if(typeof(temp_gene_list$Locus[1])!='integer' &
     typeof(temp_gene_list$Locus[1])!='numeric'){
    warning("gene_list Locus column contains neither integer nor numeric types.")

  }
  if(mean(nchar(temp_gene_list$Locus))<3){
    warning("gene_list Locus number is unusually small. Column may represent chromosome instead.")
  }
  #####

   colnames(MarkerList)<-c("ID","Chromosome","Base")
  #MarkerList testing
  #####
  #tests Chromosome column before allowing a full run
  if(typeof(MarkerList$Chromosome[1])=='double'){
    MarkerList$Chromosome<-as.integer(MarkerList$Chromosome)
  }
  if(typeof(MarkerList$Chromosome[1])!='integer' &
     typeof(MarkerList$Chromosome[1])!='numeric'){
    warning("MarkerList Chromosome column contains neither integer nor numeric types.")

  }
  if(mean(nchar(as.character(MarkerList$Chromosome)))>4){
    warning("MarkerList Chromosome number is unusually large. Column may represent base instead.")
  }

  #tests base column before allowing a full run
  if(typeof(MarkerList$Base[1])=='double'){
    MarkerList$Base<-as.integer(MarkerList$Base)
  }
  if(typeof(MarkerList$Base[1])!='integer' &
     typeof(MarkerList$Base[1])!='numeric'){
    warning("MarkerList Base column contains neither integer nor numeric types.")

  }
  if(mean(nchar(as.character(MarkerList$Base)))<3){
    warning("MarkerList Base number is unusually small. Column may represent chromosome instead.")
  }

  #####

  #Trait testing
  #####
  if(typeof(Trait)!=typeof(gene_list$Trait[1])){
    stop('Trait input is not the same type as gene_list Trait. Check both type and
         gene_list column order.')
  }
  if(typeof(Trait)=='character'&typeof(temp_gene_list$Trait[1])=='character'){
    Trait<-tolower(Trait)
    temp_gene_list$Trait<-tolower(temp_gene_list$Trait)
    gene_list$Trait<-tolower(gene_list$Trait)
  }
  #####
   ###########
   ###########
   ###########

  temp_gene_list<-temp_gene_list[grep(Trait,temp_gene_list$Trait),]


  if(length(temp_gene_list$Trait)==0){

    return('No genes for this trait.')
  }
  if(length(temp_gene_list$Trait)==1){
    return(temp_gene_list)
  }
  if(Placement_Type=='centered'){
    Gene_List<-temp_gene_list
  }
  if(Placement_Type=='extension'){

    colnames(MarkerList)<-c('GeneID',"Chromosome","Locus")
    MarkerList$Trait<-"Marker"
    MarkerList<-MarkerList[,c(colnames(temp_gene_list))]
    LociList<-rbind(MarkerList,temp_gene_list)
    LociList<-LociList[order(LociList$Chromosome,LociList$Locus),]
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
#' @param chromosome_size A matrix containing the number of each chromosome
#'  and the first and last markers used for the original mapping experiment.
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

  ############################
  ######Checking inputs#######
  ############################

  #Placement_Type Checks
  #####
  if(typeof(Placement_Type) != 'character'){
    stop("Placement_Type must be either 'extension' or 'centered'.")
  }
  Placement_Type<-tolower(Placement_Type)
  if(Placement_Type!='extension' & Placement_Type!='centered'){
    stop("Placement_Type must be either 'extension' or 'centered'.")
  }
  #####

  #trait type checks
  #######
  if(typeof(Trait)=='character'&typeof(QTLofInterest$Trait[1])=='character'){
    Trait<-tolower(Trait)
    QTLofInterest$Trait<-tolower(QTLofInterest$Trait)
  }


  #######
  colnames(chromosome_size)<-c("Number","First.Marker","Last.Marker")

  #Checking chromosome_size
  #####
  if(typeof(chromosome_size$First.Marker[1])=='double'){
    chromosome_size$First.Marker<-as.integer(chromosome_size$First.Marker)
  }
  if(typeof(chromosome_size$First.Marker[1])!='integer' &
     typeof(chromosome_size$First.Marker[1])!='numeric'){
    warning("chromosome_size First.Marker column contains neither integer nor numeric types.")
  }
  if(typeof(chromosome_size$Last.Marker[1])=='double'){
    chromosome_size$Last.Marker<-as.integer(chromosome_size$Last.Marker)
  }
  if(typeof(chromosome_size$Last.Marker[1])!='integer' &
     typeof(chromosome_size$Last.Marker[1])!='numeric'){
    warning("chromosome_size Last.Marker column contains neither integer nor numeric types.")
  }
  if(typeof(chromosome_size$Number[1])=='double'){
    chromosome_size$Number<-as.integer(chromosome_size$Number)
  }
  if(typeof(chromosome_size$Number[1])!='integer' &
     typeof(chromosome_size$Number[1])!='numeric'){
    stop("Please make chromosome_size Number column numeric.")
  }
  if(mean(chromosome_size$First.Marker)>mean(chromosome_size$Last.Marker)){
    stop("chromosome_size Last.Marker column contains lower values than
            First.Marker column.")
  }
  #####


  chromosome_size<-chromosome_size[order(as.numeric(chromosome_size$Number)),]
  number_chromosomes<-as.numeric(length(chromosome_size[,1]))


  colnames(QTLofInterest)[1:3]<-c('Chromosome',"Leftmost_Marker",
                                  "Rightmost_Marker")
  #Checking QTLofInterest
  #####
  if(typeof(QTLofInterest$Leftmost_Marker[1])=='double'){
    QTLofInterest$Leftmost_Marker<-as.integer(QTLofInterest$Leftmost_Marker)
  }
  if(typeof(QTLofInterest$Leftmost_Marker[1])!='integer' &
     typeof(QTLofInterest$Leftmost_Marker[1])!='numeric'){
    warning("QTLofInterest Leftmost_Marker column contains neither integer nor numeric types.")
  }
  if(typeof(QTLofInterest$Rightmost_Marker[1])=='double'){
    QTLofInterest$Rightmost_Marker<-as.integer(QTLofInterest$Rightmost_Marker)
  }
  if(typeof(QTLofInterest$Rightmost_Marker[1])!='integer' &
     typeof(QTLofInterest$Rightmost_Marker[1])!='numeric'){
    warning("QTLofInterest Rightmost_Marker column contains neither integer nor numeric types.")
  }
  if(typeof(QTLofInterest$Chromosome[1])=='double'){
    QTLofInterest$Chromosome<-as.integer(QTLofInterest$Chromosome)
  }
  if(typeof(QTLofInterest$Chromosome[1])!='integer' &
     typeof(QTLofInterest$Chromosome[1])!='numeric'){
    warning("QTLofInterest Chromosome column contains neither integer nor numeric types.")
  }

  #####
  QTLofInterest$Chromosome<-as.numeric(QTLofInterest$Chromosome)
  QTLofInterest$Leftmost_Marker<-as.numeric(QTLofInterest$Leftmost_Marker)
  QTLofInterest$Rightmost_Marker<-as.numeric(QTLofInterest$Rightmost_Marker)

  if(length(QTLofInterest$Length)==0){
    QTLofInterest$Length<-as.numeric(QTLofInterest$Rightmost_Marker)-as.numeric(QTLofInterest$Leftmost_Marker)

  }
  if(mean(QTLofInterest$Length)==0){
    QTLofInterest$Length<-as.numeric(QTLofInterest$Rightmost_Marker)-as.numeric(QTLofInterest$Leftmost_Marker)

  }
  PossiblePositions<-as.data.frame(matrix(nrow=length(QTLofInterest$Chromosome),ncol=1+number_chromosomes))
  colnames(PossiblePositions)<-c("Length_QTL",1:number_chromosomes)
  PossiblePositions$Length_QTL<-QTLofInterest$Length

  if(Placement_Type=='centered'){
    for (QTL in 1:length(QTLofInterest$Length)){
      half_QTL_Length<-as.numeric(as.character(QTLofInterest$Length[QTL]))/2

      for (Chromosome in  1:number_chromosomes){
        FirstMarker<-as.numeric(as.character(chromosome_size$First.Marker[Chromosome]))
        LastMarker<-as.numeric(as.character(chromosome_size$Last.Marker[Chromosome]))

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
        FirstMarker<-as.numeric(as.character(chromosome_size$First.Marker[Chromosome]))
        LastMarker<-as.numeric(as.character(chromosome_size$Last.Marker[Chromosome]))

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
#' Columns should include "GeneID', "Trait", "Chromosome", and "Locus".
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
  #Trait testing
  #####
  if(typeof(Trait)!=typeof(gene_list$Trait[1])){
    stop('Trait input is not the same type as gene_list Trait. Check both type and
         gene_list column order.')
  }
  if(typeof(Trait)=='character'&typeof(gene_list$Trait[1])=='character'){
    Trait<-tolower(Trait)
    gene_list$Trait<-tolower(gene_list$Trait)
    QTL_with_Metadata$Trait<-tolower(QTL_with_Metadata$Trait)
  }
  ######
  colnames(QTL_with_Metadata)<-c('Chromosome',"Leftmost_Marker",
                                 "Rightmost_Marker","Trait","Treatment",
                                 "QTL_Type")

  QTL_with_Metadata<-QTL_with_Metadata[,c(1:6)]
  #Everything except QTL_with_Metadata gets checked here:
  temp_gene_list<-Gene_List_Groomer(MarkerList,gene_list,Placement_Type,Trait)
  if(length(temp_gene_list)< 2 ){

    return(temp_gene_list)

  }

  # Checking QTL_with_Metadata
  #####
  if(typeof(QTL_with_Metadata$Leftmost_Marker[1])=='double'){
    QTL_with_Metadata$Leftmost_Marker<-as.integer(QTL_with_Metadata$Leftmost_Marker)
  }
  if(typeof(QTL_with_Metadata$Leftmost_Marker[1])!='integer' &
     typeof(QTL_with_Metadata$Leftmost_Marker[1])!='numeric'){
    warning("QTL_with_Metadata Leftmost_Marker column contains neither integer nor numeric types.")
  }
  if(typeof(QTL_with_Metadata$Rightmost_Marker[1])=='double'){
    QTL_with_Metadata$Rightmost_Marker<-as.integer(QTL_with_Metadata$Rightmost_Marker)
  }
  if(typeof(QTL_with_Metadata$Rightmost_Marker[1])!='integer' &
     typeof(QTL_with_Metadata$Rightmost_Marker[1])!='numeric'){
    warning("QTL_with_Metadata Rightmost_Marker column contains neither integer nor numeric types.")
  }
  if(typeof(QTL_with_Metadata$Chromosome[1])=='double'){
    QTL_with_Metadata$Chromosome<-as.integer(QTL_with_Metadata$Chromosome)
  }
  if(typeof(QTL_with_Metadata$Chromosome[1])!='integer' &
     typeof(QTL_with_Metadata$Chromosome[1])!='numeric'){
    warning("QTL_with_Metadata Chromosome column contains neither integer nor numeric types.")
  }

  if(mean(nchar(QTL_with_Metadata$Chromosome))>4){
    warning("QTL_with_Metadata Chromosome number is unusually large. Column may represent confidence interval instead.")
  }

  if(mean(nchar(QTL_with_Metadata$Leftmost_Marker))<4){
    warning("QTL_with_Metadata Leftmost_Marker number is unusually small. Column may represent Chromosome instead.")
  }

  if(mean(nchar(QTL_with_Metadata$Rightmost_Marker))<4){
    warning("QTL_with_Metadata Rightmost_Marker number is unusually small. Column may represent Chromosome instead.")
  }
  if(mean(QTL_with_Metadata$Leftmost_Marker)>mean(QTL_with_Metadata$Rightmost_Marker)){
    warning("QTL_with_Metadata Leftmost_Marker is larger than Rightmost_Marker. Check column order.")
  }
  if(typeof(QTL_with_Metadata$QTL_Type)!='character'){
    warning("QTL_with_Metadata$QTL_Type should contain character values to ensure appropriate column use.")
  }
  if(typeof(QTL_with_Metadata$Trait)!='character'){
    warning("QTL_with_Metadata$Trait should contain character values to ensure appropriate column use.")
  }
  if(typeof(QTL_with_Metadata$Treatment)!='character'){
    warning("QTL_with_Metadata$Treatment should contain character values to ensure appropriate column use.")
  }


  #####


  if(as.numeric(as.character(length(QTL_with_Metadata[,1])))>1){
    TraitsandTreatments<-subset(QTL_with_Metadata, select= c(Trait, Treatment))
    sepTreatments_QTL_Count<-plyr::ddply(TraitsandTreatments,.(TraitsandTreatments$Trait,TraitsandTreatments$Treatment),nrow)
    colnames(sepTreatments_QTL_Count)<-c("Trait",'Treatment','Frequency')



    QTLData<-QTL_with_Metadata[grep(Trait, TraitsandTreatments$Trait),]
    if(length(QTLData$Chromosome)==0){
      return("No QTL found for this trait.")
    }

    QTLData$N_QTL<-c(0)
    count<-c()

    for(identified in 1:length(QTLData$Chromosome)){
      sep_Treat_QTL_4_u<-sepTreatments_QTL_Count$Frequency[which(sepTreatments_QTL_Count$Trait==QTLData$Trait[identified] &
                                                                   sepTreatments_QTL_Count$Treatment==QTLData$Treatment[identified])]
      count<-c(count,sep_Treat_QTL_4_u)
    }

    QTLData$N_QTL<-count
    QTLData$N_Genes<-c(0)

    identified_genes<-data.frame(matrix(ncol=length(QTLData)))
    colnames(identified_genes)<-colnames(QTLData)


    for (i in 1:length(QTLData$Chromosome)){
      QTLchromosome<-as.numeric(as.character(QTLData$Chromosome[i]))
      QTLLCI<-as.numeric(as.character(QTLData$Leftmost_Marker[i]))
      QTLRCI<-as.numeric(as.character(QTLData$Rightmost_Marker[i]))

      for (Gene in 1: length(temp_gene_list$GeneID)){
        if(as.numeric(as.character(temp_gene_list$Chromosome[Gene]))==as.numeric(as.character(QTLchromosome))){
          GeneStartSite<-as.numeric(temp_gene_list$Locus[Gene])
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
    identified_genes<-identified_genes[,c("Chromosome", "Leftmost_Marker",
                         "Rightmost_Marker", "Trait",'Treatment',
                         "Length","QTL_Type","N_Genes",'N_QTL')]
    colnames(identified_genes)[9]<-'Number_Trait_QTL'
    colnames(QTLData)[7]<-'Number_Trait_QTL'
    identified_genes2<-identified_genes

    if(length(identified_genes[,1])>1){
      for(possible_Duplicates in 1:(length(identified_genes[,1])-1)){
        if(all(identified_genes[possible_Duplicates,c(1:6)]==identified_genes[possible_Duplicates+1,c(1:6)])){
          identified_genes2[possible_Duplicates,]<-0
        }
      }
      if(length(which(identified_genes2$N_Genes==0))>0){
        CountedIdentifiedGenes<-identified_genes2[-which(identified_genes2$N_Genes==0),]
      } else{
        CountedIdentifiedGenes<-identified_genes2
      }

      CountedIdentifiedGenesOutput<-CountedIdentifiedGenes
      NoGenes<-dplyr::setdiff(QTLData[1:7],CountedIdentifiedGenes[1:7])

      if(length(NoGenes$Chromosome)>0){
        NoGenes$Length<-as.numeric(as.character(NoGenes$Rightmost_Marker))-as.numeric(as.character(NoGenes$Leftmost_Marker))
        NoGenes$N_Genes<-0
        NoGenesOutput<-NoGenes
        colnames(NoGenesOutput)<-c("Chromosome", "Leftmost_Marker",
                                   "Rightmost_Marker", "Trait",
                                   "Treatment","QTL_Type","Number_Trait_QTL",
                                   "Length", "N_Genes")
        colnames(CountedIdentifiedGenesOutput)<-c("Chromosome", "Leftmost_Marker",
                                                  "Rightmost_Marker", "Trait", "Treatment", "QTL_Type",
                                                  "Number_Trait_QTL",
                                                  "N_Genes", "Length")


        CountedIdentifiedGenesOutput2<-rbind(CountedIdentifiedGenesOutput,NoGenesOutput)
      }else{
        CountedIdentifiedGenesOutput2<-CountedIdentifiedGenesOutput
      }

      return(CountedIdentifiedGenesOutput2)
    }


    if(length(identified_genes[,1])==1){


      CountedIdentifiedGenesOutput<-identified_genes[,c(1:5,7,9,8,6)]
      toCompareID<-CountedIdentifiedGenesOutput[,c(1:7)]
      toCompareQTLData<-QTLData[,c(1:7)]
      NoGenes<-dplyr::setdiff(toCompareQTLData,toCompareID)


      if(length(NoGenes$Chromosome)>0){
        NoGenes$N_Genes<-c(0)
        NoGenes$Length<-as.numeric(as.character(NoGenes$Rightmost_Marker))-as.numeric(as.character(NoGenes$Leftmost_Marker))
        NoGenesOutput<-NoGenes

        colnames(NoGenesOutput)<-c("Chromosome", "Leftmost_Marker",
                                   "Rightmost_Marker", "Trait",
                                   "Treatment","QTL_Type","Number_Trait_QTL",
                                    "N_Genes","Length")
        colnames(CountedIdentifiedGenesOutput)<-c("Chromosome", "Leftmost_Marker",
                                                  "Rightmost_Marker", "Trait", "Treatment", "QTL_Type",
                                                  "Number_Trait_QTL",
                                                  "N_Genes", "Length")


        CountedIdentifiedGenesOutput2<-rbind(CountedIdentifiedGenesOutput,NoGenesOutput)
        CountedIdentifiedGenesOutput2<-CountedIdentifiedGenesOutput2[,c(1:6,9,7:8)]
        return(CountedIdentifiedGenesOutput2)
      }else{
        CountedIdentifiedGenesOutput2<-CountedIdentifiedGenesOutput
        return(CountedIdentifiedGenesOutput2)
      }
      return(CountedIdentifiedGenesOutput2)
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
        GeneStartSite<-as.numeric(temp_gene_list$Locus[Gene])
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
      CountedIdentifiedGenes<-identified_genes2[,c("Chromosome", "Leftmost_Marker",
                                        "Rightmost_Marker", "Trait",'Treatment',
                                                   "Length","QTL_Type","N_Genes",'N_QTL')]
      colnames(CountedIdentifiedGenes)[9]<-'Number_Trait_QTL'
    }
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
#'  Confidence interval's location,the Rightmost Confidence interval's location,
#'  the trait for which the QTL were discovered, the treatment that was applied,
#'  the mapping method, the overall type of experiment, and the lengths of the QTL.
#' @param number_repetitions An integer defining the number of simulations to run.
#' @param gene_list A 4 column matrix containing a list of known genes.
#' Columns should include "GeneID', "Trait", "Chromosome", and "Locus".
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
#' @param Simulation_Environment A new environment in which a dataframe containing the
#' results of each simulation is stored. This dataframe can be accessed as
#' Your_Environment$Simulation_Dataframe.
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
                     Placement_Type,
                     Simulation_Environment){

  ############################
  ######Checking inputs#######
  ############################

  #Trait testing
  #####
  if(typeof(Trait)!=typeof(gene_list$Trait[1])){
    stop('Trait input is not the same type as gene_list Trait. Check both type and
         gene_list column order.')
  }
  if(typeof(Trait)=='character'&typeof(gene_list$Trait[1])=='character'){
    Trait<-tolower(Trait)
    gene_list$Trait<-tolower(gene_list$Trait)
    QTL_with_Metadata$Trait<-tolower(QTL_with_Metadata$Trait)
  }
  ######

  #Checking WholeGenomeGeneDistribution
  #####
  colnames(WholeGenomeGeneDistribution)<-c("Chromosome","GeneStart","GeneEnd","GeneMiddle")
  for(i in 1:length(WholeGenomeGeneDistribution)){
    if(typeof(WholeGenomeGeneDistribution[,i])!='integer'&
       typeof(WholeGenomeGeneDistribution[,i])!='numeric'&
       typeof(WholeGenomeGeneDistribution[,i])!='double'){
      to_warn<-paste(colnames(WholeGenomeGeneDistribution[,i])," must contain
            numeric, integer, or double values.")
      stop(to_warn)
    }
  }

  if(mean(nchar(WholeGenomeGeneDistribution$Chromosome))>4){
    warning("WholeGenomeGeneDistribution Chromosome number is unusually large. Column may represent base pair position instead.")
  }

  if(mean(nchar(WholeGenomeGeneDistribution$GeneStart))<4){
    warning("WholeGenomeGeneDistribution GeneStart number is unusually small. Column may represent Chromosome instead.")
  }

  if(mean(nchar(WholeGenomeGeneDistribution$GeneEnd))<4){
    warning("WholeGenomeGeneDistribution GeneEnd number is unusually small. Column may represent Chromosome instead.")
  }
  if(mean(WholeGenomeGeneDistribution$GeneStart)>mean(WholeGenomeGeneDistribution$GeneEnd)){
    warning("WholeGenomeGeneDistribution GeneStart is larger than GeneEnd. Check column order.")
  }
  #####

  #Checking Simulation_Environment
  #####
  if(typeof(Simulation_Environment)!='environment'){
    stop("Simulation_Environment must be type 'environment'. Use the function
         new.env() to produce an environment.")
  }
  #####
  #Checking number_repetitions
  #####
  if(typeof(number_repetitions)=='character'){
    if(typeof(as.numeric(number_repetitions)) =='numeric'|typeof(as.numeric(number_repetitions)) =='double'){
      number_repetitions<-as.numeric(number_repetitions)
    } else{
      stop("number_repetitions must be a numeric or integer value (i.e. 1000).")
    }
  }
  #####
  #####

  # Checking QTL_with_Metadata
  #####
  if(typeof(QTL_with_Metadata$Leftmost_Marker[1])=='double'){
    QTL_with_Metadata$Leftmost_Marker<-as.integer(QTL_with_Metadata$Leftmost_Marker)
  }
  if(typeof(QTL_with_Metadata$Leftmost_Marker[1])!='integer' &
     typeof(QTL_with_Metadata$Leftmost_Marker[1])!='numeric'){
    warning("QTL_with_Metadata Leftmost_Marker column contains neither integer nor numeric types.")
  }
  if(typeof(QTL_with_Metadata$Rightmost_Marker[1])=='double'){
    QTL_with_Metadata$Rightmost_Marker<-as.integer(QTL_with_Metadata$Rightmost_Marker)
  }
  if(typeof(QTL_with_Metadata$Rightmost_Marker[1])!='integer' &
     typeof(QTL_with_Metadata$Rightmost_Marker[1])!='numeric'){
    warning("QTL_with_Metadata Rightmost_Marker column contains neither integer nor numeric types.")
  }
  if(typeof(QTL_with_Metadata$Chromosome[1])=='double'){
    QTL_with_Metadata$Chromosome<-as.integer(QTL_with_Metadata$Chromosome)
  }
  if(typeof(QTL_with_Metadata$Chromosome[1])!='integer' &
     typeof(QTL_with_Metadata$Chromosome[1])!='numeric'){
    warning("QTL_with_Metadata Chromosome column contains neither integer nor numeric types.")
  }

  if(mean(nchar(QTL_with_Metadata$Chromosome))>4){
    warning("QTL_with_Metadata Chromosome number is unusually large. Column may represent confidence interval instead.")
  }

  if(mean(nchar(QTL_with_Metadata$Leftmost_Marker))<4){
    warning("QTL_with_Metadata Leftmost_Marker number is unusually small. Column may represent Chromosome instead.")
  }

  if(mean(nchar(QTL_with_Metadata$Rightmost_Marker))<4){
    warning("QTL_with_Metadata Rightmost_Marker number is unusually small. Column may represent Chromosome instead.")
  }
  if(mean(QTL_with_Metadata$Leftmost_Marker)>mean(QTL_with_Metadata$Rightmost_Marker)){
    warning("QTL_with_Metadata Leftmost_Marker is larger than Rightmost_Marker. Check column order.")
  }

  if(typeof(QTL_with_Metadata$Trait)!='character'){
    warning("QTL_with_Metadata$Trait should contain character values to ensure appropriate column use.")
  }
  if(typeof(QTL_with_Metadata$Treatment)!='character'){
    warning("QTL_with_Metadata$Treatment should contain character values to ensure appropriate column use.")
  }

  #Placement_Type Checks
  #####
  if(typeof(Placement_Type) != 'character'){
    stop("Placement_Type must be either 'extension' or 'centered'.")
  }
  Placement_Type<-tolower(Placement_Type)
  if(Placement_Type!='extension' & Placement_Type!='centered'){
    stop("Placement_Type must be either 'extension' or 'centered'.")
  }
  #####


  TrueGeneList<-gene_list
  QTLofInterest<-GeneCounter(QTL_with_Metadata,gene_list,Trait,Placement_Type,MarkerList)

  if(length(QTLofInterest)<=1){
    return(QTLofInterest)
  }
  gene_list<-Gene_List_Groomer(MarkerList,gene_list,Placement_Type,Trait)
  gene_list<-subset(gene_list,select=c(GeneID,Chromosome,Locus))
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
          gene_list$Locus[Gene]<-as.numeric(as.character(WholeGenomeGeneDistribution$GeneMiddle[SelectedLocus]))
          gene_list$Chromosome[Gene]<-as.numeric(as.character(WholeGenomeGeneDistribution$Chromosome[SelectedLocus]))
          BP_locus<-gene_list$Locus[Gene]
          Chr_number<-gene_list$Chromosome[Gene]

          MarkersOnChromosome<-as.data.frame(as.matrix(eval(as.name(Sectioned_List[Chr_number]))))

          chromosomeMarkerPositions<-as.integer(as.character(MarkersOnChromosome$Base))

          knownLikelihood<-as.numeric(Probabilities[QTL])


          RangeL<- BP_locus-half_QTL_Length
          RangeR<- BP_locus+half_QTL_Length
          MarkersL<-MarkersOnChromosome[which(chromosomeMarkerPositions + half_QTL_Length <= as.numeric(chromosome_size$Last.Marker[Chr_number]) &
                                                chromosomeMarkerPositions <= BP_locus &
                                                chromosomeMarkerPositions >= RangeL),3]
          LengthMarkersL<-as.numeric(length(MarkersL))

          MarkersR<-MarkersOnChromosome[which(chromosomeMarkerPositions - half_QTL_Length >= as.numeric(chromosome_size$First.Marker[Chr_number]) &
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
          gene_list$Locus[Gene]<-as.numeric(as.character(WholeGenomeGeneDistribution$GeneMiddle[SelectedLocus]))
          gene_list$Chromosome[Gene]<-as.numeric(as.character(WholeGenomeGeneDistribution$Chromosome[SelectedLocus]))
          BP_locus<-gene_list$Locus[Gene]
          Chr_number<-gene_list$Chromosome[Gene]

          MarkersOnChromosome<-as.data.frame(as.matrix(eval(as.name(Sectioned_List[Chr_number]))))

          chromosomeMarkerPositions<-as.integer(as.character(MarkersOnChromosome$Base))

          knownLikelihood<-as.numeric(Probabilities[QTL])


          RangeL<- BP_locus-QTL_Length
          RangeR<- BP_locus+QTL_Length
          MarkersL<-MarkersOnChromosome[which(chromosomeMarkerPositions + QTL_Length <= as.numeric(chromosome_size$Last.Marker[Chr_number]) &
                                                chromosomeMarkerPositions <= BP_locus &
                                                chromosomeMarkerPositions >= RangeL),3]
          LengthMarkersL<-as.numeric(length(MarkersL))

          MarkersR<-MarkersOnChromosome[which(chromosomeMarkerPositions - QTL_Length >= as.numeric(chromosome_size$First.Marker[Chr_number]) &
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
  upper<-c()
  lower<-c()

  for(QTL in 1: length(QTLofInterest$Length)){
    Zvalue<-qnorm(.025/as.numeric(QTLofInterest$Number_Trait_QTL[QTL]),lower.tail=FALSE)
    mu<-rowMeans(toOutput)[QTL]
    upper<-c(upper, mu +Zvalue*SDs[QTL])
    lower<-c(lower,mu -Zvalue*SDs[QTL])
  }


  CIs<-as.data.frame(as.matrix(cbind(QTLofInterest$Length,QTLofInterest$N_Genes,rowMeans(toOutput),SDs,lower,upper)))
  toOutput$QTL_Length<-QTLofInterest$Length
  toOutput$Observed_Value<-QTLofInterest$N_Genes
  toOutput<-toOutput[,c((length(toOutput)-1),length(toOutput),1:(length(toOutput)-2))]
  colnames(toOutput)[3:length(toOutput)]<-paste0("Simulation_Round_",1:(length(toOutput)-2))
  Simulation_Environment$Simulation_DataFrame<-toOutput


  colnames(CIs)<-c("QTL","Observed Value","Mean","SEM","Lower 95% CI","Upper 95% CI")
  return(CIs)


}



