
validateDf <- function(df, req_cols) {
  req_names <- c()
  for (col in req_cols) {
    req_name <- col[1]
    req_type <- col[2]
    if (!(req_name %in% colnames(df))) {
      stop(
        sprintf(
          "Dataframe does not contain required column: %s. \nExpected columns: \n%s \nActual columns: \n%s\n",
          req_name,
          paste(req_cols, collapse = ", "),
          paste(colnames(df), collapse = ", ")
        )
      )
    }
    if (typeof(df[, req_name]) != req_type) {
      stop(
        sprintf(
          "Column %s has wrong type. \nExpected: %s \tActual: %s\n",
          req_name,
          req_type,
          typeof(df[, req_name])
        )
      )
    }
    req_names <- c(req_names, req_name)
  }
  return(df[, req_names])
}


#' Create chromosome specific matrices of QTL mapping marker locations
#'
#'
#' This function splits marker_list into a list pf matrices for each chromosome.
#'
#' @param marker_list A 3-column matrix containing the marker loci used in a QTL mapping experiment
#'   (ID, Chromosome, and Base)).
#'
#' @param num_chromosomes Numeric value refering the the number of chromosomes in the organism of
#'   interest.
#'
#' @return a list of N matrices of marker loci
#' @export
SectionMarkers <- function(marker_list, num_chromosomes) {
  marker_list <- validateDf(marker_list, list(
    c("ID", "character"),
    c("Chromosome", "integer"),
    c("Base", "integer")
  ))

  if (!(typeof(num_chromosomes) %in% c('integer','double','numeric'))) {
    stop("Input num_chromosomes must be of type 'integer','double', or 'numeric'.")
  }

  #Testing num_chromosomes
  apparent_num_chr <- max(marker_list$Chromosome)
  if (apparent_num_chr != num_chromosomes) {
    warning(
      sprintf(
        'Input value %d for num_chromosomes not congruent with marker_list (max chr number %d).',
        num_chromosomes,
        apparent_num_chr
      )
    )
  }

  sectioned_list <- vector("list", num_chromosomes)

  for (chromosome in 1:num_chromosomes) {
    marker_df <- marker_list[which(marker_list$Chromosome == chromosome), ]
    sectioned_list[[chromosome]] <- marker_df
  }
  return(sectioned_list)
}

#' Create known gene lists appropriate for specific traits and your style of QTL mapping
#'
#' This function takes a gene list and refines it to only specific trait you would like to examine.
#' Will remove duplicate rows. Additionally, the function will act to reduce the final list if
#' drop_tandem is true, so as to avoid the artificial inflation of identified genes that linkage
#' causes. This reduction is performed such that a single randomly chosen gene from any set of genes
#' that lie between two consecutive markers remains.
#'
#' @param trait The trait for which you would like to produce a gene list.
#' @param gene_list A matrix containing a list of known genes. Columns must include "ID', "Trait",
#'   "Chromosome", and "Base".
#' @param marker_list A matrix containing the marker loci Required columns: ID, Chromosome, and
#'   Base.
#' @param drop_tandem If true, remove "tandem repeats" - multiple genes with no markers between
#'   them.
#'
#' @return A 4 column matrix containing "ID', "Trait", "Chromosome", and "Base" for each remaining
#'   gene in the trait specified.
#' @export
FilterGeneList <-
  function(trait,
           gene_list,
           marker_list,
           drop_tandem = TRUE) {
    if (typeof(trait) != 'character') {
      stop("Trait input must be of type 'character'.")
    }
    trait <- tolower(trait)
    gene_list$Trait <- tolower(gene_list$Trait)

    gene_list <- validateDf(gene_list, list(
      c("ID", "character"),
      c("Trait", "character"),
      c("Chromosome", "integer"),
      c("Base", "integer")
    ))

    marker_list <- validateDf(marker_list, list(
      c("ID", "character"),
      c("Chromosome", "integer"),
      c("Base", "integer")
    ))

    trait_gene_list <- gene_list[which(gene_list$Trait==trait), ]
    trait_gene_list <- unique(trait_gene_list)

    if (nrow(trait_gene_list) == 0) {
      stop(sprintf('No genes for trait: %s.', trait))
    }
    if ( !drop_tandem | nrow(trait_gene_list) == 1) {
      return(trait_gene_list)
    }

    marker_list$Trait <- "Marker"
    loci_list <- rbind(marker_list, trait_gene_list)
    loci_list <-
      loci_list[order(loci_list$Chromosome, loci_list$Base), ]
    gene_spots<-which(loci_list$Trait==trait)
    tandem_arrays <- c()
    # TODO is equiv to saving prev marker, shuffling, and then dropping dups?
    for (i in 2:length(gene_spots)) {
      if (gene_spots[i] == gene_spots[i - 1] + 1) {
        tandem_arrays[i] <- 1
        tandem_arrays[i - 1] <- 1
      } else{
        tandem_arrays[i] <- 0
      }
    }
    if (is.na(tandem_arrays[1])) {
      tandem_arrays[1] <- 0
    }
    sequential_increases <- rle(tandem_arrays)
    vals <- sequential_increases$values
    longs <- sequential_increases$lengths
    presence_absence <- c()
    for (seq_i in 1:length(vals)) {
      if (vals[seq_i] == 0) {
        presence_absence <- c(presence_absence, rep(1, longs[seq_i]))
      }
      else{
        pos_selected <- sample(1:longs[seq_i], 1)
        low_substitute_zeroes <- rep(0, (pos_selected - 1))
        high_substitute_zeroes <- rep(0, longs[seq_i] - pos_selected)
        all_subs <- c(low_substitute_zeroes, 1, high_substitute_zeroes)
        presence_absence <- c(presence_absence, all_subs)
      }
    }
    return(unique(trait_gene_list[grep(1, presence_absence), ]))
  }


#' Determine the probability of placing a particular QTL on any marker in the genome
#'
#' This function determines the likelihood of placing a QTL on any individual marker in the genome
#' that is reachable by that QTL. If placement_type is "extension", the likelihood reflects the
#' probability of choosing that marker and either L->R or R->L mapping. It takes into account the
#' sizes of your identified QTL, the lengths of each chromosome, the locations of markers, and the
#' type of mapping that originally produced the QTL.
#'
#' @param qtl_list A matrix containing (at minimum) the chromosome on which the QTL was originally
#'   found ('Chromosome'), the Leftmost Confidence interval's location, and the Rightmost Conficence
#'   interval's location.
#' @param placement_type Either 'extension' or 'centered'. This depends on the mapping style of the
#'   original QTL experiment.
#' @param sectioned_marker_list The list of the names of the N matrices of marker loci resulting
#'   from use of the Marker_Sectioner function.
#' @param chromosome_size A matrix containing the number of each chromosome and the first and last
#'   markers used for the original mapping experiment.
#'
#' @return A list with length nrow(qtl_list), representing the probability of choosing any
#'   particular marker reachable with that QTL
#' @export
QTLPlacementProbabilities <-
  function(qtl_list,
           placement_type,
           sectioned_marker_list,
           chromosome_size) {

    qtl_list <- validateDf(qtl_list, list(
      c("Chromosome", "integer"),
      c("LeftmostMarker", "integer"),
      c("RightmostMarker", "integer"),
      c("Trait", "character"),
      c("Length", "integer")
    ))

    chromosome_size <- validateDf(chromosome_size, list(
      c("Chromosome", "integer"),
      c("LeftmostMarker", "integer"),
      c("RightmostMarker", "integer"),
      c("Length", "integer")
    ))

    placement_types <- c("extension", "centered")
    if (!(placement_type %in% placement_types)) {
      stop(
        sprintf(
          "Unrecognized placement_type: %s. \nMust be one of: \n%s \n",
          placement_type,
          paste(placement_types, collapse = ", ")
        )
      )
    }

    # Sort chromosomes
    chromosome_size <- chromosome_size[order(as.numeric(chromosome_size$Chromosome)), ]
    num_chromosomes <- as.numeric(nrow(chromosome_size))

    if (nrow(qtl_list) == 0 | sum(qtl_list$Length) == 0) {
      qtl_list$Length <-
        qtl_list$RightmostMarker - qtl_list$LeftmostMarker

    }

    poss_positions <-
      as.data.frame(matrix(
        nrow = nrow(qtl_list),
        ncol = 1 + num_chromosomes
      ))
    colnames(poss_positions) <- c("QTLLength", 1:num_chromosomes)
    poss_positions$QTLLength <- qtl_list$Length


    for (qtl_i in 1:nrow(qtl_list)) {

      if (placement_type == 'centered') {
        qtl_ext_length = qtl_list$Length[qtl_i] / 2
      } else if (placement_type == 'extension') {
        qtl_ext_length = qtl_list$Length[qtl_i]
      }

      for (chr_i in  1:num_chromosomes) {
        first_marker <-
          chromosome_size$LeftmostMarker[chr_i]
        last_marker <-
          chromosome_size$RightmostMarker[chr_i]

        first_avail_base <- first_marker + qtl_ext_length
        last_avail_base <- last_marker - qtl_ext_length

        chr_markers <-sectioned_marker_list[chr_i]
        chr_markers<-chr_markers[[1]]

        remaining_markers_lr <-
          chr_markers[which(chr_markers$Base <= last_avail_base), ]
        remaining_markers_rl <-
          chr_markers[which(chr_markers$Base >= first_avail_base), ]

        if (placement_type == "extension") {
          remaining_markers <-
            rbind(remaining_markers_lr, remaining_markers_rl)
        } else if (placement_type == "centered") {
          # only use markers that are OK in both directions, and only count each once
          remaining_markers <-
            remaining_markers_lr[which(remaining_markers_lr$Base >= first_avail_base), ]
        }

        length_remaining_markers <- nrow(remaining_markers)

        poss_positions[qtl_i, chr_i + 1] <-
          length_remaining_markers
      }
    }
    marker_probs <- 1 / rowSums(poss_positions[, 2:ncol(poss_positions)])

    return(marker_probs)
  }


#' Count the number of known genes under a particular QTL
#'
#' This function determines the number of genes found by a QTL in the original
#' mapping experiment. Relies on the FilterGeneList function to ensure
#' the appropriate genes are being used.
#'
#' @param qtl_list A matrix containing the chromosome
#' on which the QTL was originally found ('Chromosome'), the Leftmost
#'  Confidence interval's location,the Rightmost Conficence interval's location,
#'  the trait for which the QTL were discovered, and the lengths of the QTL.
#'  Required columns: (Chromosome, LeftmostMarker, RightmostMarker, Trait, Length)
#' @param trait The trait for which you would like to produce a gene list.
#' @param trait_gene_list A 4 column matrix containing a list of known genes.
#' Columns should include "GeneID', "Trait", "Chromosome", and "Locus".
#' @param marker_list A 3 column matrix containing the marker loci
#' (ID, Chromosome, and Base)).
#'
#' TODO fix doc
#'
#' @return A matrix containing the number of genes identified within a QTL
#' for a specific trait, as well as the number of QTL that were identified
#' for that trait
#' @export
CountGenesFound <-
  function(qtl_list,
           trait,
           trait_gene_list,
           marker_list) {

    qtl_list  <- validateDf(qtl_list, list(
      c("Chromosome", "integer"),
      c("LeftmostMarker", "integer"),
      c("RightmostMarker", "integer"),
      c("Trait", "character"),
      c("Length", "integer")
    ))

    if  (typeof(trait) != 'character') {
      stop("Trait input must be of type 'character'.")
    }
    trait <- tolower(trait)
    trait_gene_list$Trait <- tolower(trait_gene_list$Trait)
    qtl_list$Trait <- tolower(qtl_list$Trait)
    if (! all(trait_gene_list$Trait == trait)) {
      stop(sprintf("Gene list should already be filtered to contain only trait %s; instead found %s",
                   trait,
                   unique(trait_gene_list$Trait)
          ))
    }

    if (nrow(trait_gene_list) <= 1) {
      return(trait_gene_list)
    }

    trait_qtl_list <-
      qtl_list[which(qtl_list$Trait==trait), ]
    if (nrow(trait_qtl_list) == 0) {
      stop("No QTL found for this trait.")
    }


    trait_qtl_list$NumGenes <- 0
    trait_qtl_list$FoundGeneIDs <- 0

    for (qtl_i in 1:nrow(trait_qtl_list)) {
      qtl_chr <- trait_qtl_list$Chromosome[qtl_i]
      qtl_marker_l <- trait_qtl_list$LeftmostMarker[qtl_i]
      qtl_marker_r <- trait_qtl_list$RightmostMarker[qtl_i]

      for (gene_i in 1:nrow(trait_gene_list)) {
        if (trait_gene_list$Chromosome[gene_i] != qtl_chr) {
          next
        }
        gene_locus <- trait_gene_list$Base[gene_i]
        if (gene_locus < qtl_marker_l | gene_locus > qtl_marker_r) {
          next
        }
        trait_qtl_list$NumGenes[qtl_i] <- trait_qtl_list$NumGenes[qtl_i] + 1
        trait_qtl_list$FoundGeneIDs[qtl_i] <-
          paste0(trait_qtl_list$FoundGeneIDs[qtl_i], " and ", trait_gene_list$ID[gene_i])
      }
    }

    # Remove original 0 in FoundGeneIDs
    trait_qtl_list$FoundGeneIDs<- sapply(
      trait_qtl_list$FoundGeneIDs,
      function(x) gsub("0 and ", "", x)  # TODO make a list instead of string
    )

    return(trait_qtl_list)
  }


#' TODO write docs
#'
#' @param gene_chr TODO
#' @param gene_locus TODO
#' @param qtl_ext_length TODO
#' @param placement_type TODO
#' @param per_marker_likelihood TODO
#' @param sectioned_marker_list TODO
#' @param chromosome_size TODO
#'
#' @return TODO
#' @export
GeneFoundLikelihood <- function(
  gene_chr, gene_locus,
  qtl_ext_length, placement_type,
  per_marker_likelihood, sectioned_marker_list,
  chromosome_size
  ) {
  markers_on_chromosome <- sectioned_marker_list[[gene_chr]]
  chr_marker_positions <- markers_on_chromosome$Base

  range_l <- gene_locus - qtl_ext_length
  range_r <- gene_locus + qtl_ext_length

  reachable_markers_on_l <-
    markers_on_chromosome[which(
      chr_marker_positions <= gene_locus &
        chr_marker_positions >= range_l &
        chr_marker_positions + qtl_ext_length <= chromosome_size$RightmostMarker[gene_chr]
    ),]
  reachable_markers_on_r <-
    markers_on_chromosome[which(
      # It's not >= so we don't double-count locus
      chr_marker_positions > gene_locus &
        chr_marker_positions <= range_r &
        chr_marker_positions - qtl_ext_length >= chromosome_size$LeftmostMarker[gene_chr]
    ),]
  if (placement_type == "centered") {
    # check that the other end of the QTL also won't hang off the end of the chromosome
    reachable_markers_on_l <- reachable_markers_on_l[
      which(reachable_markers_on_l$Base - qtl_ext_length >= chromosome_size$LeftmostMarker[gene_chr]),
      ]
    reachable_markers_on_r <- reachable_markers_on_r[
      which(reachable_markers_on_r$Base + qtl_ext_length <= chromosome_size$RightmostMarker[gene_chr]),
      ]
  }
  num_markers_avail <- nrow(reachable_markers_on_l) + nrow(reachable_markers_on_r)
  gene_found_likelihood <- num_markers_avail * per_marker_likelihood

  # Sometimes one of the multiplicands is of type 'levels' for some reason and then gene_found_likelihood is NA
  #jfc I hate levels
  #I think I had an as.numeric() in here somewhere to resolve that
  if (length(gene_found_likelihood) == 0) {
    gene_found_likelihood <- 0
  }

  return(gene_found_likelihood)
}



#' Validate QTL using a list of previously known genes
#'
#' This function determines the number of genes expected to be identified by a QTL
#' of a particular length, given the distribution of genes and mapping markers in
#' the genome. Relies on the Gene_List_Groomer, GeneCounter, MarkerSectioner,
#' and QTL_Placement_Probabilities functions.
#'
#' TODO fix doc
#'
#' @param qtl_list A 7 column matrix containing the chromosome
#' on which the QTL was originally found ('Chromosome'), the Leftmost
#'  Confidence interval's location,the Rightmost Confidence interval's location,
#'  the trait for which the QTL were discovered, the treatment that was applied,
#'  the mapping method, the overall type of experiment, and the lengths of the QTL.
#' @param num_repetitions An integer defining the number of simulations to run.
#' @param gene_list A 4 column matrix containing a list of known genes.
#' Columns should include "GeneID', "Trait", "Chromosome", and "Locus".
#' @param placement_type Either 'extension' or 'centered'. This depends on the
#' mapping style of the original QTL experiment.
#' @param trait The trait for which you would like to analyze your QTL.
#' @param marker_list A 3 column matrix containing the marker loci
#' (ID, Chromosome, and Base)).
#' @param chromosome_size A matrix containing the number of each chromosome,
#'  their total lengths, and the first and last markers used for the original
#'  mapping experiment
#' @param whole_genome_gene_dist A 4 column matrix of all genes in the
#'  genome. Include Chromosome, GeneStart, GeneEnd, and GeneMiddle. The SPQValidate
#'  function uses GeneMiddle to simulate gene placement.
#' @param simulation_env A new environment in which a dataframe containing the
#' results of each simulation is stored. This dataframe can be accessed as
#' Your_Environment$Simulation_Dataframe.
#' @param progress_bar If true, show a progress bar while running.
#' @return A matrix containing the lengths of QTL (as an identifier), the true observed
#' value for identified gene number, the mean value resulting from the simulations,
#' and the upper and lower 95% confidence intervals for the simulation distributions.
#' @export
SPQValidate <- function(qtl_list,
                        trait,
                        num_repetitions,
                        placement_type,
                        gene_list,
                        marker_list,
                        whole_genome_gene_dist,
                        chromosome_size,
                        simulation_env,
                        progress_bar=TRUE) {

  ######Checking inputs#######

  # Split qtl_list into 2 dfs - metadata only needed for grouping
  qtl_metadata <- validateDf(qtl_list, list(
    c("Trait", "character"),
    c("Treatment", "character"),
    c("Method", "character"),
    c("ExptType", "character")
  ))
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

  whole_genome_gene_dist <- validateDf(whole_genome_gene_dist, list(  # rename as whole_genome_gene_list #eh it's a df though  # wait gene_list, marker_list etc aren't dfs?? WHAT IS R EVEN
    c("Chromosome", "integer"),
    c("GeneMiddle", "integer")  # see above - rename to locus?
  ))

  # TODO check that loci are in bp and not in cM?

  if (typeof(trait) != 'character') {
    stop("Trait input must be of type 'character'.")
  }
  trait <- tolower(trait)
  gene_list$Trait <- tolower(gene_list$Trait)
  qtl_list$Trait <- tolower(qtl_list$Trait)

  placement_types <- c("extension", "centered")
  if (!(placement_type %in% placement_types)) {
    stop(
      sprintf(
        "Unrecognized placement_type: %s. \nMust be one of: \n%s \n",
        placement_type,
        paste(placement_types, collapse = ", ")
      )
    )
  }

  if (typeof(simulation_env) != 'environment') {
    stop(
      "Simulation_Environment must be type 'environment'. Use the function
         new.env() to produce an environment."
    )
  }

  if (!(typeof(num_repetitions) %in% c('integer','double','numeric'))) {
    stop("num_repetitions must be an numeric,integer, or double value (i.e. 1000).")
  }

  # end validation #########

  trait_gene_list <- FilterGeneList(
    trait = trait,
    gene_list = gene_list,
    marker_list = marker_list,
    drop_tandem = (placement_type != "centered")
  )

  qtl_gene_counts <- CountGenesFound(
    trait = trait,
    qtl_list = qtl_list,
    trait_gene_list = trait_gene_list,
    marker_list = marker_list
  )
  if (nrow(qtl_gene_counts) <= 1) {
    # TODO Why does everything become "double" if there's only one row?
    qtl_gene_counts$Chromosome <- as.integer(qtl_gene_counts$Chromosome)
    qtl_gene_counts$Length <- as.integer(qtl_gene_counts$Length)
    qtl_gene_counts$LeftmostMarker <- as.integer(qtl_gene_counts$LeftmostMarker)
    qtl_gene_counts$RightmostMarker <- as.integer(qtl_gene_counts$RightmostMarker)
  }


  num_qtl <- nrow(qtl_gene_counts)
  if (sum(qtl_gene_counts$Length) == 0) {
    qtl_gene_counts$Length <-
      qtl_gene_counts$RightmostMarker - qtl_gene_counts$LeftmostMarker
  }
  num_genes <- nrow(trait_gene_list)
  num_chromosomes <- nrow(chromosome_size)
  output <-
    as.data.frame(matrix(
      nrow = num_qtl,
      ncol = num_repetitions,
      data = 0
    ))
  if (progress_bar) {
    pb <- utils::txtProgressBar(0, 1, style = 3)
  }


  sectioned_marker_list <-
    SectionMarkers(marker_list = marker_list,
                    num_chromosomes = num_chromosomes)

  qtl_probs <- QTLPlacementProbabilities(
    placement_type = placement_type,
    qtl_list = qtl_gene_counts,
    sectioned_marker_list = sectioned_marker_list,
    chromosome_size = chromosome_size
  )

  for (rep_i in 1:num_repetitions) {
    random_gene_list <- dplyr::sample_n(whole_genome_gene_dist, num_genes, replace=FALSE)

    for (qtl_i in 1:num_qtl) {
      random_gene_list$EGN <- 0

      if (placement_type == 'centered') {
        qtl_ext_length = qtl_gene_counts$Length[qtl_i] / 2
      } else if (placement_type == 'extension') {
        qtl_ext_length = qtl_gene_counts$Length[qtl_i]
      }

      per_marker_likelihood <- qtl_probs[qtl_i]

      for (gene_i in 1:num_genes) {
        sampled_gene_locus <- random_gene_list$GeneMiddle[gene_i]
        sampled_chr_number <- random_gene_list$Chromosome[gene_i]

        gene_found_likelihood <- GeneFoundLikelihood(sampled_chr_number, sampled_gene_locus,
                                                     qtl_ext_length, placement_type,
                                                     per_marker_likelihood, sectioned_marker_list,
                                                     chromosome_size)

        random_gene_list$EGN[gene_i] <-
          random_gene_list$EGN[gene_i] + gene_found_likelihood
      }

      sim_egn <- sum(random_gene_list$EGN)
      output[qtl_i, rep_i] <- sim_egn

      if (progress_bar) {
        utils::setTxtProgressBar(pb, qtl_i + ((rep_i-1)*num_qtl) / (num_qtl * num_repetitions))
      }
    }
  }

  colnames(output) <- paste0("Simulation_Round_", 1:num_repetitions)

  simulation_env$SimulationDataFrame <- output

  std_devs <- apply(output, MARGIN = 1, FUN = stats::sd)  # 1 means row-wise
  upper <- c()
  center <- c()
  lower <- c()

  for (qtl_i in 1:nrow(qtl_gene_counts)) {
    z_value <-
      stats::qnorm(.025,
            lower.tail = FALSE)
    mu <- rowMeans(output)[qtl_i]

    center <- c(center, mu)
    upper <- c(upper, mu + z_value * std_devs[qtl_i])
    lower <- c(lower, mu - z_value * std_devs[qtl_i])
  }

  #dealing with multiple testing
  #gotta do the triathalon math
  #as in https://stats.stackexchange.com/questions/223924/how-to-add-up-partial-confidence-intervals-to-create-a-total-confidence-interval?fbclid=IwAR0mz3ypdVGQjopYytJnfrXA6o50MJyZ0OAxnUHcwE7kBRKbGqGiEUY6mSY
  upper_sum_of_CIstoSum<-c()
  lower_sum_of_CIstoSum<-c()

  # group by ('Trait','Treatment',"Method","ExptType") by hand bc plyr wasn't cooperating
  metadata_qtl_list<-qtl_metadata[,c('Trait','Treatment',"Method","ExptType")]
  if (! all(rownames(metadata_qtl_list) == rownames(qtl_gene_counts))) {
    stop("Probably runtime error: Rownames for metadata_qtl_list don't match qtl_gene_counts.")
  }
  groups_of_qtl<-unique(metadata_qtl_list)
  for(grouping_i in 1:nrow(groups_of_qtl)){
    tr <- groups_of_qtl$Treatment[grouping_i]
    meth <- groups_of_qtl$Method[grouping_i]
    ET <- groups_of_qtl$ExptType[grouping_i]
    tra <- groups_of_qtl$Trait[grouping_i]

    indices_for_grouping <- rownames(metadata_qtl_list[metadata_qtl_list$Treatment==tr &
                                                       metadata_qtl_list$Trait==tra &
                                                       metadata_qtl_list$ExptType==ET &
                                                       metadata_qtl_list$Method==meth,])
    qtl_gene_counts[indices_for_grouping,'QTLGroup'] <- paste0('Group',grouping_i)
  }

  for(grouping_i in unique(qtl_gene_counts$QTLGroup)){
    CIstoSum_indices<-which(qtl_gene_counts$QTLGroup==grouping_i)

    dist_ctr <- center[CIstoSum_indices]
    dist_radius<-upper[CIstoSum_indices] - dist_ctr
    sqr_rad<-dist_radius^2

    adjusted_center<-mean(dist_ctr)
    adjusted_radius<-sqrt(sum(sqr_rad))

    upper_sum_of_CIstoSum[CIstoSum_indices]<- adjusted_center + adjusted_radius
    lower_sum_of_CIstoSum[CIstoSum_indices]<- adjusted_center - adjusted_radius
  }


  conf_ints <-
    as.data.frame(as.matrix(
      cbind(
        qtl_gene_counts$Length,
        qtl_gene_counts$NumGenes,
        rowMeans(output),
        std_devs,
        lower,
        upper,
        lower_sum_of_CIstoSum,
        upper_sum_of_CIstoSum

     )
  ))

  colnames(conf_ints) <-
    c("QTL",
      "Observed Value",
      "Mean",
      "SEM",
      "Lower 95% CI",
      "Upper 95% CI",
      "Additive Lower 95% CI",
      "Additive Upper 95% CI")
  return(conf_ints)
}

