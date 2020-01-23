
validateDf <- function(df, req_cols, return_stripped = True) {
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


FilterGeneList <-
  function(trait,
           placement_type,
           gene_list,
           marker_list) {
    if (typeof(trait) != 'character') {
      stop("Trait input must be of type 'character'.")
    }
    trait <- tolower(trait)
    gene_list$Trait <- tolower(gene_list$Trait)

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
      stop('No genes for this trait.')
    }
    if (placement_type == 'centered' | nrow(trait_gene_list) == 1) {
      return(trait_gene_list)
    }
    if (placement_type == 'extension') {
      marker_list$Trait <- "Marker"
      loci_list <- rbind(marker_list, trait_gene_list)
      loci_list <-
        loci_list[order(loci_list$Chromosome, loci_list$Base), ]
      gene_spots<-which(loci_list$Trait==trait)
      tandem_arrays <- c()
      # TODO ask whats happening here, leaving rest of func alone bc i dont know
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
      # TODO orig return was here, moved it bc assuming that was unintentional?
      #probs ok, testing isn't showing problems
    }
    return(unique(trait_gene_list[grep(1, presence_absence), ]))
  }

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
      c("Treatment", "character"),
      c("Method", "character"),
      c("ExptType", "character"),
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

        first_avail_marker <- first_marker + qtl_ext_length
        last_avail_marker <- last_marker - qtl_ext_length

        #this line is the bug - can't just pick a df out of a list like this
        chr_markers <-sectioned_marker_list[chr_i]
        chr_markers<-chr_markers[[1]]

        remaining_markers_lr <-
          chr_markers[which(chr_markers$Base < last_avail_marker), ]
        remaining_markers_rl <-
          chr_markers[which(chr_markers$Base > first_avail_marker), ]

        remaining_markers <-
          rbind(remaining_markers_lr, remaining_markers_rl)
        length_remaining_markers <- nrow(remaining_markers)

        poss_positions[qtl_i, chr_i + 1] <-
          length_remaining_markers

        probs <- 1 / rowSums(poss_positions[, 2:10])
      }
    }

    return(probs)
  }


CountGenes <-
  function(qtl_list,
           trait,
           placement_type,
           gene_list,
           marker_list) {

    qtl_list  <- validateDf(qtl_list, list(
      c("Chromosome", "integer"),
      c("LeftmostMarker", "integer"),
      c("RightmostMarker", "integer"),
      c("Trait", "character"),
      c("Treatment", "character"),
      c("Method", "character"),
      c("ExptType", "character"),
      c("Length", "integer")
    ))

    if  (typeof(trait) != 'character') {
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

    trait_gene_list <-
      FilterGeneList(
        trait = trait,
        placement_type = placement_type,
        gene_list = gene_list,
        marker_list = marker_list
      )

    if (nrow(trait_gene_list) <= 1) {
      return(trait_gene_list)
    }

    trait_qtl_list <-
      qtl_list[which(qtl_list$Trait==trait), ]
    if (nrow(trait_qtl_list) == 0) {
      stop("No QTL found for this trait.")
    }
    
    
    metadata_qtl_list<-trait_qtl_list[,c('Trait','Treatment',"Method","ExptType")]
    groups_of_qtl<-unique(metadata_qtl_list)
   
    
    for(grouping_i in 1:nrow(groups_of_qtl)){
      tr<-groups_of_qtl$Treatment[grouping_i]
      meth<-groups_of_qtl$Method[grouping_i]
      ET<-groups_of_qtl$ExptType[grouping_i]
      tra<-groups_of_qtl$Trait[grouping_i]
   
     indices_for_grouping<-rownames(trait_qtl_list[trait_qtl_list$Treatment==tr & 
                                                     trait_qtl_list$Trait==tra &
                                                     trait_qtl_list$ExptType==ET &
                                                     trait_qtl_list$Method==meth,])
     trait_qtl_list[indices_for_grouping,'NumQTL']<-paste0('Group',grouping_i)
   
    }
    
    

    trait_qtl_list$NumGenes <- 0
    trait_qtl_list$FoundGeneIDs <- 0

    identified_genes <- data.frame(matrix(ncol = ncol(trait_qtl_list)))
    colnames(identified_genes) <- colnames(trait_qtl_list)

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
        # TODO I stopped here for now, went a bit farther but
        # got confused by identified_genes and then caught sight of the clock D: oops
        trait_qtl_list$NumGenes[qtl_i] <- trait_qtl_list$NumGenes[qtl_i] + 1
        trait_qtl_list$FoundGeneIDs[qtl_i] <-
          paste0(trait_qtl_list$FoundGeneIDs[qtl_i], " and ", trait_gene_list$ID[gene_i])
        identified_genes <- rbind(identified_genes, trait_qtl_list[qtl_i, ])

      }
    }
    identified_genes <- identified_genes[-1, ]  # remove initial NA row


    identified_genes <-
      identified_genes[, c(
        colnames(trait_qtl_list)
      )]  # double checking - keep
    names(identified_genes)[names(identified_genes) == "NumQTL"] <- 'QTLGroup' # TODO change to rename
    names(trait_qtl_list)[names(trait_qtl_list) == "NumQTL"] <- 'QTLGroup'
    dedup_identified_genes <- identified_genes

    if (nrow(identified_genes) == 0) {
      print(sprintf("No identified genes for trait: %s", trait))
      return(trait_qtl_list)

    } else if (nrow(identified_genes) == 1) {
      dedup_identified_genes <- identified_genes

    } else {
      for (possible_Duplicates in 1:(nrow(identified_genes) - 1)) {
        if (all(
          identified_genes[possible_Duplicates, c(1:6)] ==  # TODO replace with col names
          identified_genes[possible_Duplicates + 1, c(1:6)]
        )) {
          dedup_identified_genes[possible_Duplicates, ] <- 0
        }
      }
      if (length(which(dedup_identified_genes$NumGenes == 0)) > 0) {
        dedup_identified_genes <-
          dedup_identified_genes[-which(dedup_identified_genes$NumGenes == 0), ]
      }
    }

    empty_qtl <-
      dplyr::setdiff(trait_qtl_list, dedup_identified_genes)

    dedup_identified_genes$FoundGeneIDs<- sapply(
      dedup_identified_genes$FoundGeneIDs,
      function(x) gsub("0 and ", "", x)  # TODO make a list instead of string
    )
    if (nrow(empty_qtl) > 0) {
      qtl_gene_counts <-
        rbind(dedup_identified_genes, empty_qtl)
    } else{
      qtl_gene_counts <- dedup_identified_genes
    }

    return(qtl_gene_counts)
  }



# TODO add back documentation
SPQValidate <- function(qtl_list,
                        trait,
                        num_repetitions,
                        placement_type,
                        gene_list,
                        marker_list,
                        whole_genome_gene_dist,
                        chromosome_size,
                        simulation_env) {
  ############################
  ######Checking inputs#######
  ############################

  qtl_list <- validateDf(qtl_list, list(
    c("Chromosome", "integer"),
    c("LeftmostMarker", "integer"),
    c("RightmostMarker", "integer"),
    c("Trait", "character"),
    c("Treatment", "character"),
    c("Method", "character"),
    c("ExptType", "character"),
    c("Length", "integer")
  ))

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

  chromosome_size <- validateDf(chromosome_size, list(
    c("Chromosome", "integer"),
    c("LeftmostMarker", "integer"),
    c("RightmostMarker", "integer"),
    c("Length", "integer")
  ))

  whole_genome_gene_dist <- validateDf(whole_genome_gene_dist, list(
    c("Chromosome", "integer"),
    c("GeneStart", "integer"),
    c("GeneEnd", "integer"),
    c("GeneMiddle", "integer")
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
  #####
  #####
  #####

  qtl_of_interest <- CountGenes(
    trait = trait,
    placement_type = placement_type,
    qtl_list = qtl_list,
    gene_list = gene_list,
    marker_list = marker_list
  )
  if (nrow(qtl_of_interest) <= 1) {
    return(qtl_of_interest)
  }


  gene_list <- FilterGeneList(
    trait = trait,
    placement_type = placement_type,
    gene_list = gene_list,
    marker_list = marker_list
  )
  gene_list <-
    subset(gene_list, select = c(ID, Chromosome, Base))


  #qtl_of_interest <- qtl_of_interest[which(qtl_of_interest$Trait==trait),]
  num_qtl <- nrow(qtl_of_interest)
  if (sum(qtl_of_interest$Length) == 0) {
    qtl_of_interest$Length <-
      qtl_of_interest$RightmostMarker - qtl_of_interest$LeftmostMarker
  }
  num_genes <- nrow(gene_list)
  num_chromosomes <- nrow(chromosome_size)
  output <-
    as.data.frame(matrix(
      nrow = num_qtl,
      ncol = num_repetitions,
      data = 0
    ))
  pb <- txtProgressBar(0, 1, style = 3)

  sectioned_marker_list <-
    SectionMarkers(marker_list = marker_list,
                    num_chromosomes = num_chromosomes)

  qtl_probs <- QTLPlacementProbabilities(
    placement_type = placement_type,
    qtl_list = qtl_of_interest,
    sectioned_marker_list = sectioned_marker_list,
    chromosome_size = chromosome_size
  )

  for (rep_i in 1:num_repetitions) {
    random_gene_list <- dplyr::sample_n(whole_genome_gene_dist, num_genes, replace=FALSE)

    for (qtl_i in 1:num_qtl) {
      random_gene_list$EGN <- 0
      if (placement_type == 'centered') {
        qtl_ext_length = qtl_of_interest$Length[qtl_i] / 2
      } else if (placement_type == 'extension') {
        qtl_ext_length = qtl_of_interest$Length[qtl_i]
      }

      known_likelihood <- qtl_probs[qtl_i]  # TODO call qtl_chr_likelihood?

      for (gene_i in 1:num_genes) {
        sampled_gene_locus <- random_gene_list$GeneMiddle[gene_i]
        sampled_chr_number <- random_gene_list$Chromosome[gene_i]

        # TODO was this instead - keep in case my change doesn't work
        # markers_on_chromosome <-
        #   as.data.frame(as.matrix(eval(as.name(
        #     sectioned_marker_list[chr_number]
        #   ))))
        markers_on_chromosome <- sectioned_marker_list[[sampled_chr_number]]
        chr_marker_positions <- markers_on_chromosome$Base

        range_l <- sampled_gene_locus - qtl_ext_length
        range_r <- sampled_gene_locus + qtl_ext_length

        markers_l <-
          markers_on_chromosome[which(
            chr_marker_positions + qtl_ext_length <= chromosome_size$RightmostMarker[sampled_chr_number] &
              chr_marker_positions <= sampled_gene_locus &
              chr_marker_positions >= range_l
          ),]
        length_markers_l <- nrow(markers_l)

        markers_r <-
          markers_on_chromosome[which(
            chr_marker_positions - qtl_ext_length >= chromosome_size$LeftmostMarker[sampled_chr_number] &
              # It's not >= so we don't double-count locus
              chr_marker_positions > as.numeric(sampled_gene_locus) &
              chr_marker_positions <= range_r
          ),]
        length_markers_r <- nrow(markers_r)

        num_markers_avail <- length_markers_l + length_markers_r
        num_markers_expect <- num_markers_avail * known_likelihood
        # Sometimes one of the multiplicands is of type 'levels' for some reason and then num_markers_expect is NA
        if (length(num_markers_expect) == 0) {
          num_markers_expect <- 0
        }
        random_gene_list$EGN[gene_i] <-
          random_gene_list$EGN[gene_i] + num_markers_expect
      }

      sim_egn <- sum(random_gene_list$EGN)
      output[qtl_i, rep_i] <- sim_egn

      setTxtProgressBar(pb, qtl_i + ((rep_i-1)*num_qtl) / (num_qtl * num_repetitions))
    }
  }

  colnames(output) <- paste0("Simulation_Round_", 1:num_repetitions)

  simulation_env$SimulationDataFrame <- output

  std_devs <- apply(output, MARGIN = 1, FUN = sd)
  upper <- c()
  lower <- c()

  for (qtl_i in 1:nrow(qtl_of_interest)) {
    z_value <-
      qnorm(.025,
            lower.tail = FALSE)
    mu <- rowMeans(output)[qtl_i]
    upper <- c(upper, mu + z_value * std_devs[qtl_i])
    lower <- c(lower, mu - z_value * std_devs[qtl_i])
  }
  #dealing with multiple testing
  for(grouping_i in unique(qtl_of_interest$QTLGroup)){
    CIstoSum_indices<-which(qtl_of_interest$QTLGroup==grouping_i)
    upper_sum_of_CIStoSum<-sum(upper[CIstoSum_indices])
    lower_sum_of_CIStoSum<-sum(lower[CIstoSum_indices])
  }
  

  conf_ints <-
    as.data.frame(as.matrix(
      cbind(
        qtl_of_interest$Length,
        qtl_of_interest$NumGenes,
        rowMeans(output),
        std_devs,
        lower,
        lower_sum_of_CIStoSum,
        upper,
        upper_sum_of_CIStoSum
        
     )
  ))

  colnames(conf_ints) <-
    c("QTL",
      "Observed Value",
      "Mean",
      "SEM",
      "Lower 95% CI",
      "Additive Lower 95% CI",
      "Upper 95% CI",
      "Additive Upper 95% CI")
  return(conf_ints)
}

