# Naming case conventions:
#  - exported functions: UpperCamelCase
#  - private functions: lowerCamelCase
#  - variables: snake_case, no caps
#  - constants (rare): SNAKE_CASE, all caps
#  - Expected dataframe cols: UpperCamelCase

validateDf <- function(df, req_cols, return_stripped = True) {
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
  }
}


# TODO add back docs
SectionMarkers <- function(marker_list, num_chromosomes) {
  validateDf(marker_list, list(
    c("Id", "character"),
    c("Chromosome", "integer"),
    c("Base", "integer")
  ))

  if (typeof(num_chromosomes) != 'integer') {
    stop("Input num_chromosomes must be of type 'integer'.")
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
    marker_df <- MarkerList[which(MarkerList$Chromosome == chromosome), ]
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

    validateDf(gene_list, list(
      c("ID", "character"),
      c("Trait", "character"),
      c("Chromosome", "integer"),
      c("Base", "integer")
    ))

    validateDf(marker_list, list(
      c("ID", "character"),
      c("Chromosome", "integer"),
      c("Base", "integer")
    ))

    trait_gene_list <- gene_list[which(gene_list$Trait==trait), ]
    trait_gene_list <- unique(trait_gene_list)

    if (length(trait_gene_list$Trait) == 0) {
      stop('No genes for this trait.')
    }
    if (placement_type == 'centered' |
        length(trait_gene_list$Trait) == 1) {
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
    }
    return(unique(trait_gene_list[grep(1, presence_absence), ]))
  }

QTLPlacementProbabilities <-
  function(qtl_list,
           placement_type,
           sectioned_marker_list,
           chromosome_size) {

    validateDf(qtl_list, list(
      c("Chromosome", "integer"),
      c("LeftmostMarker", "integer"),
      c("RightmostMarker", "integer"),
      c("Trait", "character"),
      c("Treatment", "character"),
      c("Method", "character"),
      c("ExptType", "character"),
      c("Length", "integer")
    ))

    validateDf(chromosome_size, list(
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
    chromosome_size <- chromosome_size[order(as.numeric(chromosome_size$Number)), ]
    num_chromosomes <- as.numeric(length(chromosome_size[, 1]))

    if (length(qtl_list$Length) == 0 | sum(qtl_list$Length) == 0) {
      qtl_list$Length <-
        qtl_list$RightmostMarker - qtl_list$LeftmostMarker

    }

    poss_positions <-
      as.data.frame(matrix(
        nrow = length(qtl_list$Chromosome),
        ncol = 1 + num_chromosomes
      ))
    colnames(poss_positions) <- c("QTLLength", 1:num_chromosomes)
    poss_positions$QTLLength <- qtl_list$Length


    for (qtl_i in 1:length(qtl_list$Length)) {

      if (placement_type == 'centered') {
        qtl_ext_length = qtl_list$Length[qtl_i] / 2
      } else if (placement_type == 'extension') {
        qtl_ext_length = qtl_list$Length[qtl_i]
      }

      for (chr_i in  1:num_chromosomes) {
        first_marker <-
          chromosome_size$LeftmostMarker[chr_i])
        last_marker <-
          chromosome_size$RightmostMarker[chr_i])

        first_avail_marker <- first_marker + qtl_ext_length
        last_avail_marker <- last_marker - qtl_ext_length

        chr_markers <- sectioned_list[chr_i]

        remaining_markers_lr <-
          chr_markers[which(chr_markers$Base < last_avail_marker), ]
        remaining_markers_rl <-
          chr_markers[which(chr_markers$Base > first_avail_marker), ]

        remaining_markers <-
          rbind(remaining_markers_lr, remaining_markers_rl)
        length_remaining_markers <- length(remaining_markers[,1])

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

    validateDf(qtl_list, list(
      c("Chromosome", "integer"),
      c("LeftmostMarker", "integer"),
      c("RightmostMarker", "integer"),
      c("Trait", "character"),
      c("Treatment", "character"),
      c("Method", "character"),
      c("ExptType", "character"),
      c("Length", "integer")
    ))

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

    #Everything except QTL_with_Metadata gets checked here:
    trait_gene_list <-
      FilterGeneList(
        trait = trait,
        placement_type = placement_type,
        gene_list = gene_list,
        marker_list = marker_list
      )

    if (length(trait_gene_list) <= 1) {
      return(trait_gene_list)

    }

    # Checking  QTL_with_Metadata
    #####
    if (typeof(QTL_with_Metadata$Leftmost_Marker[1]) == 'double') {
      QTL_with_Metadata$Leftmost_Marker <-
        as.integer(QTL_with_Metadata$Leftmost_Marker)
    }
    if (typeof(QTL_with_Metadata$Leftmost_Marker[1]) != 'integer' &
        typeof(QTL_with_Metadata$Leftmost_Marker[1]) != 'numeric') {
      warning(
        "QTL_with_Metadata Leftmost_Marker column contains neither integer nor numeric types."
      )
    }
    if (typeof(QTL_with_Metadata$Rightmost_Marker[1]) == 'double') {
      QTL_with_Metadata$Rightmost_Marker <-
        as.integer(QTL_with_Metadata$Rightmost_Marker)
    }
    if (typeof(QTL_with_Metadata$Rightmost_Marker[1]) != 'integer' &
        typeof(QTL_with_Metadata$Rightmost_Marker[1]) != 'numeric') {
      warning(
        "QTL_with_Metadata Rightmost_Marker column contains neither integer nor numeric types."
      )
    }
    if (typeof(QTL_with_Metadata$Chromosome[1]) == 'double') {
      QTL_with_Metadata$Chromosome <-
        as.integer(QTL_with_Metadata$Chromosome)
    }
    if (typeof(QTL_with_Metadata$Chromosome[1]) != 'integer' &
        typeof(QTL_with_Metadata$Chromosome[1]) != 'numeric') {
      warning("QTL_with_Metadata Chromosome column contains neither integer nor numeric types.")
    }

    if (mean(nchar(QTL_with_Metadata$Chromosome)) > 4) {
      warning(
        "QTL_with_Metadata Chromosome number is unusually large. Column may represent confidence interval instead."
      )
    }

    if (mean(nchar(QTL_with_Metadata$Leftmost_Marker)) < 4) {
      warning(
        "QTL_with_Metadata Leftmost_Marker number is unusually small. Column may represent Chromosome instead."
      )
    }

    if (mean(nchar(QTL_with_Metadata$Rightmost_Marker)) < 4) {
      warning(
        "QTL_with_Metadata Rightmost_Marker number is unusually small. Column may represent Chromosome instead."
      )
    }
    if (mean(QTL_with_Metadata$Leftmost_Marker) > mean(QTL_with_Metadata$Rightmost_Marker)) {
      warning(
        "QTL_with_Metadata Leftmost_Marker is larger than Rightmost_Marker. Check column order."
      )
    }
    if (typeof(QTL_with_Metadata$QTL_Type) != 'character') {
      warning(
        "QTL_with_Metadata$QTL_Type should contain character values to ensure appropriate column use."
      )
    }
    if (typeof(QTL_with_Metadata$Trait) != 'character') {
      warning(
        "QTL_with_Metadata$Trait should contain character values to ensure appropriate column use."
      )
    }
    if (typeof(QTL_with_Metadata$Treatment) != 'character') {
      warning(
        "QTL_with_Metadata$Treatment should contain character values to ensure appropriate column use."
      )
    }


    #####


    if (as.numeric(as.character(length(QTL_with_Metadata[, 1]))) > 1) {
      TraitsandTreatments <-
        subset(QTL_with_Metadata, select = c(Trait, Treatment))
      sepTreatments_QTL_Count <-
        plyr::ddply(
          TraitsandTreatments,
          .(
            TraitsandTreatments$Trait,
            TraitsandTreatments$Treatment
          ),
          nrow
        )
      colnames(sepTreatments_QTL_Count) <-
        c("Trait", 'Treatment', 'Frequency')



      QTLData <-
        QTL_with_Metadata[grep(Trait, TraitsandTreatments$Trait), ]
      if (length(QTLData$Chromosome) == 0) {
        return("No QTL found for this trait.")
      }

      QTLData$N_QTL <- c(0)
      count <- c()

      for (identified in 1:length(QTLData$Chromosome)) {
        sep_Treat_QTL_4_u <-
          sepTreatments_QTL_Count$Frequency[which(
            sepTreatments_QTL_Count$Trait == QTLData$Trait[identified] &
              sepTreatments_QTL_Count$Treatment ==
              QTLData$Treatment[identified]
          )]
        count <- c(count, sep_Treat_QTL_4_u)
      }

      QTLData$N_QTL <- count
      QTLData$N_Genes <- c(0)

      identified_genes <- data.frame(matrix(ncol = length(QTLData)))
      colnames(identified_genes) <- colnames(QTLData)


      for (i in 1:length(QTLData$Chromosome)) {
        QTLchromosome <- as.numeric(as.character(QTLData$Chromosome[i]))
        QTLLCI <- as.numeric(as.character(QTLData$Leftmost_Marker[i]))
        QTLRCI <- as.numeric(as.character(QTLData$Rightmost_Marker[i]))

        for (Gene in 1:length(trait_gene_list$GeneID)) {
          if (as.numeric(as.character(trait_gene_list$Chromosome[Gene])) == as.numeric(as.character(QTLchromosome))) {
            GeneStartSite <- as.numeric(trait_gene_list$Locus[Gene])
            if (GeneStartSite >= QTLLCI & GeneStartSite <= QTLRCI) {
              QTLData$N_Genes[i] <-
                paste0(QTLData$N_Genes[i], " and ", trait_gene_list$GeneID[Gene])
              identified_genes <- rbind(identified_genes, QTLData[i, ])

            }
          }
        }
      }
      identified_genes <- identified_genes[-1, ]
      identified_genes$Length <-
        as.numeric(as.character(identified_genes$Rightmost_Marker)) - as.numeric(as.character(identified_genes$Leftmost_Marker))


      identified_genes$N_Genes <-
        stringr::str_count(identified_genes$N_Genes, 'and')
      identified_genes <-
        identified_genes[, c(
          "Chromosome",
          "Leftmost_Marker",
          "Rightmost_Marker",
          "Trait",
          'Treatment',
          "Length",
          "QTL_Type",
          "N_Genes",
          'N_QTL'
        )]
      colnames(identified_genes)[9] <- 'Number_Trait_QTL'
      colnames(QTLData)[7] <- 'Number_Trait_QTL'
      identified_genes2 <- identified_genes

      if (length(identified_genes[, 1]) > 1) {
        for (possible_Duplicates in 1:(length(identified_genes[, 1]) - 1)) {
          if (all(identified_genes[possible_Duplicates, c(1:6)] == identified_genes[possible_Duplicates +
                                                                                    1, c(1:6)])) {
            identified_genes2[possible_Duplicates, ] <- 0
          }
        }
        if (length(which(identified_genes2$N_Genes == 0)) > 0) {
          CountedIdentifiedGenes <-
            identified_genes2[-which(identified_genes2$N_Genes == 0), ]
        } else{
          CountedIdentifiedGenes <- identified_genes2
        }

        CountedIdentifiedGenesOutput <- CountedIdentifiedGenes
        NoGenes <-
          dplyr::setdiff(QTLData[1:7], CountedIdentifiedGenes[1:7])

        if (length(NoGenes$Chromosome) > 0) {
          NoGenes$Length <-
            as.numeric(as.character(NoGenes$Rightmost_Marker)) - as.numeric(as.character(NoGenes$Leftmost_Marker))
          NoGenes$N_Genes <- 0
          NoGenesOutput <- NoGenes
          colnames(NoGenesOutput) <- c(
            "Chromosome",
            "Leftmost_Marker",
            "Rightmost_Marker",
            "Trait",
            "Treatment",
            "QTL_Type",
            "Number_Trait_QTL",
            "Length",
            "N_Genes"
          )
          colnames(CountedIdentifiedGenesOutput) <-
            c(
              "Chromosome",
              "Leftmost_Marker",
              "Rightmost_Marker",
              "Trait",
              "Treatment",
              "QTL_Type",
              "Number_Trait_QTL",
              "N_Genes",
              "Length"
            )


          CountedIdentifiedGenesOutput2 <-
            rbind(CountedIdentifiedGenesOutput, NoGenesOutput)
        } else{
          CountedIdentifiedGenesOutput2 <- CountedIdentifiedGenesOutput
        }

        return(CountedIdentifiedGenesOutput2)
      }


      if (length(identified_genes[, 1]) == 1) {
        CountedIdentifiedGenesOutput <- identified_genes[, c(1:5, 7, 9, 8, 6)]
        toCompareID <- CountedIdentifiedGenesOutput[, c(1:7)]
        toCompareQTLData <- QTLData[, c(1:7)]
        NoGenes <- dplyr::setdiff(toCompareQTLData, toCompareID)


        if (length(NoGenes$Chromosome) > 0) {
          NoGenes$N_Genes <- c(0)
          NoGenes$Length <-
            as.numeric(as.character(NoGenes$Rightmost_Marker)) - as.numeric(as.character(NoGenes$Leftmost_Marker))
          NoGenesOutput <- NoGenes

          colnames(NoGenesOutput) <- c(
            "Chromosome",
            "Leftmost_Marker",
            "Rightmost_Marker",
            "Trait",
            "Treatment",
            "QTL_Type",
            "Number_Trait_QTL",
            "N_Genes",
            "Length"
          )
          colnames(CountedIdentifiedGenesOutput) <-
            c(
              "Chromosome",
              "Leftmost_Marker",
              "Rightmost_Marker",
              "Trait",
              "Treatment",
              "QTL_Type",
              "Number_Trait_QTL",
              "N_Genes",
              "Length"
            )


          CountedIdentifiedGenesOutput2 <-
            rbind(CountedIdentifiedGenesOutput, NoGenesOutput)
          CountedIdentifiedGenesOutput2 <-
            CountedIdentifiedGenesOutput2[, c(1:6, 9, 7:8)]
          return(CountedIdentifiedGenesOutput2)
        } else{
          CountedIdentifiedGenesOutput2 <- CountedIdentifiedGenesOutput
          return(CountedIdentifiedGenesOutput2)
        }
        return(CountedIdentifiedGenesOutput2)
      } else{
        return("No identified genes for this trait")
      }
    }
    else{
      i = 1
      QTL_with_Metadata$Length <-
        as.numeric(as.character(QTL_with_Metadata$Rightmost_Marker)) - as.numeric(as.character(QTL_with_Metadata$Leftmost_Marker))
      QTL_with_Metadata$N_Genes <-
        as.character(QTL_with_Metadata$N_Genes)
      QTLchromosome <-
        as.numeric(as.character(QTL_with_Metadata$Chromosome[i]))
      QTLLCI <-
        as.numeric(as.character(QTL_with_Metadata$Leftmost_Marker[i]))
      QTLRCI <-
        as.numeric(as.character(QTL_with_Metadata$Rightmost_Marker[i]))

      identified_genes <-
        data.frame(matrix(ncol = length(QTL_with_Metadata)))
      colnames(identified_genes) <- colnames(QTL_with_Metadata)

      for (Gene in 1:length(gene_list$GeneID)) {
        if (as.numeric(as.character(trait_gene_list$Chromosome[Gene])) == as.numeric(as.character(QTLchromosome))) {
          GeneStartSite <- as.numeric(trait_gene_list$Locus[Gene])
          if (GeneStartSite >= QTLLCI & GeneStartSite <= QTLRCI) {
            QTL_with_Metadata$N_Genes[i] <-
              paste0(QTL_with_Metadata$N_Genes[i],
                     " and ",
                     trait_gene_list$GeneID[Gene])
            identified_genes <-
              rbind(identified_genes, QTL_with_Metadata[i, ])
          }
        }
      }
      identified_genes <- identified_genes[-1, ]
      identified_genes$N_Genes <-
        stringr::str_count(identified_genes$N_Genes, 'and')
      identified_genes2 <- identified_genes
      if (length(as.numeric(as.character(identified_genes[, 1]))) > 1) {
        for (possible_Duplicates in 1:(length(identified_genes[, 1]) - 1)) {
          if (all(identified_genes[possible_Duplicates, c(1:6)] == identified_genes[possible_Duplicates +
                                                                                    1, c(1:6)])) {
            identified_genes2[possible_Duplicates, c(1:7)] <- 0
          }
        }

      }
      if (length(which(identified_genes2$N_Genes == 0)) > 0) {
        CountedIdentifiedGenes <-
          identified_genes2[-which(identified_genes2$N_Genes == 0), ]
      } else{
        CountedIdentifiedGenes <-
          identified_genes2[, c(
            "Chromosome",
            "Leftmost_Marker",
            "Rightmost_Marker",
            "Trait",
            'Treatment',
            "Length",
            "QTL_Type",
            "N_Genes",
            'N_QTL'
          )]
        colnames(CountedIdentifiedGenes)[9] <- 'Number_Trait_QTL'
      }
      CountedIdentifiedGenes$N_QTL <- c(1)

      colnames(CountedIdentifiedGenes) <-
        c(
          "Chromosome",
          "Leftmost_Marker",
          "Rightmost_Marker",
          "Trait",
          'Treatment',
          "Length",
          "QTL_Type",
          "N_Genes",
          "Number_Trait_QTL"
        )
      return(CountedIdentifiedGenes)
    }
  }



# TODO add back documentation
SPQValidate <- function(qtl_with_metadata,
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

  validateDf(qtl_list, list(
    c("Chromosome", "integer"),
    c("LeftmostMarker", "integer"),
    c("RightmostMarker", "integer"),
    c("Trait", "character"),
    c("Treatment", "character"),
    c("Method", "character"),
    c("ExptType", "character"),
    c("Length", "integer")
  ))

  validateDf(gene_list, list(
    c("GeneId", "character"),
    c("Trait", "character"),
    c("Chromosome", "integer"),
    c("Base", "integer")
  ))

  validateDf(marker_list, list(
    c("Id", "character"),
    c("Chromosome", "integer"),
    c("Base", "integer")
  ))

  validateDf(chromosome_size, list(
    c("Chromosome", "integer"),
    c("LeftmostMarker", "integer"),
    c("RightmostMarker", "integer"),
    c("Length", "integer")
  ))

  validateDf(whole_genome_gene_dist, list(
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

  if (typeof(num_repetitions) != 'integer') {
    stop("num_repetitions must be an integer value (i.e. 1000).")
  }
  #####
  #####
  #####

  qtl_of_interest <- GeneCounter(
    trait = trait,
    placement_type = placement_type,
    qtl_list = qtl_list,
    gene_list = gene_list,
    marker_list = marker_list
  )
  if (length(qtl_of_interest) <= 1) {
    return(qtl_of_interest)
  }


  gene_list <- FilterGeneList(
    trait = trait,
    placement_type = placement_type,
    gene_list = gene_list,
    marker_list = marker_list
  )
  gene_list <-
    subset(gene_list, select = c(GeneID, Chromosome, Locus))


  qtl_of_interest <-
    qtl_of_interest[grep(trait, qtl_of_interest$Trait),]
  num_qtl <- length(qtl_of_interest[, 1])
  if (sum(qtl_of_interest$Length) == 0) {
    qtl_of_interest$Length <-
      qtl_of_interest$RightmostMarker - qtl_of_interest$LeftmostMarker
  }
  num_genes <- length(gene_list[, 1])
  num_chromosomes <- length(chromosome_size[, 1])
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
    random_gene_list <- sample_n(whole_genome_gene_dist, num_genes, replace=FALSE)
    random_gene_list$EGN <- 0

    for (qtl_i in 1:num_qtl) {
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
        markers_on_chromosome <- sectioned_marker_list[chr_number]
        chr_marker_positions <- markers_on_chromosome$Base

        range_l <- sampled_gene_locus - qtl_ext_length
        range_r <- sampled_gene_locus + qtl_ext_length

        markers_l <-
          markers_on_chromosome[which(
            chr_marker_positions + qtl_ext_length <= chromosome_size$RightmostMarker[sampled_chr_number] &
              chr_marker_positions <= sampled_gene_locus &
              chr_marker_positions >= range_l
          ),]
        length_markers_l <- length(markers_l[,1])

        markers_r <-
          markers_on_chromosome[which(
            chr_marker_positions - qtl_ext_length >= chromosome_size$LeftmostMarker[sampled_chr_number] &
              # It's not >= so we don't double-count locus
              chr_marker_positions > as.numeric(sampled_gene_locus) &
              chr_marker_positions <= range_r
          ),]
        length_markers_r <- length(markers_r[,1])

        num_markers_avail <- length_markers_l + length_markers_r
        num_markers_expect <- num_markers_avail * known_likelihood
        if (length(num_markers_expect) == 0) {  # TODO is this necessary?
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
  output$QTLLength <- qtl_of_interest$Length
  output$ObservedValue <- qtl_of_interest$N_Genes
  simulation_env$SimulationDataFrame <- output

  std_devs <- apply(output, MARGIN = 1, FUN = sd)
  upper <- c()
  lower <- c()

  for (qtl_i in 1:length(qtl_of_interest$Length)) {
    z_value <-
      qnorm(.025 / qtl_of_interest$Number_Trait_QTL[QTL],
            lower.tail = FALSE)
    mu <- rowMeans(output)[QTL]
    upper <- c(upper, mu + z_value * std_devs[qtl_i])
    lower <- c(lower, mu - z_value * std_devs[qtl_i])
  }

  conf_ints <-
    as.data.frame(as.matrix(
      cbind(
        qtl_of_interest$Length,
        qtl_of_interest$N_Genes,
        rowMeans(output),
        std_devs,
        lower,
        upper
     )
  ))
  colnames(conf_ints) <-
    c("QTL",
      "Observed Value",
      "Mean",
      "SEM",
      "Lower 95% CI",
      "Upper 95% CI")
  return(conf_ints)
}
