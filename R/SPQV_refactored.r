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
# getOnlyCols <- function(df, req_names, req_types=NULL) {
#   if ()
#
# }

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

  validate_df(qtl_list, list(
    c("Chromosome", "integer"),
    c("LeftmostMarker", "integer"),
    c("RightmostMarker", "integer"),
    c("Trait", "character"),
    c("Treatment", "character"),
    c("Method", "character"),
    c("ExptType", "character"),
    c("Length", "integer")
  ))

  validate_df(gene_list, list(
    c("GeneId", "character"),
    c("Trait", "character"),
    c("Chromosome", "integer"),
    c("Base", "integer")
  ))

  validate_df(marker_list, list(
    c("Id", "character"),
    c("Chromosome", "integer"),
    c("Base", "integer")
  ))

  validate_df(chromosome_size, list(
    c("Chromosome", "integer"),
    c("LeftmostMarker", "integer"),
    c("RightmostMarker", "integer"),
    c("Length", "integer")
  ))

  validate_df(whole_genome_gene_dist, list(
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


  gene_list <- GeneListGroomer(
    trait = trait,
    placement_type = placement_type,
    gene_list = gene_list,
    marker_list = marker_list
  )
  gene_list <-
    subset(gene_list, select = c(GeneID, Chromosome, Locus))
  gene_list$EGN <- 0


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
    MarkerSectioner(marker_list = marker_list,
                    num_chromosomes = num_chromosomes)
  qtl_probs <- QTLPlacementProbabilities(
    placement_type = placement_type,
    qtl_list = qtl_of_interest,
    sectioned_marker_list = sectioned_marker_list,
    chromosome_size = chromosome_size
  )


  for (qtl_i in 1:num_qtl) {
    if (placement_type == 'centered') {
      qtl_ext_length = qtl_of_interest$Length[qtl_i] / 2
    } else if (placement_type == 'extension') {
      qtl_ext_length = qtl_of_interest$Length[qtl_i]
    }
    half_qtl_length <- as.numeric(qtl_of_interest$Length[qtl_i]) / 2
    for (rep_i in 1:num_repetitions) {
      for (gene_i in 1:num_genes) {
        sampled_gene_row <- sample(1:length(whole_genome_gene_dist[, 1]), 1)
        gene_list$Locus[gene_i] <- whole_genome_gene_dist$GeneMiddle[sampled_gene_row]
        gene_list$Chromosome[gene_i] <-  whole_genome_gene_dist$Chromosome[sampled_gene_row]
        sampled_gene_locus <- gene_list$Locus[gene_i]
        sampled_chr_number <- gene_list$Chromosome[gene_i]

        # TODO was this instead
        # markers_on_chromosome <-
        #   as.data.frame(as.matrix(eval(as.name(
        #     sectioned_marker_list[chr_number]
        #   ))))
        markers_on_chromosome <- sectioned_marker_list[chr_number]

        chr_marker_positions <-
          as.integer(as.character(markers_on_chromosome$Base))

        known_likelihood <- qtl_probs[qtl_i]


        range_l <- sampled_gene_locus - half_qtl_length
        range_r <- sampled_gene_locus + half_qtl_length

        markers_l <-
          markers_on_chromosome[which(
            chr_marker_positions + half_qtl_length <= chromosome_size$RightmostMarker[sampled_chr_number] &
              chr_marker_positions <= sampled_gene_locus &
              chr_marker_positions >= range_l
          ), 3]
        length_markers_l <- length(markers_l)

        markers_r <-
          markers_on_chromosome[which(
            chr_marker_positions - half_qtl_length >= chromosome_size$LeftmostMarker[sampled_chr_number] &
              # It's not >= so we don't double-count locus
              chr_marker_positions > as.numeric(sampled_gene_locus) &
              chr_marker_positions <= range_r
          ), 3]
        length_markers_r <- length(markers_r)

        num_markers_avail <- length_markers_l + length_markers_r
        num_markers_expect <- num_markers_avail * known_likelihood
        if (length(num_markers_expect) == 0) {  # TODO is this necessary?
          num_markers_expect <- 0
        }
        gene_list$EGN[gene_i] <-
          gene_list$EGN[gene_i] + num_markers_expect
      }
      sim_egn <- sum(gene_list$EGN)
      output[qtl_i, rep_i] <- sim_egn
      gene_list$EGN <- 0
    }
    setTxtProgressBar(pb, qtl_i / length(num_qtl))
  }





  SDs <- apply(toOutput, 1, sd)
  upper <- c()
  lower <- c()

  for (QTL in 1:length(QTLofInterest$Length)) {
    Zvalue <-
      qnorm(.025 / as.numeric(QTLofInterest$Number_Trait_QTL[QTL]),
            lower.tail = FALSE)
    mu <- rowMeans(toOutput)[QTL]
    upper <- c(upper, mu + Zvalue * SDs[QTL])
    lower <- c(lower, mu - Zvalue * SDs[QTL])
  }


  CIs <-
    as.data.frame(as.matrix(
      cbind(
        QTLofInterest$Length,
        QTLofInterest$N_Genes,
        rowMeans(toOutput),
        SDs,
        lower,
        upper
      )
    ))
  toOutput$QTL_Length <- QTLofInterest$Length
  toOutput$Observed_Value <- QTLofInterest$N_Genes
  toOutput <-
    toOutput[, c((length(toOutput) - 1), length(toOutput), 1:(length(toOutput) -
                                                                2))]
  colnames(toOutput)[3:length(toOutput)] <-
    paste0("Simulation_Round_", 1:(length(toOutput) - 2))
  Simulation_Environment$Simulation_DataFrame <- toOutput


  colnames(CIs) <-
    c("QTL",
      "Observed Value",
      "Mean",
      "SEM",
      "Lower 95% CI",
      "Upper 95% CI")
  return(CIs)


}
