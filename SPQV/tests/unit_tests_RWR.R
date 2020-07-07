library(testthat)
devtools::document("..")

# Two chromosomes, with this structure (markers are numbers, genes are letters, each is 10 bp apart):
# Chr 1: a  1  2  b  c  3  d  4  E  5  (has: gene w/o marker on end (a); tandem array (bc); different trait (E))
#        10 20 30 40 50 60 70 80 90 120
# Chr 2: 1  a  2  b  3  c  4 (regular structure)
#        10 20 30 40 50 60 70

# Test data init ####

TRAIT <- "test_trait"

MARKER_LIST <- data.frame(
  ID =         c("1.1", "1.2", "1.3", "1.4", "1.5", "2.1", "2.2", "2.3", "2.4"),
  Chromosome = c( 1,     1,     1,     1,     1,     2,     2,     2,     2     ),
  Base =       c( 20,    30,    60,    80,    120,   10,    30,    50,    70    ),
  stringsAsFactors = FALSE
)

GENE_LIST <- data.frame(
  ID =         c("1.a", "1.b", "1.c", "1.d", "1.E", "2.a", "2.b", "2.c"),
  Chromosome = c( 1,     1,     1,     1,     1,     2,     2,     2     ),
  Base =       c( 10,    40,    50,    70,    90,    20,    40,    60    ),
  stringsAsFactors = FALSE
)
GENE_LIST$Trait <- TRAIT
GENE_LIST[which(GENE_LIST$ID == "1.E"), "Trait"] = "other_trait"

CHROMOSOME_SIZE <- data.frame(
  Chromosome =      c( 1,     2),
  LeftmostMarker =  c( 20,    10),
  RightmostMarker = c( 120,   70),
  Length =          c( 110,   80)
)

QTL_LIST <- data.frame(
  Length =          c( 5,    10,    30,    80), # don't actually care if on markers or if overlap
  Chromosome =      c( 1,     1,     1,     1),
  LeftmostMarker =  c( 20,    60,    20,    10),
  RightmostMarker = c( 25,    70,    50,    90)
)
QTL_LIST$Trait <- TRAIT
QTL_LIST$Treatment <- "UNUSED_FOR_NOW"
QTL_LIST$Method <- "UNUSED_FOR_NOW"
QTL_LIST$ExptType <- "UNUSED_FOR_NOW"

for (col in c("Chromosome", "Base")) {
  MARKER_LIST[[col]] <- as.integer(MARKER_LIST[[col]])
  GENE_LIST[[col]] <- as.integer(GENE_LIST[[col]])
}
for (col in c("Chromosome", "Length", "LeftmostMarker", "RightmostMarker")) {
  CHROMOSOME_SIZE[[col]] <- as.integer(CHROMOSOME_SIZE[[col]])
  QTL_LIST[[col]] <- as.integer(QTL_LIST[[col]])
}

WGD <- GENE_LIST
WGD$GeneMiddle <- WGD$Base

# (
#   markers_only,  # vs. random placement
#   skip_hangovers,  # vs bounceback
#   bidirectional,  # vs left-to-right only
#   drop_tandem,
#   qtl_list,
#   gene_list,
#   marker_list,
#   chromosome_size
# )
QTL_LIST <- CountGenesFound(
  QTL_LIST,
  TRAIT,
  FilterGeneList(trait = TRAIT, gene_list = GENE_LIST, marker_list = MARKER_LIST),
  MARKER_LIST
)


test_that("RWR runs", {
    output <- RWR(
      FALSE,
      FALSE,
      FALSE,
      FALSE,
      QTL_LIST,
      GENE_LIST,
      MARKER_LIST,
      CHROMOSOME_SIZE,
      10
    )
    expect_length(output, nrow(QTL_LIST))
  }
)

