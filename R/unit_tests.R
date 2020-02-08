source('./SPQV_refactored.r', chdir=TRUE)
library(testthat)

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

# SectionMarkers ####
test_that("SectionMarkers creates correct number of lists", {
  output <- SectionMarkers(MARKER_LIST, nrow(CHROMOSOME_SIZE))
  expect_length(output, 2)
})

# FilterGenes ####
test_that("FilterGenes filters by trait", {
  output <- FilterGeneList("other_trait", GENE_LIST, MARKER_LIST)
  expect_setequal(
    output$ID,
    c("1.E")
  )
})

test_that("FilterGenes doesn't remove tandem arrays if drop_tandem is false", {
  output <- FilterGeneList(TRAIT, GENE_LIST, MARKER_LIST, drop_tandem=FALSE)
  expect_setequal(
    output$ID,
    c("1.a", "1.b", "1.c", "1.d", "2.a", "2.b", "2.c")
  )
})

test_that("FilterGenes removes tandem arrays by default", {
  output <- FilterGeneList(TRAIT, GENE_LIST, MARKER_LIST)
  expect_true(
    all(output$ID == c("1.a", "1.b", "1.d", "2.a", "2.b", "2.c")) |
    all(output$ID == c("1.a", "1.c", "1.d", "2.a", "2.b", "2.c"))
  )
})


# QTLPlacementProbabilities ####
test_that("QTLPlacementProbabilities gets correct probabilites for 'extension'", {
  output <- QTLPlacementProbabilities(
    QTL_LIST, "extension",
    SectionMarkers(MARKER_LIST, nrow(CHROMOSOME_SIZE)), CHROMOSOME_SIZE)
  expect_equal(
    output,
    c(1/14, 1/14, 1/11, 1/3)
  )
})

test_that("QTLPlacementProbabilities gets correct probabilites for 'centered', and counts each marker only once", {
  output <- QTLPlacementProbabilities(
    QTL_LIST, "centered",
    SectionMarkers(MARKER_LIST, nrow(CHROMOSOME_SIZE)), CHROMOSOME_SIZE)
  expect_equal(
    output,
    # 5,   10,  30,  80
    c(1/5, 1/5, 1/4, 1/2)
  )
})

# CountGenesFound ####
test_that("CountGenesFound doesn't allow unfiltered gene lists", {
  expect_error(
    CountGenesFound(QTL_LIST, TRAIT, GENE_LIST, MARKER_LIST)
  )
})

test_that("CountGenesFound finds correct genes", {
  output <- CountGenesFound(
    QTL_LIST, TRAIT, 
    FilterGeneList(TRAIT, GENE_LIST, MARKER_LIST, drop_tandem = FALSE), 
    MARKER_LIST
    )
  expect_equal(
    output$NumGenes, c(0, 1, 2, 4)
  )
  expect_equal(
    output$FoundGeneIDs, c("0", "1.d", "1.b and 1.c", "1.a and 1.b and 1.c and 1.d") 
  )
})


# GeneFoundLikelihood ####
test_that("GeneFoundLikelihood gets correct likelihood for 'extension'", {
  output <- GeneFoundLikelihood(
    gene_chr = 1, 
    gene_locus = 40, 
    qtl_ext_length = 30, 
    placement_type = "extension",
    per_marker_likelihood = 0.1, 
    sectioned_marker_list = SectionMarkers(MARKER_LIST, nrow(CHROMOSOME_SIZE)),
    chromosome_size = CHROMOSOME_SIZE
  )
  expect_equal(
    output, 0.1 * 3
  )
})

test_that("GeneFoundLikelihood gets correct likelihood for 'centered'", {
  output <- GeneFoundLikelihood(
    gene_chr = 1, 
    gene_locus = 70, 
    qtl_ext_length = 30, 
    placement_type = "centered",
    per_marker_likelihood = 0.1, 
    sectioned_marker_list = SectionMarkers(MARKER_LIST, nrow(CHROMOSOME_SIZE)),
    chromosome_size = CHROMOSOME_SIZE
  )
  expect_equal(
    output, 0.1 * 2
  )
})

test_that("GeneFoundLikelihood allows overhang in other direction for 'extension'", {
  output <- GeneFoundLikelihood(
    gene_chr = 1, 
    gene_locus = 20, 
    qtl_ext_length = 30, 
    placement_type = "extension",
    per_marker_likelihood = 0.1, 
    sectioned_marker_list = SectionMarkers(MARKER_LIST, nrow(CHROMOSOME_SIZE)),
    chromosome_size = CHROMOSOME_SIZE
  )
  expect_equal(
    output, 0.1 * 1
  )
})

test_that("GeneFoundLikelihood doesn't allow overhang in other direction for 'centered'", {
  output <- GeneFoundLikelihood(
    gene_chr = 1, 
    gene_locus = 20, 
    qtl_ext_length = 30, 
    placement_type = "centered",
    per_marker_likelihood = 0.1, 
    sectioned_marker_list = SectionMarkers(MARKER_LIST, nrow(CHROMOSOME_SIZE)),
    chromosome_size = CHROMOSOME_SIZE
  )
  expect_equal(
    output, 0.1 * 0
  )
})

test_that("GeneFoundLikelihood returns 0 when gene not reachable", {
  output <- GeneFoundLikelihood(
    gene_chr = 1, 
    gene_locus = 50, 
    qtl_ext_length = 2, 
    placement_type = "centered",
    per_marker_likelihood = 0.1, 
    sectioned_marker_list = SectionMarkers(MARKER_LIST, nrow(CHROMOSOME_SIZE)),
    chromosome_size = CHROMOSOME_SIZE
  )
  expect_equal(
    output, 0
  )
})

test_that("GeneFoundLikelihood returns 0 when QTL can't be placed", {
  output <- GeneFoundLikelihood(
    gene_chr = 2, 
    gene_locus = 30, 
    qtl_ext_length = 100, 
    placement_type = "centered",
    per_marker_likelihood = 0.1, 
    sectioned_marker_list = SectionMarkers(MARKER_LIST, nrow(CHROMOSOME_SIZE)),
    chromosome_size = CHROMOSOME_SIZE
  )
  expect_equal(
    output, 0
  )
})

# SPQValidate ####



