library(stringr)
context("Test formatting of single-Ensembl IDs")

id <- "ENSMUST0000011043_ENSMUSG00000111044_mm9_chr1"

test_that("A single Ensembl ID remains unchanged", {
    expect_equal(format_multi_ensembl_ids(id), id)
})

test_that("A vector of single Ensembl IDs remains unchanged", {
    expect_equal(format_multi_ensembl_ids(c(id, id)), c(id, id))
})

context("Test formatting of multi-Ensembl IDs")

id <- "ENSMUST00000111043_ENSMUSG00000048482,ENSMUST00000111044_ENSMUSG00000048482_mm9_chr1"
expected <- "ENSMUST00000111043,ENSMUST00000111044_ENSMUSG00000048482_mm9_chr1"

test_that("A suffix is added for single UTR", {
    expect_equal(format_multi_ensembl_ids(id), expected)
})

test_that("A vector of multi-Ensembl IDs is re-formatted", {
    expect_equal(format_multi_ensembl_ids(c(id, id)), c(expected, expected))
})

context("Test vector of mixed (single and multi) Ensembl IDs")

id <- c("ENSMUST00000111043_ENSMUSG00000048482,ENSMUST00000111044_ENSMUSG00000048482_hg19_chr1",
        "ENSMUST00000100011_ENSMUSG00000048481_hg19_chr2",
        "ENSMUST00000111043_ENSMUSG00000048480,ENSMUST00000111044_ENSMUSG00000048482_hg19_chr1")

test_that("A vector of mixed Ensembl IDs", {
    expect_equal(format_multi_ensembl_ids(id), 
                 c("ENSMUST00000111043,ENSMUST00000111044_ENSMUSG00000048482_hg19_chr1",
                   id[2],
                   "ENSMUST00000111043,ENSMUST00000111044_ENSMUSG00000048480,ENSMUSG00000048482_hg19_chr1")
    )
})


context("Test non-Ensembl IDs")

test_that("Non-Ensembl transcript ID is accepted", {
    id <- c("XY.00000027036_000..2_mm10_chr1",
            "FF.22_1.z,0101_1.z_mm10_chr1",
            "1_.,2_.,3_._mm10_chr1")
    expect_equal(format_multi_ensembl_ids(id),
                 c("XY.00000027036_000..2_mm10_chr1",
                   "FF.22,0101_1.z_mm10_chr1",
                   "1,2,3_._mm10_chr1"))
})

test_that("Underscore at beginning of ID will fail", {
    expect_error(format_multi_ensembl_ids("_ENSMUST00000111043_ENSMUSG00000048482,ENSMUST00000111044_ENSMUSG00000048482_hg19_chr1"))
})

context("Test chromosomes without chr prefix")
test_that("Chromosome without chr prefix is accepted", {
    id <- c("ENSMUST00000111043_ENSMUSG00000048482,ENSMUST00000111044_ENSMUSG00000048482_hg19_1",
            "ENSMUST00000100011_ENSMUSG00000048481_hg19_2",
            "ENSMUST00000111043_ENSMUSG00000048480,ENSMUST00000111044_ENSMUSG00000048482_hg19_z")
    expect_equal(format_multi_ensembl_ids(id),
                 c("ENSMUST00000111043,ENSMUST00000111044_ENSMUSG00000048482_hg19_1",
                   id[2],
                   "ENSMUST00000111043,ENSMUST00000111044_ENSMUSG00000048480,ENSMUSG00000048482_hg19_z")
    )
})

context("Test unk species")
test_that("Non-hg19 and non-mm10 species are allowed as unk", {
    expect_equal(format_multi_ensembl_ids("1_0,2_0_unk_chr1"),
                 "1,2_0_unk_chr1")
})

test_that("Complex unknown species and non-standard chr is accepted", {
    expect_equal(format_multi_ensembl_ids("1_0,2_0_ut_z"),
                 "1,2_0_ut_z")
})

