context("Test splitting of underscore-delimited input ensembl field")

set.seed(123)
df <- data.table(Transcript=
                     c("ENSMUST00000027036_Lypla1_mm10_chr1_4844962_4846739_+_utr_4845016_4846739::chr1:4844962-4846739(+)",
                       "ENSMUST00000081551_Tcea1,ENSMUST00000165720_Tcea1_mm10_chr1_4896355_4897910_+_utr_4896364_4897910::chr1:4896355-4897910(+)"),
                 Length=c(100, 1000),
                 SampleA=runif(2),
                 SampleB=runif(2),
                 SampleC=rnorm(2),
                 stringsAsFactors=FALSE
)

test_that("'Transcript' field is split into components for 1-row data frame", {
    exp.df <- df[1,]
    separate_ensembl_field(exp.df)
    # expect_equal(ncol(exp.df), 12)
    # expect_match(exp.df$Transcript[1], "^ENS.*\\d$")
    expected_cols <- union(c("Transcript", "Gene", "Chr", "LastExon.Start",
                             "LastExon.End", "Strand", "UTR3.Start", "UTR3.End"),
                           colnames(df)) 
    expect_true(all(expected_cols %in% colnames(exp.df)))
})

test_that("'Transcript' field is split into components for multi-line data frame", {
    exp.df <- copy(df)
    separate_ensembl_field(exp.df)
    # expect_equal(ncol(exp.df), 10)
    expect_is(exp.df$Transcript, "character")
    expected_cols <- union(c("Transcript", "Gene", "Chr", "LastExon.Start",
                             "LastExon.End", "Strand", "UTR3.Start", "UTR3.End"),
                           colnames(df)) 
    expect_true(all(expected_cols %in% colnames(exp.df)))
})

context("Test non-Ensembl IDs")
set.seed(123)
df <- data.table(Transcript=
                     c("XY.00000027036_000..2_mm10_chr1_4844962_4846739_+_utr_4845016_4846739::chr1:4844962-4846739(+)",
                       "FF.22_1.z,0101_1.z_mm10_chr1_4896355_4897910_+_utr_4896364_4897910::chr1:4896355-4897910(+)"),
                 Length=c(100, 1000),
                 SampleA=runif(2),
                 SampleB=runif(2),
                 SampleC=rnorm(2),
                 stringsAsFactors=FALSE
)

test_that("Non-Ensembl transcript ID is accepted", {
    exp.df <- copy(df)
    separate_ensembl_field(exp.df)
    expect_equal(exp.df$Gene, c("000..2", "1.z"), "Genes do not match")
    expect_equal(exp.df$Transcript, c("XY.00000027036", "FF.22,0101"),
                 "Transcripts do not match")
})

