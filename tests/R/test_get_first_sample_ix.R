context("Test getting the index of the first sample in merged data")

set.seed(123)
nr <- 10
N <- 3
df <- data.frame(APA_ID = letters[1:nr],
                 Ensembl_Gene = sample(letters, nr),
                 Gene_Name = rainbow(nr),
                 Chr = "chr1",
                 Start = round(runif(nr, 1000, 2000)),
                 End = round(runif(nr, 3000, 5000)),
                 Strand = "+",
                 Length = round(runif(nr, 100, 500)),
                 SampleA = rnorm(nr),
                 SampleB = rnorm(nr),
                 SampleC = rnorm(nr)
)

test_that("Gets index of first sample column in data frame", {
    expected <- 9
    expect_equal(get_first_sample_ix(df), expected)
})


context("Test getting the number of samples from merged data")

test_that("Can calculate number of samples in data frame", {
    expected <- N
    expect_equal(get_num_samples(df), expected)
})
