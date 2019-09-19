context("Test APA_ID generation")

x <- c("ENSG00000001", "ENSG00000001", "ENSG00000999", "ENSG00000001")

test_that("A number is added to Ensembl ID", {
    expect_equal(apa_id(x), paste(x, c(1,2,1,1), sep="_"))
    expect_equal(apa_id(sort(x)), paste(sort(x), c(1:3,1), sep="_"))
    expect_equal(apa_id(x[3]), paste(x[3], "1", sep="_"))
})

context("Test APA_ID suffix update")

test_that("A suffix is added for single UTR", {
    expect_equal(update_apa_id(x[3], 5, 10), paste(x[3], "S", sep="_"))
    expect_equal(update_apa_id(x[1:2], c(5, 5), c(10, 9)), paste(x[1:2], c("D", "P"), sep="_"))
    expect_equal(update_apa_id(x[1:2], c(5, 5), c(4, 3)), paste(x[1:2], c("P", "D"), sep="_"))
})
