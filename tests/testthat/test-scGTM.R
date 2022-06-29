library(scGTM)


test_that("scGTM works", {
  data("df")
  res1 <- scGTM::scGTM(gene_index = 1, t=df$Time, y1=df$Gene1)
  expect_equal(length(res1), 22)

  res2 <- scGTM::scGTM(gene_index = 1, t=df$Time, y1=df$Gene1, marginal="Poisson")
  expect_equal(length(res2), 22)

  res3 <- scGTM::runscGTM(gene.vec = 1:3, t=df$Time, y1=df[,3:5],gene_name=colnames(df[,3:5]))
  expect_equal(dim(res3)[2], 3)
})
