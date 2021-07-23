context("model selection")

## Simulation of an input data matrix

npc <- c(50, 40) # nodes per class
KG <- c(2, 3) # classes
nm <- npc * KG # nodes
Z <- diag( KG[1]) %x% matrix(1, npc[1], 1)
W <- diag(KG[2]) %x% matrix(1, npc[2], 1)
L <- 70 * matrix( runif( KG[1] * KG[2]), KG[1], KG[2])
M_in_expectation <- Z %*% L %*% t(W)
size <- 50
M <- matrix(
  rnbinom(
    n = length(as.vector(M_in_expectation)),
    mu = as.vector(M_in_expectation), size = size)
  , nm[1], nm[2])
rownames(M) <- paste("OTU", seq_len(nrow(M)), sep = "_")
colnames(M) <- paste("S", seq_len(ncol(M)), sep = "_")

res <- cobiclust(M, K = 2, G = 3, nu_j = rep(1,120), a = 1/size, cvg_lim = 1e-5, akg = FALSE)

#crit <- selection_criteria(res, K = 2, G = 3)

##################################################################
#---------------------selection-----------------------------
##################################################################
test_that("Check that selection is running and robust",  {

  expect_error(selection_criteria(res, G = 1))
  expect_error(selection_criteria(res, K = 1))
 # expect_type(selection_criteria(res, K = 2, G = 3),"vector")
  expect_equal(length(selection_criteria(res, K = 2, G = 3)), 7)
})
