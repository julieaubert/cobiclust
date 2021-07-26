context("testing VEM steps")

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

##################################################################
#---------------------INITIALISATION -----------------------------
##################################################################
K <- KG[1]
G <- KG[2]
x <- M
a <- 1/size
nu_j <- rep(1,120)
res_init <- init_pam(x = x, nu_j = nu_j, a = a, K = K, G = G, akg = FALSE)
#res_init_akg <- init_pam(x = M, nu_j = nu_j, a = 1/size, K = KG[1], G = KG[2], akg = TRUE)
# a doit etre coherent avec akg si TRUE a doit etre de la dimension K*G


test_that("The initialisation is giving  elements of right size", {

  expect_type(res_init,"list")
  expect_true(is.null(res_init$classification))

  ########### Dimension and type of the objects

  expect_type(res_init$info,"list")
  expect_type(res_init$parameters,"list")
  expect_type(res_init$strategy,"list")

  expect_equal(dim(res_init$info$s_ik), c(nrow(x), K))
  expect_equal(dim(res_init$info$t_jg), c(ncol(x), G))
  expect_equal(dim(res_init$info$exp_logutilde), dim(x))
  expect_equal(dim(res_init$info$exp_utilde), dim(x))

  expect_equal(dim(res_init$parameters$alpha), c(K, G))
  expect_equal(length(res_init$parameters$pi), K)
  expect_equal(length(res_init$parameters$rho), G)
  expect_equal(length(res_init$parameters$mu_i), nrow(x))

  expect_true(res_init$parameters$a > 0)
})


##################################################################
#---------------------cobiclust ARGUMENTS -----------------------------
##################################################################
test_that("Invalid arguments throw error", {

  expect_error(cobiclust(M, K = 2, G = 3, nu_j = rep(1,120), a = 1/size, cvg_lim = 1e-5, akg = TRUE))
  expect_error(cobiclust(M, K = 2, G = 3, nu_j = rep(1,120), a = c(1,2), cvg_lim = 1e-5, akg = FALSE))
  expect_error(cobiclust(M, K = 0, G = 3, nu_j = rep(1,120), a = 1/size, cvg_lim = 1e-5, akg = FALSE))
  expect_error(cobiclust(M, K = 2, G = "3 groupes", nu_j = rep(1,120), a = 1/size, cvg_lim = 1e-5, akg = FALSE))
  expect_error(cobiclust(M, K = 2, G = 3, nu_j = rep(1,100), a = 1/size, cvg_lim = 1e-5, akg = FALSE))
  expect_error(cobiclust(M, K = 2, G = 3, nu_j = matrix(ncol = 2, rep(1,240)), a = 1/size, cvg_lim = 1e-5, akg = FALSE))


})

##################################################################
#---------------------cobiclust RESULTS -----------------------------
##################################################################
test_that("Check that cobiclust is running and robust",  {

  res <- cobiclust(M, K = 2, G = 3, nu_j = rep(1,120), a = NULL, cvg_lim = 1e-5, akg = TRUE)

 # expect_is(res, "cobiclustering")
  expect_type(res, "list")
  expect_s3_class(res, "cobiclustering")
  expect_equal(length(res$parameters$pi), K)
  expect_equal(length(res$parameters$rho), G)
})
