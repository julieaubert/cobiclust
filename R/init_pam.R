# cobiclust R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
####################################################################################
#' Initialisation of the co-clusters by partitioning around medoids method.
#'
#' @param x The output of the cobiclust function.
#' @param nu_j a vector of . The length is equal to the number of colums.
#' @param a an numeric.
#' @param K an integer specifying the number of groups in rows.
#' @param G an integer specifying the number of groups in columns.
#' @param akg a logical variable indicating whether to use a common dispersion parameter (akg = FALSE) or a dispersion parameter per cocluster (akg = TRUE).
#' @return A list of
#' \describe{
#' \item{\code{nu_j}}{nu_j.}
#' \item{\code{mu_i}}{mu_i.}
#' \item{\code{t_jg}}{t_jg.}
#' \item{\code{s_ik}}{s_ik.}
#' \item{\code{pi_c}}{pi.}
#' \item{\code{rho_c}}{rho.}
#' \item{\code{a}}{a.}
#' \item{\code{exp_utilde}}{exp_utilde.}
#' \item{\code{exp_logutilde}}{exp_logutilde.}
#' \item{\code{alpha_c}}{alpha.}
#'}
#' @import cluster
#' @importFrom cluster pam
#' @importFrom stats median
#' @export
#' @keywords internal
#'
init_pam <- function(x, nu_j = NULL, a = NULL, K = K, G = G, akg = FALSE){
  n <- nrow(x)
  m <- ncol(x)
  # Initialisation of mu_i ----------------------------------------------------------------
  mu_i <- rowMeans(x)
  # Calculation of nu_j -------------------------------------------------------------------
  if (is.null(nu_j)) {
    nu_j <- colSums(x) / mean(colSums(x))
  } else nu_j <- nu_j

  y <- sqrt(x + 0.5)

  # Initialisation of the co-clusters by partitioning around medoids method ---------------
  init.rows <- pam(y, k = K)
  pi_c <- tabulate(init.rows$clustering) / n

  init.cols <- pam(t(y), k = G)
  rho_c <- tabulate(init.cols$clustering) / m

  if (is.matrix(nu_j)) {
    tmp0 <- sapply(1:m, FUN = function(j) sapply(1:n, FUN = function(i) x[i,j] / (mu_i[i] * nu_j[i,j])))
  } else {
  tmp0 <- t(scale(t(scale(x,center=FALSE,scale=nu_j)),center = FALSE, scale = mu_i))
  }
  alpha_c <- sapply(1:G, FUN = function(y) sapply(1:K, FUN = function(x)
    mean(tmp0[which(init.rows$clustering == x), which(init.cols$clustering == y)])))

  t_jg <- matrix(nrow = m, ncol = G, 0)
  for (i in 1:m){
    t_jg[i, init.cols$clustering[i]] <- 1
  }

  s_ik <- matrix(nrow = n, ncol = K, 0)
  for (i in 1:n){
    s_ik[i, init.rows$clustering[i]] <- 1
  }


  if (is.null(a)){
    # Initialisation of a (gamma parameter) ----------------------------------------------
    m_co <- (t(s_ik) %*% x %*% t_jg) / (t(s_ik) %*% matrix(1, nrow = n, ncol = m) %*% t_jg)
    tmp <-   sapply(1:G, FUN = function(g) sapply(1:K, FUN = function(k)
      sum((x[init.rows$clustering == k, init.cols$clustering == g] - m_co[k, g])^2)))
    var_co <- tmp/((t(s_ik)%*%matrix(1,nrow=n,ncol=m)%*%t_jg)-1)

    a <- (var_co - m_co) / (m_co * m_co)

  }
  if (akg == FALSE){
    a <- max(median(a[is.finite(a)]), 0.1)
    if (is.matrix(nu_j)){
      b_tilde <- a + t(sapply(1:length(mu_i), FUN=function(j)
        (s_ik[j,] * mu_i[j]) %*% tcrossprod(alpha_c, t_jg * nu_j[j,])))
    } else {
    b_tilde <- a + (s_ik * mu_i) %*% tcrossprod(alpha_c, t_jg * nu_j)
    }

    a_tilde <- a + rowSums(s_ik) %*% diag(nrow = 1, ncol = 1, 1) %*% rowSums(t_jg) * x

  } else {
    a <- apply(a, c(1,2), FUN=function(y) max(y[is.finite(y)], 0.1))
    if (is.matrix(nu_j)){
    b_tilde <- tcrossprod(s_ik %*% a, t_jg) + t(sapply(1:length(mu_i), FUN = function(j) (s_ik[j,] * mu_i[j]) %*% tcrossprod(alpha_c, t_jg * nu_j[j, ])))
    } else {
      b_tilde <- tcrossprod(s_ik %*% a, t_jg) + t(sapply(1:length(mu_i), FUN = function(j) (s_ik[j,] * mu_i[j]) %*% tcrossprod(alpha_c, t_jg *nu_j)))
    }
    a_tilde <- tcrossprod(s_ik %*% a, t_jg) + rowSums(s_ik) %*% diag(nrow = 1, ncol = 1, 1) %*% rowSums(t_jg) * x

  }
  # }
  exp_utilde <- a_tilde / b_tilde
  exp_logutilde <- digamma(a_tilde) - log(b_tilde)
  # Output
  strategy <- list(akg = akg, cvg_lim = NULL)
  parameters <- list(alpha = alpha_c, pi = pi_c, rho = rho_c, mu_i = mu_i, nu_j = nu_j, a = a)
  info <- list(s_ik = s_ik, t_jg = t_jg, exp_logutilde = exp_logutilde, exp_utilde = exp_utilde)
  output <- list(data = x, K = K, G = G, classification = NULL,
                 strategy = strategy, parameters = parameters, info = info)
  class(output) <- append(class(output), "cobiclustering")
  invisible(output)

}
