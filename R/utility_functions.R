# cobiclust R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
####################################################################################
#' Calculate E_q(U_{ijr}) and E_q(log U_{ijr}) when U_ijr follows a Gamma(a,a)
#'
#' @param s_ik s_ik.
#' @param t_jg t_jg.
#' @param x a matrix of observations. Columns correspond to biological samples and rows to microorganisms.
#' @param a a.
#' @param mu_i mu_i.
#' @param nu_j nu_j.
#' @param alpha_c alpha_c.
#' @return A list of 4 elements.
#' \describe{
#' \item{\code{a_tilde}}{a_tilde.}
#' \item{\code{b_tilde}}{b_tilde.}
#' \item{\code{exp_utilde}}{exp_utilde.}
#' \item{\code{exp_logutilde}}{exp_logutilde.}
#'}
#' @export
#' @keywords internal

qu_calculation <- function(s_ik = s_ik, t_jg = t_jg, x = x, mu_i = mu_i, nu_j = nu_j, alpha_c = alpha_c, a = a){
  a_tilde <- a + rowSums(s_ik) %*% diag(nrow = 1, ncol = 1, 1)%*%rowSums(t_jg)*x
  if (is.matrix(nu_j)){
     b_tilde <- a + t(sapply(1:length(mu_i),FUN=function(j) (s_ik[j, ] * mu_i[j]) %*% tcrossprod(alpha_c, t_jg * nu_j[j,])))

  } else {
    b_tilde <- a + (s_ik * mu_i) %*% tcrossprod(alpha_c, t_jg * nu_j)
  }
  exp_utilde <- a_tilde / b_tilde
  exp_logutilde <- digamma(a_tilde) - log(b_tilde)
  return(list(a_tilde = a_tilde, b_tilde = b_tilde, exp_utilde = exp_utilde,
              exp_logutilde = exp_logutilde))
}

#' Calculate E_q(U_{ijr}) and E_q(log U_{ijr}) when U_ijr follows a Gamma(a_{k,g},a_{k,g})
#'
#' @param s_ik s_ik.
#' @param t_jg t_jg.
#' @param x a matrix of observations. Columns correspond to biological samples and rows to microorganisms.
#' @param a a0.
#' @param mu_i mu_i.
#' @param nu_j nu_j.
#' @param alpha_c alpha_c.
#' @return A list of 4 elements.
#' \describe{
#' \item{\code{a_tilde}}{a_tilde.}
#' \item{\code{b_tilde}}{b_tilde.}
#' \item{\code{exp_utilde}}{exp_utilde.}
#' \item{\code{exp_logutilde}}{exp_logutilde.}
#'}
#' @export
#' @keywords internal

qukg_calculation <- function(s_ik = s_ik, t_jg = t_jg, x = x, mu_i = mu_i, nu_j = nu_j, alpha_c = alpha_c, a = a){
  a_tilde <- tcrossprod(s_ik %*% a, t_jg) + rowSums(s_ik) %*%
    diag(nrow = 1, ncol = 1, 1) %*% rowSums(t_jg) * x
  if (is.matrix(nu_j)){
    #b_tilde <- a + t(sapply(1:length(mu_i), FUN=function(j)
    #  (s_ik[j,] * mu_i[j]) %*% tcrossprod(alpha_c, t_jg * nu_j[j,])))
    b_tilde <- tcrossprod(s_ik %*% a, t_jg) + t(sapply(1:length(mu_i), FUN=function(j)(s_ik[j,]*mu_i[j]) %*%
      tcrossprod(alpha_c, t_jg * nu_j[j,])))
  } else {
    b_tilde <- tcrossprod(s_ik %*% a, t_jg) + (s_ik*mu_i) %*%
      tcrossprod(alpha_c, t_jg * nu_j)
  }

  exp_utilde <- a_tilde / b_tilde
  exp_logutilde <- digamma(a_tilde) - log(b_tilde)
  return(list(a_tilde = a_tilde, b_tilde = b_tilde, exp_utilde = exp_utilde,
              exp_logutilde = exp_logutilde))
}

#' Calculate the matrix (alpha_{k,g})
#'
#' @param s_ik s_ik.
#' @param t_jg t_jg.
#' @param nu_j nu_j.
#' @param mu_i mu_i.
#' @param K K.
#' @param G G.
#' @param x a matrix of observations. Columns correspond to biological samples and rows to microorganisms.
#' @param exp_utilde exp_utilde.
#' @return a matrix of dimension (\code{K},\code{G}) of the terms of interactions.
#' @export
#' @keywords internal
#'

alpha_calculation <- function(s_ik = s_ik, t_jg = t_jg, nu_j = nu_j, mu_i = mu_i,
                              K = K, G = G, x = x, exp_utilde = exp_utilde){
if (is.matrix(nu_j)){
  denum <- crossprod((s_ik * mu_i), ((nu_j * exp_utilde) %*% t_jg))
  #denum <- crossprod((s_ik * mu_i), t(sapply(1:length(mu_i), FUN=function(j) exp_utilde[j,] %*% (t_jg * mat_nuj[j,]) )))
} else {
  denum <- crossprod((s_ik * mu_i), (exp_utilde %*% (t_jg * nu_j)))
}
alpha_tmp <- crossprod(s_ik, x %*% t_jg) / denum
alpha_c <- K * G * alpha_tmp / sum(alpha_tmp)
return(alpha_c)
}

#' Calculate the lower bound
#'
#' @param x a matrix of observations. Columns correspond to biological samples and rows to microorganisms.
#' @param qu.param qu.param.
#' @param s_ik s_ik.
#' @param pi_c pi_c.
#' @param t_jg t_jg.
#' @param rho_c rho_c.
#' @param mu_i mu_i.
#' @param nu_j nu_j.
#' @param alpha_c a matrix the terms of interactions.
#' @param a a.
#' @param akg a logical variable indicating whether to use a common dispersion parameter (akg = FALSE) or a dispersion parameter per cocluster (akg = TRUE).
#' @return a list of 2 elements.
#' \describe{
#' \item{\code{lb}}{value of the lower bound.}
#' \item{\code{ent}}{value of the entropy term.}
#'}
#'@importFrom stats dnbinom
#' @export
#' @keywords internal
#'

lb_calculation <- function(x = x, qu_param = qu_param, s_ik = s_ik, pi_c = pi_c,
                           t_jg = t_jg, rho_c = rho_c, mu_i = mu_i, nu_j = nu_j,
                           alpha_c = alpha_c, a = a, akg = TRUE) {
  term1 <- sum(s_ik %*% log(pi_c))
  term2 <- sum(t_jg %*% log(rho_c))

  alpha_matrix <- s_ik%*%tcrossprod(alpha_c,t_jg)
  if (akg == FALSE) {
    if (is.matrix(nu_j)) {
      term3  <- sum(sapply(1:ncol(x), FUN = function(j) dnbinom(x[,j], size = a,
                                                                mu = mu_i*nu_j[,j]*alpha_matrix[,j], log = TRUE))) # verif 1/a ou a.

    } else {
      term3  <- sum(sapply(1:ncol(x), FUN = function(j) dnbinom(x[,j], size = a,
                                                                mu = mu_i*nu_j[j]*alpha_matrix[,j], log = TRUE))) # verif 1/a ou a.

    }

  }
  else {
    aa <- tcrossprod(s_ik %*% a, t_jg)
    alpha_matrix <- s_ik%*%tcrossprod(alpha_c,t_jg)
    if (is.matrix(nu_j)){
      term3  <- sum(sapply(1:ncol(x), FUN = function(j) dnbinom(x[,j], size = aa[,j],
                                                                mu = mu_i*nu_j[,j]*alpha_matrix[,j], log = TRUE)))

    } else {
      term3  <- sum(sapply(1:ncol(x), FUN = function(j) dnbinom(x[,j], size = aa[,j],
                                                                mu = mu_i*nu_j[j]*alpha_matrix[,j], log = TRUE)))

    }


  }
  ent_ZW <- - sum(s_ik * log(s_ik)) - sum(t_jg * log(t_jg))

  lb <- term1 + term2 + term3 + ent_ZW

  return(list(lb = lb, ent_ZW = ent_ZW))

}


#' Useful function to estimate the parameter a
#'
#' @param x x.
#' @param nb nb.
#' @param left_bound left_bound.
#' @param right_bound right_bound.
#' @return a numeric.
#' @export
#' @keywords internal
foo_a <- function(x, nb, left_bound, right_bound){
  nb * (1 + log(x) - digamma(x)) + left_bound - right_bound
}

#' Internal function - estimation of the parameter a
#'
#' @param x a matrix of observations. Columns correspond to biological samples and rows to microorganisms.
#' @param y y.
#' @param threshold threshold.
#' @param nb nb.
#' @param left_bound left_bound.
#' @param right_bound right_bound.
#' @return a numeric.
#' @export
#' @keywords internal
dicho <- function(x, y, threshold, nb, left_bound, right_bound){
  while (abs(y - x) >= threshold){
    z <- (x + y) / 2
    if (foo_a(z, nb = nb, left_bound = left_bound, right_bound = right_bound) == 0)
      break
    else if (foo_a(x, nb = nb, left_bound = left_bound, right_bound = right_bound)
             * foo_a(z, nb = nb, left_bound = left_bound, right_bound = right_bound) <= 0) {
      y <- z
    } else {
        x <- z
        }
  }
  return(z)
}

#' Calculate the BIC penalty
#'
#' @param x an object of class biclustering.
#' @return the value of the BIC penalty.
#' @export
#' @keywords internal
penalty <-  function(x){
  K <- x$K
  G <- x$G
  n <- length(x$parameters$mu_i)
  m <- length(x$parameters$nu_j)

  if (x$strategy$akg == TRUE){
    penalty <- 0.5 * (K - 1) * log(n) + 0.5 * (G - 1) * log(m) + 0.5 * (2 * K * G - 1) * log(m * n)
  } else {
    penalty <- 0.5 * (K - 1) * log(n) + 0.5 * (G - 1) * log(m) + 0.5 * (K * G) * log(m * n)
  }
  return(penalty)
}
