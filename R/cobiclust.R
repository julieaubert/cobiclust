# cobiclust R package Copyright INRAE 2020 Universite Paris-Saclay,
# AgroParisTech, INRAE, UMR MIA-Paris, 75005, Paris, France
#' Perform a biclustering adapted to overdispersed count data.
#'
#' @param x the input matrix of observed data.
#' @param K an integer specifying the number of groups in rows.
#' @param G an integer specifying the number of groups in columns.
#' @param nu_j a vector of . The length is equal to the number of colums.
#' @param a an numeric.
#' @param akg a logical variable indicating whether to use a common dispersion parameter (akg = FALSE) or a dispersion parameter per cocluster (akg = TRUE).
#' @param cvg_lim a number specifying the threshold used for convergence criterion (cvg_lim = 1e-05 by default).
#' @param nbiter the maximal number of iterations for the global loop of variational EM algorithm (nbiter = 5000  by default).
#' @param tol the level of relative iteration convergence tolerance (tol = 1e-04  by default).
#' @return An object of class cobiclustering
#' @examples
#' npc <- c(50, 40) # nodes per class
#' KG <- c(2, 3) # classes
#' nm <- npc * KG # nodes
#' Z <- diag( KG[1]) %x% matrix(1, npc[1], 1)
#' W <- diag(KG[2]) %x% matrix(1, npc[2], 1)
#' L <- 70 * matrix( runif( KG[1] * KG[2]), KG[1], KG[2])
#' M_in_expectation <- Z %*% L %*% t(W)
#' size <- 50
#' M<-matrix(
#'  rnbinom(
#'    n = length(as.vector(M_in_expectation)),
#'    mu = as.vector(M_in_expectation), size = size)
#'  , nm[1], nm[2])
#' rownames(M) <- paste('OTU', 1:nrow(M), sep = '_')
#' colnames(M) <- paste('S', 1:ncol(M), sep = '_')
#' res <- cobiclust(M, K = 2, G = 3, nu_j = rep(1,120), a = 1/size, cvg_lim = 1e-5)
#' @seealso \code{\link{cobiclustering}} for the cobiclustering class.
#' @export
#'

cobiclust <- function(x, K = 2, G = 3, nu_j = NULL, a = NULL, akg = FALSE, cvg_lim = 1e-05,
    nbiter = 5000, tol = 1e-04) {

    #---------------------  Some preliminary tests ---------------------

    #--- on the input data
    assertthat::assert_that(is.matrix(x))
    # if (!is.matrix(x)){ stop('x should be a matrix.') }

    #--- on the LBM parameters
    assertthat::assert_that((K < nrow(x)) | (K > 1), msg = "Inadequate number of groups in row. K should be between 1 and the total number of rows of the input data.")

    # if ((K > nrow(x)) | (K < 1)) { stop('Inadequate number of groups in row. K
    # should be between 1 and the total number of rows of the input data.') }

    assertthat::assert_that((G < ncol(x)) | (G > 1), msg = "Inadequate number of groups in columns. G should be between 1 and the total number of rows of the input data.")
    # if ((G > ncol(x)) | (G < 1)) { stop('Inadequate number of groups in columns. G
    # should be between 1 and the total number of rows of the input data.') }

    # --- on the emission parameters
    if (!is.null(nu_j)) {
        if (is.vector(nu_j)) {
            assertthat::assert_that(length(nu_j) == ncol(x), msg = "The dimensions of nu_j and x should be coherent.")
            # if (length(nu_j) != ncol(x)){ stop('The dimensions of nu_j and x should be
            # coherent.') }
        }
        if (is.matrix(nu_j)) {
            assertthat::assert_that(ncol(nu_j) == ncol(x), msg = "The dimensions of nu_j and x should be coherent.")
            assertthat::assert_that(nrow(nu_j) == nrow(x), msg = "The dimensions of nu_j and x should be coherent.")
            # if (ncol(nu_j) != ncol(x)){ stop('The dimensions of nu_j and x should be
            # coherent.') } if (nrow(nu_j) != nrow(x)){ stop('The dimensions of nu_j and x
            # should be coherent.') }
        }
    }
    #---------------------  TESTS on other parameters
    # --- nbiter ---
    assertthat::assert_that(is.numeric(nbiter))
    assertthat::assert_that(nbiter > 2, msg = "Please choose a value greater than 1 for nbiter.")
    # if (nbiter < 2) { stop('Please choose an higher value for nbiter.') } if
    # (!is.integer(nbiter)) { nbiter <- as.integer(nbiter) warning('nbiter was
    # automatically converted to an integer.') }

    # --- tol ---
    assertthat::assert_that(tol > 0 | (tol < 1), msg = "Please choose a value between 0 and 1 for tol.")
    # if (tol < 0 | (tol > 1)) { stop('Please choose a reasonable value for tol.') }

    # --- akg ---
    assertthat::assert_that(is.logical(akg))
    # if (!is.logical(akg)) { stop('akg should be logical.') } --- a ---
    if (!is.null(a)) {
        assertthat::assert_that(is.logical(akg))
        if (isFALSE(akg)) {
            testthat::expect_gt(a, 0)
        } else {
              assertthat::assert_that(length(a) == length(x), msg = "Please chose a value for a coherent with akg. a must be a scalar if aFALSE and a matrix if akg is TRUE")
        }

        # if (a < 0) { stop('a should be positive.') }
    }

    # Parameter initialisation ---------------------------
    res_init <- init_pam(x = x, nu_j = nu_j, a = a, K = K, G = G, akg = akg)

    n <- nrow(x)
    m <- ncol(x)
    nu_j <- res_init$parameters$nu_j
    mu_i <- res_init$parameters$mu_i
    t_jg <- res_init$info$t_jg
    s_ik <- res_init$info$s_ik
    pi_c <- res_init$parameters$pi
    rho_c <- res_init$parameters$rho
    alpha_c <- matrix(nrow = K, ncol = G, res_init$parameters$alpha)
    a0 <- res_init$parameters$a
    exp_utilde <- res_init$info$exp_utilde
    exp_logutilde <- res_init$info$exp_logutilde
    lb <- NULL
    lbtt <- NULL

    # Global EM -----------------------------------------------------------------
    j <- 0
    crit <- 1

    while ((crit > cvg_lim) & (j < nbiter)) {
        j <- j + 1
        # cat('Iteration EM global ',j,'\\n')
        t_old <- t_jg
        s_old <- s_ik
        pi_old <- pi_c
        rho_old <- rho_c
        alpha_old <- alpha_c
        a_old <- a0
        exp_utilde_old <- exp_utilde
        exp_logutilde_old <- exp_logutilde


        # EM on the rows
        # -----------------------------------------------------------------
        i <- 0
        crit_em1 <- 1
        while (crit_em1 > cvg_lim) {
            i <- i + 1
            pi_old <- pi_c
            alpha_old1 <- alpha_c
            # s_ik ------------------
            if (is.matrix(nu_j)) {
                s_ik_tmp1 <- sapply(seq_len(K), FUN = function(k) log(pi_c[k]) +
                  rowSums(sapply(seq_len(ncol(x)), FUN = function(j) rowSums(sapply(seq_len(G),
                    FUN = function(l) t_jg[j, l] * x[, j] * (log(alpha_c[k, l]) +
                      exp_logutilde[, j]) - t_jg[j, l] * mu_i * nu_j[, j] * alpha_c[k,
                      l] * exp_utilde[, j])))))
                # s_ik_tmp1 <- sapply(seq_len(K), FUN = function(k) log(pi_c[k]) +
                # rowSums(sapply(seq_len(ncol(x)), FUN = function(j) rowSums(sapply(seq_len(G),
                # FUN = function(l) t_jg[j, l] * x[, j] * (log(alpha_c[k, l]) + exp_logutilde[,
                # j]) - t_jg[j, l] * mu_i * nu_j[, j] * alpha_c[k, l] * exp_utilde[, j])))))

            } else {
                s_ik_tmp1 <- sapply(seq_len(K), FUN = function(k) log(pi_c[k]) +
                  rowSums(sapply(seq_len(ncol(x)), FUN = function(j) rowSums(sapply(seq_len(G),
                    FUN = function(l) t_jg[j, l] * x[, j] * (log(alpha_c[k, l]) +
                      exp_logutilde[, j]) - t_jg[j, l] * mu_i * nu_j[j] * alpha_c[k,
                      l] * exp_utilde[, j])))))

            }

            s_ik_tmp2 <- s_ik_tmp1 - rowMeans(s_ik_tmp1)
            s_ik_tmp <- apply(s_ik_tmp2, 2, FUN = function(x) exp(x)/rowSums(exp(s_ik_tmp2)))
            if (sum(is.nan(s_ik)) > 0) {
                rmax <- apply(s_ik_tmp1, 1, max)
                s_ik_tmp3 <- s_ik_tmp1 - rmax
                s_ik <- apply(s_ik_tmp3, 2, FUN = function(x) exp(x)/rowSums(exp(s_ik_tmp3)))
            }
            rm(s_ik_tmp1, s_ik_tmp2, s_ik_tmp)
            s_ik <- apply(s_ik, c(1, 2), FUN = function(x) if (x <= 0.5)
                max(tol, x) else if (x > 0.5)
                min(x, 1 - tol))
            s_ik[s_ik == tol] <- tol/(ncol(s_ik) - 1)
            s_ik <- s_ik/rowSums(s_ik)

            # Update pi_c, rho_c, alpha, a ------------------
            pi_c <- colMeans(s_ik)
            # alpha_c
            alpha_c <- alpha_calculation(s_ik = s_ik, t_jg = t_jg, nu_j = nu_j, mu_i = mu_i,
                K = K, G = G, x = x, exp_utilde = exp_utilde)
            # Calculations of exp_utilde and exp_logutilde ------------------
            if (isFALSE(akg)) {
                qu_param <- qu_calculation(s_ik = s_ik, t_jg = t_jg, x = x, mu_i = mu_i,
                  nu_j = nu_j, alpha_c = alpha_c, a = a0)
            } else {
                qu_param <- qukg_calculation(s_ik = s_ik, t_jg = t_jg, x = x, mu_i = mu_i,
                  nu_j = nu_j, alpha_c = alpha_c, a = a0)
            }
            exp_utilde <- qu_param$exp_utilde
            exp_logutilde <- qu_param$exp_logutilde
            # Stopping criteria
            crit_em1 <- sum((sort(alpha_c) - sort(alpha_old1))^2)
            +sum((sort(pi_c) - sort(pi_old))^2)
        }

        # Calculation of t_{c+1} ------------------
        i <- 0
        crit_em2 <- 1
        while (crit_em2 > cvg_lim) {
            i <- i + 1
            rho_old <- rho_c
            alpha_old2 <- alpha_c
            t_old <- t_jg
            if (is.matrix(nu_j)) {
                t_jg_tmp <- sapply(seq_len(G), FUN = function(l) log(rho_c[l]) +
                  rowSums(sapply(seq_len(nrow(x)), FUN = function(i) rowSums(sapply(seq_len(K),
                    FUN = function(k) s_ik[i, k] * x[i, ] * (log(alpha_c[k, l]) +
                      exp_logutilde[i, ]) - s_ik[i, k] * mu_i[i] * nu_j[i, ] * alpha_c[k,
                      l] * exp_utilde[i, ])))))
            } else {
                t_jg_tmp <- sapply(seq_len(G), FUN = function(l) log(rho_c[l]) +
                  rowSums(sapply(seq_len(nrow(x)), FUN = function(i) rowSums(sapply(seq_len(K),
                    FUN = function(k) s_ik[i, k] * x[i, ] * (log(alpha_c[k, l]) +
                      exp_logutilde[i, ]) - s_ik[i, k] * mu_i[i] * nu_j * alpha_c[k,
                      l] * exp_utilde[i, ])))))
            }
            t_jg_tmp2 <- t_jg_tmp - rowMeans(t_jg_tmp)
            t_jg <- apply(t_jg_tmp2, 2, FUN = function(x) exp(x)/rowSums(exp(t_jg_tmp2)))
            if (sum(is.nan(t_jg)) > 0) {
                rmax <- apply(t_jg_tmp, 1, max)
                t_jg_tmp3 <- t_jg_tmp - rmax
                t_jg <- apply(t_jg_tmp3, 2, FUN = function(x) exp(x)/rowSums(exp(t_jg_tmp3)))
            }

            t_jg <- apply(t_jg, c(1, 2), FUN = function(x) if (x <= 0.5)
                max(tol, x) else if (x > 0.5)
                min(x, 1 - tol))
            t_jg[t_jg == tol] <- (tol)/(ncol(t_jg) - 1)

            t_jg <- t_jg/rowSums(t_jg)
            # Update of pi_c, rho_c, alpha, a ------------------------------------
            rho_c <- colMeans(t_jg)

            #--------------- alpha
            alpha_c <- alpha_calculation(s_ik = s_ik, t_jg = t_jg, nu_j = nu_j, mu_i = mu_i,
                K = K, G = G, x = x, exp_utilde = exp_utilde)

            # calculations of exp_utilde and exp_logutilde ------------------
            if (isFALSE(akg)) {
                qu_param <- qu_calculation(s_ik = s_ik, t_jg = t_jg, x = x, mu_i = mu_i,
                  nu_j = nu_j, alpha_c = alpha_c, a = a0)
            } else {
                qu_param <- qukg_calculation(s_ik = s_ik, t_jg = t_jg, x = x, mu_i = mu_i,
                  nu_j = nu_j, alpha_c = alpha_c, a = a0)
            }
            exp_utilde <- qu_param$exp_utilde
            exp_logutilde <- qu_param$exp_logutilde
            crit_em2 <- sum((sort(rho_c) - sort(rho_old))^2)
            +sum((sort(alpha_c) - sort(alpha_old2))^2)
        }

        # Estimation of mu_i ------------------
        if (is.matrix(nu_j)) {
            mu_i <- rowSums(x %*% t_jg)/rowSums(s_ik * (nu_j %*% tcrossprod(t_jg,
                alpha_c)))
        } else {
            mu_i <- rowSums(x %*% t_jg)/rowSums(s_ik * rowSums(matrix(nrow = K, sapply(seq_len(G),
                FUN = function(l) alpha_c[, l] * colSums(t_jg * nu_j)[l]))))
        }
        # Estimation of a and update of exp_utilde, exp_logutilde, alpha_c
        # --------------------
        lb_old <- lb
        a_old <- a0
        if (is.null(a)) {
            if (isFALSE(akg)) {
                left_bound <- sum(exp_logutilde)
                right_bound <- sum(exp_utilde)
                a0 <- stats::uniroot(f = foo_a, lower = 0.01, upper = n * m, left_bound = left_bound,
                  right_bound = right_bound, nb = n * m)$root
                # a0 <- dicho(x = 0.01, y = abs(max(left_bound, right_bound)), threshold = 1e-08,
                # nb = n * m, left_bound = left_bound, right_bound = right_bound)
            } else {
                left_bound <- crossprod(s_ik, exp_logutilde %*% t_jg)
                right_bound <- crossprod(s_ik, exp_utilde %*% t_jg)
                n_kg <- crossprod(s_ik, matrix(nrow = n, ncol = m, 1) %*% t_jg)
                # a0 <- matrix(nrow = K, ncol = G, sapply(1:(K*G), FUN = function(g) dicho(x =
                # 0.01, y = 100, threshold = 1e-08, nb = n_kg[g], left_bound = left_bound[g],
                # right_bound = right_bound[g])))
                a0 <- matrix(nrow = K, ncol = G, sapply(1:(K * G), FUN = function(g) stats::uniroot(f = foo_a,
                  lower = 0.01, upper = n_kg[g], left_bound = left_bound[g], right_bound = right_bound[g],
                  nb = n_kg[g])$root))
            }
        }
        if (isFALSE(akg)) {
            qu_param <- qu_calculation(s_ik = s_ik, t_jg = t_jg, x = x, mu_i = mu_i,
                nu_j = nu_j, alpha_c = alpha_c, a = a0)
        } else {
            qu_param <- qukg_calculation(s_ik = s_ik, t_jg = t_jg, x = x, mu_i = mu_i,
                nu_j = nu_j, alpha_c = alpha_c, a = a0)
        }
        alpha_c <- alpha_calculation(s_ik = s_ik, t_jg = t_jg, nu_j = nu_j, mu_i = mu_i,
            K = K, G = G, x = x, exp_utilde = qu_param$exp_utilde)
        # Calculation of the lower bound ----------------------------------------
        lb_out <- lb_calculation(x = x, qu_param = qu_param, s_ik = s_ik, pi_c = pi_c,
            t_jg = t_jg, rho_c = rho_c, mu_i = mu_i, nu_j = nu_j, alpha_c = alpha_c,
            a = a0, akg = akg)

        lb <- lb_out$lb
        # Stopping criterion based on the lower bound --------------------
        if (j > 1) {
            crit <- abs((lb - lb_old)/lb_old)
        }
        lbtt <- c(lbtt, lb)
    }
    colclass <- apply(t_jg, 1, which.max)
    rowclass <- apply(s_ik, 1, which.max)
    # Output
    strategy <- list(akg = akg, cvg_lim = cvg_lim)
    parameters <- list(alpha = alpha_c, pi = pi_c, rho = rho_c, mu_i = mu_i, nu_j = nu_j,
        a = a0)
    info <- list(s_ik = s_ik, t_jg = t_jg, exp_logutilde = qu_param$exp_logutilde,
        exp_utilde = qu_param$exp_utilde, lb = lb, ent_ZW = lb_out$ent_ZW, nbiter = j,
        lbtt = lbtt, a_tilde = qu_param$a_tilde, b_tilde = qu_param$b_tilde)
    output <- list(data = x, K = K, G = G, classification = list(rowclass = rowclass,
        colclass = colclass), strategy = strategy, parameters = parameters, info = info)
    class(output) <- append(class(output), "cobiclustering")
    return(output)
}
