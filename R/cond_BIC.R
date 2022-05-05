#' Neighborhood selection method with BIC method to select the tuning parameter.
#' @import gss
#' @export
#' @param data Data frame
#' @param formula Symbolic description of the model to be fit.
#' @param response Formula listing response variables.
#' @param type List specifying the type of spline for each variable.
#' @param alpha Parameter defining cross-validation score for smoothing parameter selection.
#' @param subset Optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action Function which indicates what should happen when the data contain NAs.
#' @param rho Method to construct rho function for neighborhood selection method.
#' @param ydomain Data frame specifying marginal support of conditional density in the neighborhood selection method.
#' @param yquad Quadrature for calculating integral on Y domain in the neighborhood selection method. Mandatory if response variables other than factors or numerical vectors are involved.
#' @param prec Precision requirement for internal iterations.
#' @param maxiter Maximum number of iterations allowed for internal iterations.
#' @param skip.iter Flag indicating whether to use initial values of theta and skip theta iteration in the neighborhood selection method.
#' @param M0 Upper bound
#' @param M_list List of values for tuning parameter selection
#' @param maxiteration Max number of iteration
#' @param tolerance Threshold for convergence
#' @param id.basis Index of observations to be used as "knots."
#' @param theta2 Parameters for edge selection.
#' @param w2 Optional vector to specify weights for two-way interactions
#' @param n Number of observation
#' @param p Dimension of data frame
cond_BIC <- function(data, formula, response, type = NULL, alpha, subset, na.action, rho, ydomain, yquad, prec, maxiter, skip.iter,
                     M0, M_list, maxiteration, tolerance, id.basis = NULL, theta2, w2 = NULL, n, p) {
  loop <- 1
  ltheta2 <- length(theta2)
  Theta2 <- matrix(0, ltheta2, maxiteration)
  Theta2[, loop] <- theta2
  l <- rep(0, maxiteration)
  l[loop] <- 0
  diff2 <- 1
  d <- rep(0, maxiteration)
  d[1] <- diff2
  zerodiff <- 1
  check_na <- 0
  while (diff2 >= tolerance & loop < maxiteration & (zerodiff) >= 1) {
    loop <- loop + 1
    densityfit <- sscden_selection(
      formula = formula, response = response, data = data, w2 = w2, type = type, alpha = alpha, subset = subset,
      na.action = na.action, rho = rho, ydomain = ydomain, yquad = yquad, prec = prec, maxiter = maxiter, skip.iter = skip.iter, id.basis = id.basis, p = p, theta2 = theta2
    )
    M <- cond_select_M(densityfit, w2, M0, M_list, k = 5, n, p = p)[[1]]
    theta1 <- densityfit$theta1
    theta2.1 <- densityfit$theta2
    R <- densityfit$R
    ltheta2 <- length(theta2.1)
    c1 <- densityfit$c
    d1 <- densityfit$d
    U <- densityfit$r
    U2 <- NULL
    for (i in (p + 1):(2 * p - 1)) {
      U2 <- cbind(U2, U[, , i])
    }
    u <- densityfit$rbasis
    u2 <- NULL
    for (i in (p + 1):(2 * p - 1)) {
      u2 <- cbind(u2, u[, , i])
    }
    g <- densityfit$g
    t2 <- rep(0, length(g))
    for (i in 1:length(g)) {
      t2[i] <- exp(-g[i]) / sum(exp(-g))
    }
    B2 <- densityfit$int.r[, (p + 1):(2 * p - 1)]
    la1 <- 10^densityfit$lambda / 2
    Hpart1 <- U2 %*% kronecker(diag(ltheta2), c1)
    Gmatrix <- -t(Hpart1) %*% t2 + t(B2) %*% c1 + la1 * kronecker(diag(ltheta2), t(c1)) %*% t(u2) %*% c1
    Hpart2 <- t(t2) %*% U2 %*% kronecker(diag(ltheta2), c1)
    Hpart3 <- sqrt(t2) * Hpart1
    Hmatrix <- t(Hpart3) %*% Hpart3 - t(Hpart2) %*% Hpart2
    if (is.positive.definite(Hmatrix)) {
      Hmatrix <- Hmatrix
    } else {
      Hmatrix <- Hmatrix + diag(rep(1e-6, dim(Hmatrix)[1]))
    }
    dvec <- Hmatrix %*% theta2.1 - Gmatrix
    Dmat <- Hmatrix
    Amat <- rbind(diag(ltheta2), -w2)
    bvec <- c(rep(0, ltheta2), -M)
    theta2 <- My_solve.QP(Dmat, dvec, t(Amat), bvec)
    theta2[theta2 < 1e-10] <- 0
    l[loop] <- la1
    Theta2[, loop] <- theta2
    if (loop <= 3) {
      zerodiff <- 1
    } else {
      zerodiff <- sum(theta2 > 0) - sum(Theta2[, loop - 3] > 0)
    }
    diff2 <- sqrt(sum((Theta2[, loop] - Theta2[, loop - 1])^2)) / (sqrt(sum((Theta2[, loop - 1])^2)) + 1e-6)
    # diff_c <- sqrt(sum((c_matrix[, loop] - c_matrix[, loop - 1])^2)) / (sqrt(sum((c_matrix[, loop - 1])^2)) + 1e-6)
    d[loop] <- diff2
  }

  thre_list <- c(1e-2, 1e-3, 1e-4)
  count <- 0
  thre_mat <- rep(0, length(thre_list))
  for (thre in thre_list) {
    count <- count + 1
    theta20 <- theta2
    theta20[theta20 < thre] <- 0
    R0 <- matrix(0, dim(densityfit$r)[1], dim(densityfit$r)[2])
    for (i in 1:p) {
      R0 <- R0 + 10^densityfit$theta1[i] * densityfit$r[, , i]
    }
    for (i in (p + 1):(2 * p - 1)) {
      R0 <- R0 + theta20[i - p] * densityfit$r[, , i] / w2[i - p]
    }
    cv_g <- densityfit$s %*% d1 + R0 %*% c1
    score <- (sum(exp(-cv_g)))
    thre_mat[count] <- log(score / n) + densityfit$loss_int
  }
  min_pos <- which.min(thre_mat)
  bound <- thre_list[min_pos]
  theta2[theta2 < bound] <- 0
  list(loop = loop, theta2 = theta2, densityfit = densityfit)
}


cond_select_M <- function(densityfit, w2, M0, M_list, k = 1, n, p) {
  M_list1 <- M0 * M_list
  cv_mat <- rep(0, length(M_list1))
  count <- 0
  densityfit <- densityfit
  theta1.1 <- densityfit$theta1
  theta2.1 <- densityfit$theta2
  ltheta2 <- length(theta2.1)
  c1 <- densityfit$c
  d1 <- densityfit$d
  R <- densityfit$R
  ltheta2 <- length(theta2.1)
  U <- densityfit$r
  int.r <- densityfit$int.r
  U2 <- NULL
  for (i in (p + 1):(2 * p - 1)) {
    U2 <- cbind(U2, U[, , i])
  }
  u <- densityfit$rbasis
  u2 <- NULL
  for (i in (p + 1):(2 * p - 1)) {
    u2 <- cbind(u2, u[, , i])
  }
  g <- densityfit$g
  t2 <- rep(0, length(g))
  for (i in 1:length(g)) {
    t2[i] <- exp(-g[i]) / sum(exp(-g))
  }
  B2 <- densityfit$int.r[, (p + 1):(ltheta2 + p)]
  la1 <- 10^densityfit$lambda / 2
  Hpart1 <- U2 %*% kronecker(diag(ltheta2), c1)
  Gmatrix <- -t(Hpart1) %*% t2 + t(B2) %*% c1 + la1 * kronecker(diag(ltheta2), t(c1)) %*% t(u2) %*% c1
  Hpart2 <- t(t2) %*% U2 %*% kronecker(diag(ltheta2), c1)
  Hpart3 <- sqrt(t2) * Hpart1
  Hmatrix <- t(Hpart3) %*% Hpart3 - t(Hpart2) %*% Hpart2
  if (is.positive.definite(Hmatrix)) {
    Hmatrix <- Hmatrix
  } else {
    Hmatrix <- Hmatrix + diag(rep(1e-6, dim(Hmatrix)[1]))
  }
  dvec <- Hmatrix %*% theta2.1 - Gmatrix
  Dmat <- Hmatrix
  Amat <- rbind(diag(ltheta2), -w2)
  for (M in M_list1) {
    count <- count + 1
    bvec <- c(rep(0, ltheta2), -M)
    theta2 <- My_solve.QP(Dmat, dvec, t(Amat), bvec)
    theta2[theta2 < 1e-10] <- 0
    int.r.wk <- 0
    for (i in 1:p) {
      int.r.wk <- int.r.wk + 10^theta1.1[i] * int.r[, i]
    }
    for (i in (p + 1):(2 * p - 1)) {
      int.r.wk <- int.r.wk + theta2[i - p] * int.r[, i] / w2[i - p]
    }
    r <- densityfit$r
    R0 <- densityfit$R1
    for (i in (p + 1):(2 * p - 1)) {
      R0 <- R0 + theta2[i - p] * r[, , i] / w2[i - p]
    }
    cv_g <- densityfit$s %*% d1 + R0 %*% c1
    score <- (sum(exp(-cv_g)))
    cv_mat[count] <- log(score / n) + (dot(int.r.wk, c1) + dot(densityfit$int.s, d1)) / n + log(n * sum(theta2 > 1e-10))
  }
  min_error_pos <- which.min(cv_mat)
  M_opt <- M_list1[min_error_pos]
  list(M_opt = M_opt)
}
