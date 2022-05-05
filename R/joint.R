#' Joint method for edge selection
#' @description Apply joint selection method by considering joint distribution of all variables.
#' @export
#' @param data Data frame
#' @param ... Any options can be defined.
#' \itemize{
#' \item \code{type} List specifying the type of spline for each variable.
#' \item \code{alpha} Parameter defining cross-validation score for smoothing parameter selection.
#' \item \code{subset} Optional vector specifying a subset of observations to be used in the fitting process.
#' \item \code{na.action} Function which indicates what should happen when the data contain NAs.
#' \item \code{seed} Seed to be used for the random generation of "knots."
#' \item \code{prec} Precision requirement for internal iterations.
#' \item \code{maxiter} Maximum number of iterations allowed for internal iterations.
#' \item \code{id.basis} Index of observations to be used as "knots."
#' \item \code{nbasis} Number of "knots" to be used.
#' \item \code{domain} Data frame specifying marginal support of density in the joint method.
#' \item \code{quad} Quadrature for calculating integral in the joint method. Mandatory if variables other than factors or numerical vectors are involved.
#' \item \code{w} Optional vector to specify weights for two-way interactions in the joint method.
#' }
#' @usage
#' selection.joint(data, ...)
#' @examples
#' library(gss)
#' data(NO2)
#' edge.selection(data = NO2, family = "joint", nbasis = 100)

selection.joint <- function(data, ...) {
  n <- dim(data)[1]
  p <- dim(data)[2]
  type <- NULL
  alpha <- 1.4
  subset <- NULL
  na.action <- na.omit
  seed <- 5732
  prec <- 1e-7
  maxiter <- 30
  id.basis <- NULL
  nbasis <- NULL
  domain <- as.list(NULL)
  quad <- NULL
  w <- NULL
  params <- list(...)
  for (name in names(params)) {
    assign(name, params[[name]])
  }
  if (is.null(subset)) {
    subset <- 1:n
  }
  if (p > 12) {
    stop("dimension is too high, try neighborhood selection")
  }
  formula <- names(data)[1]
  for (i in 2:p) {
    formula <- paste(formula, "+", names(data)[i], sep = "")
  }
  formula <- paste("~", "(", formula, ")", "^2", sep = "")
  if (is.null(w)) {
    w <- rep(1, p * (p - 1) / 2)
  } else {
    w <- w
  }
  if (is.null(nbasis) & is.null(id.basis)) {
    densityfit0 <- ssden2(
      formula = formula, data = data, type = type, alpha = alpha, subset = subset,
      na.action = na.action, w = w, domain = domain, quad = quad, seed = seed, prec = prec, maxiter = maxiter
    )
  } else if (is.null(id.basis)) {
    densityfit0 <- ssden2(
      formula = formula, data = data, type = type, alpha = alpha, subset = subset,
      na.action = na.action, w = w, domain = domain, quad = quad, seed = seed, prec = prec, maxiter = maxiter, nbasis = nbasis
    )
  } else {
    densityfit0 <- ssden2(
      formula = formula, data = data, type = type, alpha = alpha, subset = subset,
      na.action = na.action, w = w, domain = domain, quad = quad, prec = prec, maxiter = maxiter, id.basis = id.basis
    )
  }
  theta2 <- rep(1, (p - 1) * p / 2)
  ind <- densityfit0$id.basis
  for (i in 1:((p - 1) * p / 2)) {
    theta2[i] <- sqrt(sum((densityfit0$theta2[i] * densityfit0$r[, , i + p] %*% densityfit0$c)^2))
  }
  densityfit <- ssden2(
    formula = formula, data = data, type = type, alpha = alpha, subset = subset,
    na.action = na.action, w = w, domain = domain, quad = quad, seed = seed, prec = prec, maxiter = maxiter, theta2 = theta2, id.basis = ind
  )
  theta1.1 <- densityfit$theta1
  theta2.1 <- densityfit$theta2
  ind0 <- densityfit$id.basis
  ltheta2 <- length(theta2.1)
  c1 <- densityfit$c
  d1 <- densityfit$d
  U <- densityfit$r
  U2 <- NULL
  for (i in (p + 1):(p * (p + 1) / 2)) {
    U2 <- cbind(U2, U[, , i])
  }
  u <- densityfit$rbasis
  u2 <- NULL
  for (i in (p + 1):(p * (p + 1) / 2)) {
    u2 <- cbind(u2, u[, , i])
  }
  d1 <- densityfit$d
  g <- densityfit$g
  t1 <- rep(0, length(g))
  for (i in 1:length(g)) {
    t1[i] <- exp(-g[i]) / sum(exp(-g))
  }
  B2 <- densityfit$int$r[, (p + 1):(p * (p + 1) / 2)]
  la1 <- 10^densityfit$lambda / 2
  Hpart1 <- U2 %*% kronecker(diag(ltheta2), c1)
  Gmatrix <- -t(Hpart1) %*% t1 + t(B2) %*% c1 + la1 * kronecker(diag(ltheta2), t(c1)) %*% t(u2) %*% c1
  Hpart2 <- t(t1) %*% U2 %*% kronecker(diag(ltheta2), c1)
  Hpart3 <- sqrt(t1) * Hpart1
  Hmatrix <- t(Hpart3) %*% Hpart3 - t(Hpart2) %*% Hpart2
  dvec <- (Hmatrix %*% theta2.1 - Gmatrix)
  Dmat <- Hmatrix
  Amat <- diag(ltheta2)
  bvec <- rep(0, ltheta2)
  starttime <- Sys.time()
  theta2 <- My_solve.QP(Dmat, dvec, t(Amat), bvec)
  # the range of M to be chosen from
  prop <- c(0.6, 0.65, 0.7, 0.75, 0.8)
  M <- sum(theta2)
  obj <- edgedect(
    data = data, formula = formula, type = type, alpha = alpha, subset = subset, na.action = na.action, w = w, domain = domain, quad = quad, prec = prec, maxiter = maxiter,
    M0 = M, M_list = prop, maxiteration = 10, tolerance = 1e-4, id.basis = ind0, theta2 = theta2, n = n, dimen = p
  )
  estimatedPara <- obj$theta2
  edgeMatrix <- diag(1, p)
  count <- 1
  for (l in 1:(p - 1)) {
    for (k in (l + 1):p) {
      if (estimatedPara[count] > 0) {
        edgeMatrix[l, k] <- estimatedPara[count]
        edgeMatrix[k, l] <- estimatedPara[count]
      }
      count <- count + 1
    }
  }
  edgeMatrix
}
par_select_M <- function(densityfit, w, M0, M_list, k = 5, n) {
  M_list1 <- M0 * M_list
  cv_mat <- rep(0, length(M_list1))
  count <- 0
  densityfit <- densityfit
  theta1.1 <- densityfit$theta1
  p <- length(theta1.1)
  theta2.1 <- densityfit$theta2
  ltheta2 <- length(theta2.1)
  c1 <- densityfit$c
  d1 <- densityfit$d
  R <- densityfit$R
  ltheta2 <- length(theta2.1)
  U <- densityfit$r
  U2 <- NULL
  for (i in (p + 1):(p * (p + 1) / 2)) {
    U2 <- cbind(U2, U[, , i])
  }
  u <- densityfit$rbasis
  u2 <- NULL
  for (i in (p + 1):(p * (p + 1) / 2)) {
    u2 <- cbind(u2, u[, , i])
  }
  g <- densityfit$g
  for (M in M_list1) {
    count <- count + 1
    loglikFold <- matrix(NA, ncol = k, nrow = 1)
    grps <- cut(1:n, k, labels = FALSE)[sample(n)]
    for (kth in 1:k) {
      omit <- which(grps == kth)
      t2 <- rep(0, length(g))
      gomit <- g[omit]
      gkeep <- g[-omit]
      for (i in 1:length(g)) {
        t2[i] <- exp(-g[i]) / sum(exp(-gkeep))
      }
      for (i in omit) {
        t2[i] <- 0
      }
      B2 <- densityfit$int$r[, (p + 1):(ltheta2 + p)]
      la1 <- 10^densityfit$lambda / 2
      Gmatrix <- -kronecker(diag(ltheta2), t(c1)) %*% t(U2) %*% t2 + t(B2) %*% c1 + la1 * kronecker(diag(ltheta2), t(c1)) %*% t(u2) %*% c1
      Hpart1 <- U2 %*% kronecker(diag(ltheta2), c1)
      Hpart2 <- t(t2) %*% U2 %*% kronecker(diag(ltheta2), c1)
      Hpart3 <- sqrt(t2) * Hpart1
      Hmatrix <- (t(Hpart3) %*% Hpart3 - t(Hpart2) %*% Hpart2)
      dvec <- Hmatrix %*% theta2.1 - Gmatrix
      Dmat <- Hmatrix
      Amat <- rbind(diag(ltheta2), -w)
      bvec <- c(rep(0, ltheta2), -M)
      theta2 <- My_solve.QP(Dmat, dvec, t(Amat), bvec)
      theta2[theta2 < 1e-8] <- 0
      r <- densityfit$r
      R <- densityfit$R1
      for (i in (p + 1):(p * (p + 1) / 2)) {
        R <- R + theta2[i - p] * r[, , i]
      }
      cv_g <- densityfit$s %*% d1 + R %*% c1
      cv_g <- cv_g[omit]
      score <- (sum(exp(-cv_g)))
      loglikFold[, kth] <- score
    }
    cv_mat[count] <- log((apply(loglikFold, 1, sum)) / n)
  }
  min_error_pos <- which.min(cv_mat)
  M_opt <- M_list1[min_error_pos]
  list(M_opt = M_opt)
}

edgedect <- function(data, formula, alpha = 1.4, subset = NULL, na.action = na.omit, w = NULL, type = NULL,
                     domain = as.list(NULL), quad = NULL, prec = 1e-7, maxiter = 30, M0, M_list, maxiteration, tolerance, id.basis = NULL, theta2, n, dimen) {
  p <- dimen
  loop <- 1
  ltheta2 <- length(theta2)
  Theta2 <- matrix(0, ltheta2, maxiteration)
  Theta2[, loop] <- theta2
  l <- rep(0, maxiteration)
  l[loop] <- 0
  diff2 <- 1
  d <- rep(0, maxiteration)
  d[1] <- diff2
  zeros <- (theta2 == 0)
  zerodiff <- 1
  while (diff2 >= tolerance & loop < maxiteration & (zerodiff) >= 1) {
    loop <- loop + 1

    densityfit <- ssden2(
      formula = formula, data = data, id.basis = id.basis, type = type, alpha = alpha, subset = subset,
      na.action = na.action, domain = domain, quad = quad, maxiter = maxiter, theta2 = theta2, w = w
    )
    M <- par_select_M(densityfit, w, M0, M_list, k = 5, n)[[1]]
    theta1.1 <- densityfit$theta1
    theta2.1 <- densityfit$theta2
    R <- densityfit$R
    ltheta2 <- length(theta2.1)
    c1 <- densityfit$c
    d1 <- densityfit$d
    U <- densityfit$r
    U2 <- NULL
    for (i in (p + 1):(p * (p + 1) / 2)) {
      U2 <- cbind(U2, U[, , i])
    }
    u <- densityfit$rbasis
    u2 <- NULL
    for (i in (p + 1):(p * (p + 1) / 2)) {
      u2 <- cbind(u2, u[, , i])
    }
    g <- densityfit$g
    t2 <- rep(0, length(g))
    for (i in 1:length(g)) {
      t2[i] <- exp(-g[i]) / sum(exp(-g))
    }
    B2 <- densityfit$int$r[, (p + 1):(p * (p + 1) / 2)]
    la1 <- 10^densityfit$lambda / 2
    Gmatrix <- -kronecker(diag(ltheta2), t(c1)) %*% t(U2) %*% t2 + t(B2) %*% c1 + la1 * kronecker(diag(ltheta2), t(c1)) %*% t(u2) %*% c1
    Hpart1 <- U2 %*% kronecker(diag(ltheta2), c1)
    Hpart2 <- t(t2) %*% U2 %*% kronecker(diag(ltheta2), c1)
    Hpart3 <- sqrt(t2) * Hpart1
    Hmatrix <- (t(Hpart3) %*% Hpart3 - t(Hpart2) %*% Hpart2)
    dvec <- Hmatrix %*% theta2.1 - Gmatrix
    Dmat <- Hmatrix
    Amat <- rbind(diag(ltheta2), -w)
    bvec <- c(rep(0, ltheta2), -M)
    theta2 <- My_solve.QP(Dmat, dvec, t(Amat), bvec)
    theta2[theta2 < 1e-8] <- 0
    theta2[zeros] <- 0
    zeros <- (theta2 == 0)
    l[loop] <- la1
    Theta2[, loop] <- theta2
    if (loop <= 3) {
      zerodiff <- 1
    } else {
      zerodiff <- sum(theta2 > 0) - sum(Theta2[, loop - 3] > 0)
    }
    diff2 <- sqrt(sum((Theta2[, loop] - Theta2[, loop - 1])^2)) / (sqrt(sum((Theta2[, loop - 1])^2)) + 1e-6)
    d[loop] <- diff2
  }
  thre_list <- c(1e-1, 1e-2, 1e-3)
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
    for (i in (p + 1):(p * (p + 1) / 2)) {
      R0 <- R0 + theta20[i - p] * densityfit$r[, , i] / w[i - p]
    }
    cv_g <- densityfit$s %*% d1 + R0 %*% c1
    score <- (sum(exp(-cv_g)))
    thre_mat[count] <- log(score / n)
  }
  min_pos <- which.min(thre_mat)
  bound <- thre_list[min_pos]
  theta2[theta2 < bound] <- 0
  list(loop = loop, theta2 = theta2, densityfit = densityfit)
}
