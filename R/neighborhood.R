#' @title Neighborhood selection method
#' @description Apply neighborhood selection method by considering conditional density function on each variable.
#' @import doParallel
#' @import foreach
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
#' \item \code{rho} Method to construct rho function for neighborhood selection method.
#' \item \code{ydomain} Data frame specifying marginal support of conditional density in the neighborhood selection method.
#' \item \code{yquad} Quadrature for calculating integral on Y domain in the neighborhood selection method. Mandatory if response variables other than factors or numerical vectors are involved.
#' \item \code{skip.iter} Flag indicating whether to use initial values of theta and skip theta iteration in the neighborhood selection method.
#' \item \code{W} Optional matrix to specify weights for two-way interactions in the neighborhood selection method for each node.
#' \item \code{neighborhoodMethod} Method type in the neighborhood selection method to select tuning parameter which controls sparsity of the graph.
#' }
#' @usage
#' selection.neighborhood(data, ...)
#' @examples
#' # Parallel backend
#' library(doMC)
#' library(foreach)
#' library(huge)
#' registerDoMC(20)
#' n <- 200; p <- 20
#' # Simulate high dimension data
#' set.seed(5732)
#' z <- huge.generator(n, d = p, graph = "random", prob = .2, verbose = FALSE, vis = FALSE, v = .65)
#' data <- data.frame(z$data)
#' edge.selection(data = data, family = "neighborhood")
selection.neighborhood <- function(data, ...) {
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
  rho <- list("xy")
  ydomain <- as.list(NULL)
  yquad <- NULL
  skip.iter <- TRUE
  W <- NULL
  neighborhoodMethod <- "cv"
  params <- list(...)
  for (name in names(params)) {
    assign(name, params[[name]])
  }
  formula <- names(data)[1]
  for (i in 2:p) {
    formula <- paste(formula, "+", names(data)[i], sep = "")
  }
  estimated <- foreach(x = 1:p, .combine = cbind) %dopar% {
    if (is.null(W)) {
      w2 <- rep(1, p - 1)
    } else {
      w2 <- W[x, ]
    }
    formula <- paste("~", "(", formula, "-", colnames(data)[x], ")", "*", colnames(data)[x], sep = "")
    formula <- eval(parse(text = formula))
    response <- paste("~", colnames(data)[x], sep = "")
    response <- eval(parse(text = response))
    if (is.null(nbasis) & is.null(id.basis)) {
      densityfit0 <- sscden0(
        formula = formula, response = response, data = data, type = type, alpha = alpha, subset = subset,
        na.action = na.action, ydomain = ydomain, yquad = yquad, seed = seed, prec = prec, maxiter = maxiter, skip.iter = skip.iter
      )
    } else if (is.null(id.basis)) {
      densityfit0 <- sscden0(
        formula = formula, response = response, data = data, type = type, alpha = alpha, subset = subset,
        na.action = na.action, ydomain = ydomain, yquad = yquad, seed = seed, prec = prec, maxiter = maxiter, nbasis = nbasis, skip.iter = skip.iter
      )
    } else {
      densityfit0 <- sscden0(
        formula = formula, response = response, data = data, type = type, alpha = alpha, subset = subset,
        na.action = na.action, ydomain = ydomain, yquad = yquad, prec = prec, maxiter = maxiter, id.basis = id.basis, skip.iter = skip.iter
      )
    }
    theta0 <- rep(1, p - 1)
    for (i in 2:p) {
      theta0[i - 1] <- sqrt(sum((10^densityfit0$theta[i] * densityfit0$r[, , i] %*% densityfit0$c)^2))
    }
    if (is.null(nbasis) & is.null(id.basis)) {
      densityfit <- sscden_selection(
        formula = formula, response = response, data = data, w2 = w2, type = type, alpha = alpha, subset = subset,
        na.action = na.action, rho = rho, ydomain = ydomain, yquad = yquad, seed = seed, prec = prec, maxiter = maxiter, skip.iter = skip.iter, p = p, theta2 = theta0
      )
    } else if (is.null(id.basis)) {
      densityfit <- sscden_selection(
        formula = formula, response = response, data = data, w2 = w2, type = type, alpha = alpha, subset = subset,
        na.action = na.action, rho = rho, ydomain = ydomain, yquad = yquad, seed = seed, prec = prec, maxiter = maxiter, nbasis = nbasis, skip.iter = skip.iter, p = p, theta2 = theta0
      )
    } else {
      densityfit <- sscden_selection(
        formula = formula, response = response, data = data, w2 = w2, type = type, alpha = alpha, subset = subset,
        na.action = na.action, rho = rho, ydomain = ydomain, yquad = yquad, prec = prec, maxiter = maxiter, id.basis = id.basis, skip.iter = skip.iter, p = p, theta2 = theta0
      )
    }
    theta1 <- densityfit$theta1
    theta2.1 <- densityfit$theta2
    ind0 <- densityfit$id.basis
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
    t1 <- rep(0, length(g))
    for (i in 1:length(g)) {
      t1[i] <- exp(-g[i]) / sum(exp(-g))
    }
    B2 <- densityfit$int.r[, (p + 1):(2 * p - 1)]
    la1 <- 10^densityfit$lambda / 2
    Hpart1 <- U2 %*% kronecker(diag(ltheta2), c1)
    Gmatrix <- -t(Hpart1) %*% t1 + t(B2) %*% c1 + la1 * kronecker(diag(ltheta2), t(c1)) %*% t(u2) %*% c1
    Hpart2 <- t(t1) %*% U2 %*% kronecker(diag(ltheta2), c1)
    Hpart3 <- sqrt(t1) * Hpart1
    Hmatrix <- t(Hpart3) %*% Hpart3 - t(Hpart2) %*% Hpart2
    if (is.positive.definite(Hmatrix)) {
      Hmatrix <- Hmatrix
    } else {
      Hmatrix <- Hmatrix + diag(rep(1e-6, dim(Hmatrix)[1]))
    }
    dvec <- (Hmatrix %*% theta2.1 - Gmatrix)
    Dmat <- Hmatrix
    Amat <- diag(ltheta2)
    bvec <- rep(0, ltheta2)
    starttime <- Sys.time()
    theta2 <- My_solve.QP(Dmat, dvec, t(Amat), bvec)
    # the range of M to be chosen from
    prop <- c(0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
    M <- sum(theta2)
    if (neighborhoodMethod == "cv") {
      obj <- cond_select(
        data = data, formula = formula, response = response, w2 = w2, alpha = alpha, subset = subset, na.action = na.action, rho = rho, ydomain = ydomain, yquad = yquad,
        prec = prec, maxiter = maxiter, skip.iter = skip.iter, M0 = M, M_list = prop, maxiteration = 10, tolerance = 1e-4, id.basis = ind0, theta2 = theta2, n = n, p = p
      )
    } else {
      obj <- cond_BIC(
        data = data, formula = formula, response = response, w2 = w2, alpha = alpha, subset = subset, na.action = na.action, rho = rho, ydomain = ydomain, yquad = yquad,
        prec = prec, maxiter = maxiter, skip.iter = skip.iter, M0 = M, M_list = prop, maxiteration = 10, tolerance = 1e-4, id.basis = ind0, theta2 = theta2, n = n, p = p
      )
    }
    estimatedPara <- obj$theta2
    append(estimatedPara, 0, x - 1)
  }

  estimatedPara <- as.matrix(estimated)
  edgeMatrix <- diag(1, p)
  for (l in 1:p) {
    for (k in 1:p) {
      if (estimatedPara[l, k] > 0 & estimatedPara[k, l] > 0) {
        edgeMatrix[l, k] <- estimatedPara[l, k]
        edgeMatrix[k, l] <- estimatedPara[k, l]
      }
    }
  }
  edgeMatrix
}
