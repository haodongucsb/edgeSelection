#' @title Edge Selection
#' @description  Extension of gss package for edge selection
#' @export
#' @param data Data frame containing all variables.
#' @param method Method type to select edges.
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
#' \item \code{rho} Method to construct rho function for neighborhood selection method.
#' \item \code{ydomain} Data frame specifying marginal support of conditional density in the neighborhood selection method.
#' \item \code{yquad} Quadrature for calculating integral on Y domain in the neighborhood selection method. Mandatory if response variables other than factors or numerical vectors are involved.
#' \item \code{skip.iter} Flag indicating whether to use initial values of theta and skip theta iteration in the neighborhood selection method.
#' \item \code{W} Optional matrix to specify weights for two-way interactions in the neighborhood selection method for each node.
#' \item \code{neighborhoodMethod} Method type in the neighborhood selection method to select tuning parameter which controls sparsity of the graph.
#' \item \code{px} Dimension of variables in the semi-parametric method to be estimated using non-parametric method.
#' \item \code{maxLambda} Number used to generate the range of the tuning parameter for selection of Lambda matrix.
#' \item \code{maxTheta} Number used to generate the range of the tuning parameter for selection of Theta matrix.
#' \item \code{iterLambda} Number of iterations to find optimal Lambda matrix.
#' \item \code{iterTheta} Number of iterations to find optimal Theta matrix.
#' \item \code{cutoff} Cutoff value for squared projection.
#' \item \code{N} Number of simulations in non-parametric part to calculate the standard error.
#' \item \code{semiMethod} Method type to select optimal matrix in semi-parametric method.
#' }
#' @details
#' \code{type, alpha, subset, na.action, seed, prec, maxiter, id.basis, nbasis} are
#' arguments shared by the joint and neighborhood selection method. They also work the same as in \code{gss} package.
#' \code{domain, quad} are two arguments for the joint method and work the same as \code{ssden1} in \code{gss} package.
#' \code{w} is an option argument in the joint method.
#' \code{rho, ydomain, yquad, skip.iter} are arguments for the neighborhood selection method. They work the same as \code{sscden1} in \code{gss} package.
#' \code{W2, neighborhoodMethod} are two optional arguments in the neighborhood selection method.
#' The rest of argument options are specifically for the semi-parametric method.
#' @import stats
#' @import gss
#' @import MASS
#' @import pracma
#' @import matrixcalc
#' @import QUIC
#' @import glasso
#' @import quadprog
#' @usage
#' edge.selection(data, method, ...)
#' @examples
#' # Use joint method for edge selection.
#' library(gss)
#' data(NO2)
#' edge.selection(data = NO2, method = "joint", nbasis = 100)
#' # Use neighborhood selection method for edge selection.
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
#' edge.selection(data = data, method = "neighborhood")
#' # Use semi-parametric method for edge selection.
#' # Assume we estimate the first five variables using non-parametric method.
#' px <- 5
#' edge.selection(data = data, method ="semi", px = px)

edge.selection <- function(data, method = c("joint", "neighborhood", "semi"), ...) {
  method <- match.arg(method)
  edgeMatrix <- switch(method,
    joint = selection.joint(
      data, ...
    ),
    neighborhood = selection.neighborhood(
      data, ...
    ),
    semi = selection.semi(data, ...)
  )
  edgeMatrix
}

My_solve.QP <- function(Dmat, dvec, Amat, bvec) {
  solution <- tryCatch(solve.QP(Dmat, dvec, Amat, bvec)$solution, error = function(x) NA)
  if (is.na(solution[1])) {
    M <- solve(Dmat)
    Dmat <- t(M) %*% M
    sc <- norm(Dmat, "2")
    solution <- tryCatch(solve.QP(Dmat = Dmat / sc, dvec = dvec / sc, Amat = Amat, bvec = bvec, meq = 0, factorized = FALSE)$solution, error = function(x) NA)
    if (is.na(solution[1])) {
      Dmat <- diag(diag(Dmat))
      solution <- solve.QP(Dmat, dvec, Amat, bvec)$solution
    }
  }
  return(solution)
}
