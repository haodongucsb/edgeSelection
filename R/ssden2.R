#' Estimating joint density using smoothing splines
#' @import gss
#' @export
#' @param formula Symbolic description of the model to be fit.
#' @param type List specifying the type of spline for each variable.
#' @param data Optional data frame containing the variables in the model.
#' @param alpha Parameter defining cross-validation score for smoothing parameter selection.
#' @param weights Optional vector of bin-counts for histogram data.
#' @param subset Optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action Function which indicates what should happen when the data contain NAs.
#' @param id.basis Index of observations to be used as "knots."
#' @param nbasis Number of "knots" to be used.
#' @param seed Seed to be used for the random generation of "knots."
#' @param domain Data frame specifying marginal support of density.
#' @param quad Quadrature for calculating integral. Mandatory if variables other than factors or numerical vectors are involved.
#' @param prec Precision requirement for internal iterations.
#' @param maxiter Maximum number of iterations allowed for internal iterations.
#' @param theta2 Parameters for two-way interactions.
#' @param w Optional vector to specify weights for two-way interactions.
ssden2 <- function(formula, type = NULL, data = list(), alpha = 1.4,
                   weights = NULL, subset, na.action = na.omit,
                   id.basis = NULL, nbasis = NULL, seed = NULL,
                   domain = as.list(NULL), quad = NULL,
                   prec = 1e-7, maxiter = 30, theta2 = NULL, w = NULL) {
  ## Obtain model frame and model terms
  mf <- match.call()
  mf$type <- mf$alpha <- NULL
  mf$theta2 <- NULL
  mf$w <- NULL
  mf$id.basis <- mf$nbasis <- mf$seed <- NULL
  mf$domain <- mf$quad <- mf$quad <- NULL
  mf$prec <- mf$maxiter <- NULL
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent()))
  cnt <- model.weights(mf)
  mf$"(weights)" <- NULL
  p <- dim(mf)[2]
  ## Use ssden for 1-D estimation
  if (dim(mf)[2] == 1) stop("use ssden to estimate 1-D density")
  ## Generate sub-basis
  nobs <- dim(mf)[1]
  if (is.null(id.basis)) {
    if (is.null(nbasis)) nbasis <- max(40, ceiling(12 * nobs^(2 / 9)))
    if (nbasis >= nobs) nbasis <- nobs
    if (!is.null(seed)) set.seed(seed)
    id.basis <- sample(nobs, nbasis, prob = cnt)
  } else {
    if (max(id.basis) > nobs | min(id.basis) < 1) {
      stop("gss error in ssden1: id.basis out of range")
    }
    nbasis <- length(id.basis)
  }
  ## Set domain and type, generate rho and quadrature
  if (is.null(quad)) quad <- as.list(NULL)
  rho <- rho.log <- as.list(NULL)
  rho.int <- rho.int2 <- NULL
  for (xlab in names(mf)) {
    x <- mf[[xlab]]
    if (is.factor(x)) {
      ## factor variable
      domain[[xlab]] <- NULL
      wt <- as.numeric(table(x))
      rho[[xlab]] <- wt / sum(wt)
      quad[[xlab]] <- list(pt = unique(x), wt = rho[[xlab]])
      rho.log[[xlab]] <- log(rho[[xlab]])
      rho.int <- c(rho.int, sum(rho[[xlab]] * log(rho[[xlab]])))
      rho.int2 <- c(rho.int2, sum(rho[[xlab]] * (log(rho[[xlab]]))^2))
    }
    if (is.vector(x) & !is.factor(x)) {
      ## numerical vector
      if (is.null(domain[[xlab]])) {
        mn <- min(x)
        mx <- max(x)
        domain[[xlab]] <- c(mn, mx) + c(-1, 1) * (mx - mn) * .05
      } else {
        domain[[xlab]] <- c(min(domain[[xlab]]), max(domain[[xlab]]))
      }
      if (is.null(type[[xlab]])) {
        type[[xlab]] <- list("cubic", domain[[xlab]])
      } else {
        if (length(type[[xlab]]) == 1) {
          type[[xlab]] <- list(type[[xlab]][[1]], domain[[xlab]])
        }
      }
      form <- as.formula(paste("~", xlab))
      rho[[xlab]] <- ssden(form,
        data = mf, type = type[xlab],
        domain = data.frame(domain[xlab]),
        alpha = 2, id.basis = id.basis
      )
      qd.wk <- rho[[xlab]]$quad
      rho.wk <- dssden(rho[[xlab]], qd.wk$pt)
      qd.wk$pt <- qd.wk$pt[[1]]
      qd.wk$wt <- rho.wk * qd.wk$wt
      quad[[xlab]] <- qd.wk
      rho.log[[xlab]] <- log(rho.wk)
      rho.int <- c(rho.int, sum(log(rho.wk) * qd.wk$wt))
      rho.int2 <- c(rho.int2, sum((log(rho.wk))^2 * qd.wk$wt))
    }
    if (is.matrix(x)) {
      ## numerical matrix
      if (is.null(quad[[xlab]]) | is.null(quad)) {
        stop("gss error in ssden1: no default quadrature")
      } else {
        qd.wk <- quad[[xlab]]
        qd.wk$pt <- data.frame(I(qd.wk$pt))
        colnames(qd.wk$pt) <- xlab
        form <- as.formula(paste("~", xlab))
        rho[[xlab]] <- ssden(form,
          data = mf, type = type[xlab], quad = qd.wk,
          alpha = 2, id.basis = id.basis
        )
        rho.wk <- dssden(rho[[xlab]], qd.wk$pt)
        quad[[xlab]]$wt <- rho.wk * quad[[xlab]]$wt
        rho.log[[xlab]] <- log(rho.wk)
        rho.int <- c(rho.int, sum(log(rho.wk) * quad[[xlab]]$wt))
        rho.int2 <- c(rho.int2, sum((log(rho.wk))^2 * quad[[xlab]]$wt))
      }
    }
  }
  ## Generate terms
  term <- mkterm(mf, type)
  term$labels <- term$labels[term$labels != "1"]
  int <- mkint(mf, type, id.basis, quad, term, rho.log, rho.int, p)
  ## Generate s, r, and q
  s <- r <- NULL
  nq <- 0
  for (label in term$labels) {
    x <- mf[, term[[label]]$vlist]
    x.basis <- mf[id.basis, term[[label]]$vlist]
    nphi <- term[[label]]$nphi
    nrk <- term[[label]]$nrk
    if (nphi) {
      phi <- term[[label]]$phi
      for (i in 1:nphi) s <- cbind(s, phi$fun(x, nu = i, env = phi$env))
    }
    if (nrk) {
      if (nrk == 1) {
        rk <- term[[label]]$rk
        nq <- nq + 1
        r <- array(c(r, rk$fun(x, x.basis, nu = i, env = rk$env, out = TRUE)), c(nobs, nbasis, nq))
      } else {
        rtemp <- NULL
        rk <- term[[label]]$rk
        phi <- term[[label]]$phi
        nphi <- term[[label]]$nphi
        for (i in 1:nrk) {
          rtemp <- array(c(rtemp, rk$fun(x, x.basis, nu = i, env = rk$env, out = TRUE)), c(nobs, nbasis, i))
        }
        nq <- nq + 1
        rk0 <- matrix(0, dim(rtemp)[1], dim(rtemp)[2])
        if (nphi) {
          for (j in 1:nphi) {
            phix <- phi$fun(x, j, phi$env)
            phiy <- phi$fun(x.basis, j, phi$env)
            rk0 <- rk0 + outer(phix, phiy)
          }
        }
        rtemp <- array(c(rtemp, rk0), c(nobs, nbasis, nrk + 1))
        rk1 <- matrix(0, dim(rtemp)[1], dim(rtemp)[2])
        for (i in 1:dim(rtemp)[3]) {
          rk1 <- rk1 + rtemp[, , i]
        }
        r <- array(c(r, rk1), c(nobs, nbasis, nq))
      }
    }
  }
  if (!is.null(s)) {
    s <- s[, 1:p]
  }
  q <- r[id.basis, , ]
  if (!is.null(s)) {
    nnull <- dim(s)[2]
    ## Check s rank
    if (qr(s)$rank < nnull) {
      stop("gss error in ssden1: unpenalized terms are linearly dependent")
    }
  }
  ## Use ssden for 1-D estimation
  if (nq == 1) stop("use ssden to estimate density on 1-D continuous domain")

  rbasis <- r[id.basis, , ]
  ## Fit the model
  if (is.null(w)) {
    w <- rep(1, p * (p - 1) / 2)
  } else {
    w <- w
  }
  z <- mspdsty1(s, r, id.basis, cnt, int, prec, maxiter, alpha, theta2, p, w)
  R1 <- matrix(0, dim(r)[1], dim(r)[2])
  R2 <- matrix(0, dim(r)[1], dim(r)[2])
  w <- rep(1, (p * (p - 1) / 2))
  for (i in 1:p) {
    R1 <- R1 + 10^z$theta1[i] * r[, , i]
  }
  for (i in (p + 1):(p * (p + 1) / 2)) {
    R2 <- R2 + z$theta2[i - p] * r[, , i] / w[i - p]
  }
  R <- R1 + R2
  if (!is.null(s)) {
    estimatedg <- s %*% z$d + R %*% z$c
  } else {
    estimatedg <- R %*% z$c
  }
  ginteraction <- R1 %*% z$c
  loss_int <- dot(z$int_c, z$c) + dot(z$int_d, z$d)
  loglike <- log(mean(exp(-estimatedg))) + dot(z$int_c, z$c) + dot(z$int_d, z$d)
  for (i in (p + 1):(p * (p + 1) / 2)) {
    w[i - p] <- 1 / sum((theta2[i - p] * r[, , i] %*% z$c)^2)
    if (is.nan(w[i - p]) == TRUE || is.infinite(w[i - p]) == TRUE) w[i - p] <- 1
  }
  w <- w / sqrt(sum(w^2))
  ## Brief description of model terms
  desc <- NULL
  for (label in term$labels) {
    desc <- rbind(desc, as.numeric(c(term[[label]][c("nphi", "nrk")])))
  }
  desc <- rbind(desc, apply(desc, 2, sum))
  rownames(desc) <- c(term$labels, "total")
  colnames(desc) <- c("Unpenalized", "Penalized")
  ## Return the results
  obj <- c(list(
    call = match.call(), mf = mf, cnt = cnt, terms = term, desc = desc, alpha = alpha,
    domain = domain, rho = rho, rho.int = rho.int, rho.int2 = rho.int2, quad = quad,
    id.basis = id.basis, int = int, s = s, r = r, q = q, g = estimatedg, rbasis = rbasis, ginter = ginteraction, R = R, R1 = R1, R2 = R2, loss_int = loss_int, w = w, loglike = loglike, w = w
  ), z)
  class(obj) <- c("ssden1", "ssden")
  obj
}
## Fit multiple smoothing parameter density
mspdsty1 <- function(s, r, id.basis, cnt, int, prec, maxiter, alpha, theta2, p, w) {
  nq <- dim(r)[3]
  ## initialization
  r1 <- r[, , 1:p]
  r2 <- r[, , (p + 1):(p * (p + 1) / 2)]
  theta1.0 <- -log10(apply(r1[id.basis, , ], 3, function(x) sum(diag(x))))
  theta2.0 <- -log10(apply(r2[id.basis, , ], 3, function(x) sum(diag(x))))
  theta0 <- 10^theta2.0
  theta1 <- theta1.0
  if (is.null(theta2)) {
    theta2 <- 10^theta2.0
  } else {
    theta2 <- theta2
  }

  r1.wk <- int.r1.wk <- r2.wk <- int.r2.wk <- r.wk <- 0
  for (i in 1:p) {
    r1.wk <- r1.wk + 10^theta1[i] * r[, , i]
    r.wk <- r.wk + 10^theta1[i] * r[, , i]
    int.r1.wk <- int.r1.wk + 10^theta1[i] * int$r[, i]
  }
  for (i in (p + 1):(p * (p + 1) / 2)) {
    r2.wk <- r2.wk + theta2[i - p] * r[, , i] / w[i - p]
    r.wk <- r.wk + theta2[i - p] * r[, , i] / w[i - p]
    int.r2.wk <- int.r2.wk + theta2[i - p] * int$r[, i] / w[i - p]
  }
  int1.wk <- list(r1 = int.r1.wk, s = int$s[1:p])
  int2.wk <- list(r2 = int.r2.wk, s = int$s[1:p])
  ## theta adjustment
  z <- sspdsty1(s, r.wk, r.wk[id.basis, ], r1.wk, r1.wk[id.basis, ], r2.wk, r2.wk[id.basis, ], cnt, int1.wk, int2.wk, prec, maxiter, alpha)
  theta1 <- theta1 + z$theta
  theta2 <- theta2 * 10^(z$theta)
  r1.wk <- int.r1.wk <- r2.wk <- int.r2.wk <- r.wk <- 0
  for (i in 1:p) theta1[i] <- 2 * theta1[i] + log10(t(z$c) %*% r[id.basis, , i] %*% z$c)


  for (i in 1:p) {
    r1.wk <- r1.wk + 10^theta1[i] * r[, , i]
    r.wk <- r.wk + 10^theta1[i] * r[, , i]
    int.r1.wk <- int.r1.wk + 10^theta1[i] * int$r[, i]
  }
  int1.wk <- list(r = int.r1.wk, s = int$s[1:p])
  for (i in (p + 1):(p * (p + 1) / 2)) {
    r2.wk <- r2.wk + theta2[i - p] * r[, , i] / w[i - p]
    r.wk <- r.wk + theta2[i - p] * r[, , i] / w[i - p]
    int.r2.wk <- int.r2.wk + theta2[i - p] * int$r[, i] / w[i - p]
  }
  int2.wk <- list(r2 = int.r2.wk, s = int$s[1:p])
  ## lambda search
  z <- sspdsty1(s, r.wk, r.wk[id.basis, ], r1.wk, r1.wk[id.basis, ], r2.wk, r2.wk[id.basis, ], cnt, int1.wk, int2.wk, prec, maxiter, alpha)
  lambda <- z$lambda
  theta1 <- theta1 + z$theta
  theta2 <- theta2 * 10^(z$theta)
  cd <- c(z$c, z$d)
  scal <- NULL
  ## return
  z$theta1 <- theta1
  z$theta2 <- theta2
  z$theta0 <- theta0
  return(z)
}

## Fit single smoothing parameter density
sspdsty1 <- function(s, r, q, r1, q1, r2, q2, cnt, int1, int2, prec, maxiter, alpha) {
  nxi <- dim(r)[2]
  nobs <- dim(r)[1]
  if (!is.null(s)) {
    nnull <- dim(s)[2]
  } else {
    nnull <- 0
  }
  nxis <- nxi + nnull
  if (is.null(cnt)) cnt <- 0
  if (sum(cnt)) {
    wt <- cnt / sum(cnt)
  } else {
    wt <- 1 / nobs
  }
  if (is.null(s)) {
    int1$s <- NULL
  }
  ## cv function
  cv <- function(lambda) {
    fit <- .Fortran("dnewton10",
      cd = as.double(cd), as.integer(nxis),
      as.double(10^(lambda + theta) * q), as.integer(nxi),
      as.double(cbind(10^theta * r, s)), as.integer(nobs),
      as.integer(sum(cnt)), as.integer(cnt),
      as.double(c(10^theta * int1$r + 10^theta * int2$r, int1$s)),
      as.double(prec), as.integer(maxiter),
      as.double(.Machine$double.eps), integer(nxis),
      wk = double(2 * nobs + nxis * (nxis + 3)),
      info = integer(1), PACKAGE = "gss"
    )
    if (fit$info == 1) stop("gss error in ssden1: Newton iteration diverges")
    if (fit$info == 2) warning("gss warning in ssden1: Newton iteration fails to converge")
    aa <- fit$wk[1:nobs]
    assign("cd", fit$cd, inherits = TRUE)
    eta0 <- cbind(10^theta * r, s) %*% cd
    wwt <- wt * exp(-eta0)
    wwt <- wwt / sum(wwt)
    assign("scal", sum(wt * exp(-eta0)), inherits = TRUE)
    trc <- sum(wwt * exp(aa / (1 - aa))) - 1
    cv <- sum(c(10^theta * int1$r + 10^theta * int2$r, int1$s) * cd) + log(scal) + alpha * trc
    alpha.wk <- max(0, log.la0 - lambda - 5) * (3 - alpha) + alpha
    alpha.wk <- min(alpha.wk, 3)
    adj <- ifelse(alpha.wk > alpha, (alpha.wk - alpha) * trc, 0)
    cv + adj
  }
  ## initialization
  mu.r <- apply(wt * r, 2, sum)
  v.r <- apply(wt * r^2, 2, sum)
  if (!is.null(s)) {
    mu.s <- apply(wt * s, 2, sum)
    v.s <- apply(wt * s^2, 2, sum)
  }
  if (is.null(s)) {
    theta <- 0
  } else {
    theta <- log10(sum(v.s - mu.s^2) / nnull / sum(v.r - mu.r^2) * nxi) / 2
  }
  log.la0 <- log10(sum(v.r - mu.r^2) / sum(diag(q))) + theta
  ## lambda search
  cd <- rep(0, nxi + nnull)
  la <- log.la0
  tol <- 0
  scal <- NULL
  mn0 <- log.la0 - 6
  mx0 <- log.la0 + 6
  repeat {
    mn <- max(la - 1, mn0)
    mx <- min(la + 1, mx0)
    zz <- nlm0(cv, c(mn, mx))
    if ((min(zz$est - mn, mx - zz$est) >= 1e-1) ||
      (min(zz$est - mn0, mx0 - zz$est) < 1e-1)) {
      break
    } else {
      la <- zz$est
    }
  }
  ## return
  jk1 <- cv(zz$est)
  c <- cd[1:nxi]
  if (nnull) {
    d <- cd[nxi + (1:nnull)]
  } else {
    d <- NULL
  }
  int_c <- 10^theta * int1$r + 10^theta * int2$r
  int_d <- int1$s
  list(lambda = zz$est, theta = theta, c = c, d = d, scal = scal, cv = zz$min, int_c = int_c, int_d = int_d)
}

## Calculate integrals of phi and rk for ssden1
mkint <- function(mf, type, id.basis, quad, term, rho, rho.int, p) {
  ## Obtain model terms
  mt <- attr(mf, "terms")
  xvars <- as.character(attr(mt, "variables"))[-1]
  xfacs <- attr(mt, "factors")
  term.labels <- labels(mt)
  vlist <- xvars[as.logical(apply(xfacs, 1, sum))]
  ## Set types for marginals
  var.type <- NULL
  for (xlab in vlist) {
    x <- mf[, xlab]
    if (!is.null(type[[xlab]])) {
      ## Check consistency and set default parameters
      type.wk <- type[[xlab]][[1]]
      if
      (!(type.wk %in% c(
          "ordinal", "nominal", "cubic", "linear", "per",
          "cubic.per", "linear.per", "tp", "sphere", "custom"
        ))) {
        stop("gss error in mkint: unknown type")
      }
      if (type.wk %in% c("ordinal", "nominal")) {
        par.wk <- NULL
        if (!is.factor(x)) {
          stop("gss error in mkint: wrong type")
        }
      }
      if (type.wk %in% c("cubic", "linear")) {
        if (length(type[[xlab]]) == 1) {
          mn <- min(x)
          mx <- max(x)
          par.wk <- c(mn, mx) + c(-1, 1) * .05 * (mx - mn)
        } else {
          par.wk <- type[[xlab]][[2]]
        }
        if (is.factor(x) | !is.vector(x)) {
          stop("gss error in mkint: wrong type")
        }
      }
      if (type.wk %in% c("per", "cubic.per", "linear.per")) {
        if (type.wk == "per") type.wk <- "cubic.per"
        if (length(type[[xlab]]) == 1) {
          stop("gss error in mkint: missing domain of periodicity")
        } else {
          par.wk <- type[[xlab]][[2]]
        }
        if (is.factor(x) | !is.vector(x)) {
          stop("gss error in mkint: wrong type")
        }
      }
      if (type.wk == "tp") {
        if (length(type[[xlab]]) == 1) {
          par.wk <- list(order = 2, mesh = x, weight = 1)
        } else {
          par.wk <- par.wk1 <- type[[xlab]][[2]]
          if (length(par.wk1) == 1) {
            par.wk <- list(order = par.wk1, mesh = x, weight = 1)
          }
          if (is.null(par.wk$mesh)) par.wk$mesh <- x
          if (is.null(par.wk$weight)) par.wk$weight <- 1
        }
        if (dim(as.matrix(x))[2] != dim(as.matrix(par.wk$mesh))[2]) {
          stop("gss error in mkint: wrong dimension in normalizing mesh")
        }
      }
      if (type.wk == "sphere") {
        if (length(type[[xlab]]) == 1) {
          par.wk <- 2
        } else {
          par.wk <- type[[xlab]][[2]]
        }
        if (!(par.wk %in% (2:4))) {
          stop("gss error in mkint: spherical order not implemented")
        }
      }
      if (type.wk == "custom") par.wk <- type[[xlab]][[2]]
    } else {
      ## Set default types
      if (is.factor(x)) {
        ## categorical variable
        if (is.ordered(x)) {
          type.wk <- "ordinal"
        } else {
          type.wk <- "nominal"
        }
        par.wk <- NULL
      } else {
        ## numerical variable
        if (is.vector(x)) {
          type.wk <- "cubic"
          mn <- min(x)
          mx <- max(x)
          par.wk <- c(mn, mx) + c(-1, 1) * .05 * (mx - mn)
        } else {
          type.wk <- "tp"
          par.wk <- list(order = 2, mesh = x, weight = 1)
        }
      }
    }
    var.type[[xlab]] <- list(type.wk, par.wk)
  }
  ## Create phi and rk
  nbasis <- length(id.basis)
  nvar <- length(names(mf))
  s <- r <- s.rho <- r.rho <- NULL
  ns <- nq <- 0
  for (label in term.labels) {
    ns <- ns + term[[label]]$nphi
    nq <- nq + term[[label]]$nrk
    vlist <- xvars[as.logical(xfacs[, label])]
    x <- mf[, vlist]
    dm <- length(vlist)
    phi <- rk <- NULL
    if (dm == 1) {
      type.wk <- var.type[[vlist]][[1]]
      xx <- mf[id.basis, vlist]
      xmesh <- quad[[vlist]]$pt
      if (type.wk %in% c("nominal", "ordinal")) {
        ## factor variable
        if (type.wk == "nominal") {
          fun <- mkrk.nominal(levels(x))
        } else {
          fun <- mkrk.ordinal(levels(x))
        }
        if (nlevels(x) > 2) {
          ## rk
          rk <- fun$fun(xmesh, xx, fun$env, TRUE)
        } else {
          ## phi
          wk <- as.factor(names(fun$env$code)[1])
          phi <- fun$fun(xmesh, wk, fun$env)
        }
      }
      if (type.wk == "cubic") {
        ## cubic splines
        range <- var.type[[vlist]][[2]]
        ## phi
        phi.fun <- mkphi.cubic(range)
        phi <- phi.fun$fun(xmesh, 1, phi.fun$env)
        ## rk
        rk.fun <- mkrk.cubic(range)
        rk <- rk.fun$fun(xmesh, xx, rk.fun$env, TRUE)
      }
      if (type.wk %in% c("cubic.per", "linear", "linear.per", "sphere")) {
        ## cubic periodic, linear, and linear periodic splines
        range <- var.type[[vlist]][[2]]
        ## rk
        if (type.wk == "cubic.per") rk.fun <- mkrk.cubic.per(range)
        if (type.wk == "linear") rk.fun <- mkrk.linear(range)
        if (type.wk == "linear.per") rk.fun <- mkrk.linear.per(range)
        if (type.wk == "sphere") rk.fun <- mkrk.sphere(range)
        rk <- rk.fun$fun(xmesh, xx, rk.fun$env, TRUE)
      }
      if (type.wk == "tp") {
        ## thin-plate splines
        par <- var.type[[vlist]][[2]]
        order <- par$order
        mesh <- par$mesh
        weight <- par$weight
        if (is.vector(x)) {
          xdim <- 1
        } else {
          xdim <- dim(x)[2]
        }
        ## phi
        phi.fun <- mkphi.tp(xdim, order, mesh, weight)
        nphi <- choose(xdim + order - 1, xdim) - 1
        if (nphi > 0) {
          for (nu in 1:nphi) {
            phi <- cbind(phi, phi.fun$fun(xmesh, nu, phi.fun$env))
          }
        }
        ## rk
        rk.fun <- mkrk.tp(xdim, order, mesh, weight)
        rk <- rk.fun$fun(xmesh, xx, rk.fun$env, TRUE)
      }
      if (type.wk == "custom") {
        ## user-defined
        par <- var.type[[vlist]][[2]]
        nphi <- par$nphi
        if (nphi > 0) {
          phi.fun <- par$mkphi(par$env)
          for (nu in 1:nphi) {
            phi <- cbind(phi, phi.fun$fun(xmesh, nu, phi.fun$env))
          }
        }
        rk.fun <- par$mkrk(par$env)
        rk <- rk.fun$fun(xmesh, xx, rk.fun$env, TRUE)
      }
      wmesh <- quad[[vlist]]$wt
      if (!is.null(phi)) {
        s.rho.wk <- rho.int * sum(wmesh * phi)
        s.rho.wk[names(mf) == vlist] <- sum(wmesh * rho[[vlist]] * phi)
        s <- c(s, sum(wmesh * phi))
        s.rho <- c(s.rho, sum(s.rho.wk))
      }
      if (!is.null(rk)) {
        r.rho.wk <- outer(apply(wmesh * rk, 2, sum), rho.int)
        r.rho.wk[, names(mf) == vlist] <- apply(wmesh * rho[[vlist]] * rk, 2, sum)
        r <- cbind(r, apply(wmesh * rk, 2, sum))
        r.rho <- cbind(r.rho, apply(r.rho.wk, 1, sum))
      }
    } else {
      bin.fac <- n.phi <- phi.list <- rk.list <- NULL
      for (i in 1:dm) {
        type.wk <- var.type[[vlist[i]]][[1]]
        if (type.wk %in% c("nominal", "ordinal")) {
          ## factor variable
          if (type.wk == "nominal") {
            rk.wk <- mkrk.nominal(levels(x[[i]]))
          } else {
            rk.wk <- mkrk.ordinal(levels(x[[i]]))
          }
          phi.wk <- rk.wk
          n.phi <- c(n.phi, 0)
          bin.fac <- c(bin.fac, !(nlevels(x[[i]]) > 2))
        }
        if (type.wk == "cubic") {
          ## cubic or linear splines
          range <- var.type[[vlist[i]]][[2]]
          ## phi
          phi.wk <- mkphi.cubic(range)
          n.phi <- c(n.phi, 1)
          ## rk
          rk.wk <- mkrk.cubic(range)
          bin.fac <- c(bin.fac, 0)
        }
        if (type.wk %in% c("cubic.per", "linear", "linear.per", "sphere")) {
          ## cubic periodic, linear, or linear periodic splines
          range <- var.type[[vlist[i]]][[2]]
          n.phi <- c(n.phi, 0)
          phi.wk <- NULL
          # phi.wk<-1
          if (type.wk == "cubic.per") rk.wk <- mkrk.cubic.per(range)
          if (type.wk == "linear") rk.wk <- mkrk.linear(range)
          if (type.wk == "linear.per") rk.wk <- mkrk.linear.per(range)
          if (type.wk == "sphere") rk.wk <- mkrk.sphere(range)
          bin.fac <- c(bin.fac, 0)
        }
        if (type.wk == "tp") {
          ## thin-plate splines
          par <- var.type[[vlist[i]]][[2]]
          order <- par$order
          mesh <- par$mesh
          weight <- par$weight
          if (is.vector(x[[i]])) {
            xdim <- 1
          } else {
            xdim <- dim(x[[i]])[2]
          }
          phi.wk <- mkphi.tp(xdim, order, mesh, weight)
          n.phi <- c(n.phi, choose(xdim + order - 1, xdim) - 1)
          rk.wk <- mkrk.tp(xdim, order, mesh, weight)
          bin.fac <- c(bin.fac, 0)
        }
        if (type.wk == "custom") {
          ## user-defined
          par <- var.type[[vlist[i]]][[2]]
          n.phi <- c(n.phi, par$nphi)
          if (par$nphi > 0) {
            phi.wk <- par$mkphi(par$env)
          } else {
            phi.wk <- NULL
          }
          rk.wk <- par$mkrk(par$env)
          bin.fac <- c(bin.fac, 0)
        }
        phi.list <- c(phi.list, list(phi.wk))
        rk.list <- c(rk.list, list(rk.wk))
      }
      ## phi
      id0 <- names(mf) %in% vlist
      nphi <- term[[label]]$nphi
      iphi <- term[[label]]$iphi
      if (nphi > 0) {
        for (nu in 1:nphi) {
          ind <- nu - 1
          s.wk <- 1
          s.rho.wk <- rho.int
          s.rho.wk[id0] <- 1
          for (i in 1:dm) {
            phi.wk <- phi.list[[i]]
            xmesh <- quad[[vlist[i]]]$pt
            if (bin.fac[i]) {
              wk <- as.factor(names(phi.wk$env$code)[1])
              phi <- phi.wk$fun(xmesh, wk, phi.wk$env)
            } else {
              code <- ind %% n.phi[i] + 1
              ind <- ind %/% n.phi[i]
              phi <- phi.wk$fun(xmesh, code, phi.wk$env)
            }
            wmesh <- quad[[vlist[i]]]$wt
            s.wk <- s.wk * sum(wmesh * phi)
            id1 <- names(mf) == vlist[i]
            s.rho.wk[id1] <- s.rho.wk[id1] * sum(wmesh * rho[[vlist[i]]] * phi)
            s.rho.wk[!id1] <- s.rho.wk[!id1] * sum(wmesh * phi)
          }
          s <- c(s, s.wk)
          s.rho <- c(s.rho, sum(s.rho.wk))
        }
      }
      ## rk
      rtemp <- rtemp.rho <- NULL
      n.rk <- ifelse(n.phi, 2, 1)
      nrk <- prod(n.rk) - as.logical(nphi)
      if (nrk > 0) {
        for (nu in 1:nrk) {
          ind <- nu - !nphi
          r.wk <- 1
          r.rho.wk <- outer(rep(1, nbasis), rho.int)
          r.rho.wk[, id0] <- 1
          for (i in 1:dm) {
            code <- ind %% n.rk[i] + 1
            ind <- ind %/% n.rk[i]
            xx <- mf[id.basis, vlist[[i]]]
            xmesh <- quad[[vlist[i]]]$pt
            if (code == n.rk[i]) {
              rk.wk <- rk.list[[i]]
              rk <- rk.wk$fun(xmesh, xx, rk.wk$env, TRUE)
            } else {
              rk <- 0
              phi.wk <- phi.list[[i]]
              for (j in 1:n.phi[i]) {
                phix <- phi.wk$fun(xmesh, j, phi.wk$env)
                phiy <- phi.wk$fun(xx, j, phi.wk$env)
                rk <- rk + outer(phix, phiy)
              }
            }
            wmesh <- quad[[vlist[i]]]$wt
            r.wk <- r.wk * apply(wmesh * rk, 2, sum)
            id1 <- names(mf) == vlist[i]
            r.rho.wk[, id1] <- r.rho.wk[, id1] * apply(wmesh * rho[[vlist[i]]] * rk, 2, sum)
            r.rho.wk[, !id1] <- r.rho.wk[, !id1] * apply(wmesh * rk, 2, sum)
          }
          rtemp <- cbind(rtemp, r.wk)
          rtemp.rho <- cbind(rtemp.rho, apply(r.rho.wk, 1, sum))
        }
        if (nrk > 1 & nphi > 0) {
          r.wk <- 1
          r.rho.wk <- outer(rep(1, nbasis), rho.int)
          r.rho.wk[, id0] <- 1
          for (i in 1:dm) {
            rk <- 0
            xx <- mf[id.basis, vlist[[i]]]
            xmesh <- quad[[vlist[i]]]$pt
            phi.wk <- phi.list[[i]]
            for (j in 1:n.phi[i]) {
              phix <- phi.wk$fun(xmesh, j, phi.wk$env)
              phiy <- phi.wk$fun(xx, j, phi.wk$env)
              rk <- rk + outer(phix, phiy)
            }
            wmesh <- quad[[vlist[i]]]$wt
            r.wk <- r.wk * apply(wmesh * rk, 2, sum)
            id1 <- names(mf) == vlist[i]
            r.rho.wk[, id1] <- r.rho.wk[, id1] * apply(wmesh * rho[[vlist[i]]] * rk, 2, sum)
            r.rho.wk[, !id1] <- r.rho.wk[, !id1] * apply(wmesh * rk, 2, sum)
          }
          rtemp <- cbind(rtemp, r.wk)
          rtemp.rho <- cbind(rtemp.rho, apply(r.rho.wk, 1, sum))
        }
        r <- cbind(r, apply(rtemp, 1, sum))
        r.rho <- cbind(r.rho, apply(rtemp.rho, 1, sum))
      }
    }
  }
  list(s = s[1:p], r = r, s.rho = s.rho[1:p], r.rho = r.rho, var.type = var.type)
}
mkrk.linear <- function(range) {
  ## Create the environment
  env <- list(min = min(range), max = max(range))
  ## Create the rk function
  fun <- function(x, y, env, outer.prod = FALSE) {
    ## % Check the inputs
    if (!(is.vector(x) & is.vector(y))) {
      stop("gss error in rk: inputs are of wrong types")
    }
    if ((min(x, y) < env$min) | (max(x, y) > env$max)) {
      stop("gss error in rk: inputs are out of range")
    }
    ## % Scale the inputs
    x <- (x - env$min) / (env$max - env$min)
    y <- (y - env$min) / (env$max - env$min)
    ## % Return the result
    rk <- function(x, y) {
      k1 <- function(x) (x - .5)
      k2 <- function(x) ((x - .5)^2 - 1 / 12) / 2
      k1(x) * k1(y) + k2(abs(x - y))
    }
    if (outer.prod) {
      outer(x, y, rk)
    } else {
      rk(x, y)
    }
  }
  ## Return the function and the environment
  list(fun = fun, env = env)
}
