#' Fit conditional density model
#' @import gss
#' @export
#' @param formula Symbolic description of the model to be fit.
#' @param response Formula listing response variables.
#' @param type List specifying the type of spline for each variable.
#' @param data Optional data frame containing the variables in the model.
#' @param weights Optional vector of bin-counts for histogram data.
#' @param subset Optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action Function which indicates what should happen when the data contain NAs.
#' @param alpha Parameter defining cross-validation score for smoothing parameter selection.
#' @param id.basis Index of observations to be used as "knots."
#' @param nbasis Number of "knots" to be used.
#' @param seed Seed to be used for the random generation of "knots."
#' @param rho Method to construct rho function.
#' @param ydomain Data frame specifying marginal support of conditional density.
#' @param yquad Quadrature for calculating integral on Y domain. Mandatory if response variables other than factors or numerical vectors are involved.
#' @param prec Precision requirement for internal iterations.
#' @param maxiter Maximum number of iterations allowed for internal iterations.
#' @param skip.iter Flag indicating whether to use initial values of theta and skip theta iteration.
#' @param p Dimension of data frame.
#' @param theta2 Parameters for two-way interactions.
#' @param w2 Weights for two-way interactions.
sscden_selection <- function(formula, response, type = NULL, data = list(), weights,
                             subset, na.action = na.omit, alpha = 1.4,
                             id.basis = NULL, nbasis = NULL, seed = NULL, rho = list("xy"),
                             ydomain = as.list(NULL), yquad = NULL,
                             prec = 1e-7, maxiter = 30, skip.iter = TRUE, p = 2, theta2 = NULL, w2 = NULL) {
  ## Obtain model frame and model terms
  mf <- match.call()
  mf$response <- mf$type <- mf$alpha <- NULL
  mf$id.basis <- mf$nbasis <- mf$seed <- mf$rho <- NULL
  mf$ydomain <- mf$yquad <- mf$theta2 <- mf$w2 <- NULL
  mf$prec <- mf$maxiter <- mf$skip.iter <- mf$p <- NULL
  term.wk <- terms.formula(formula)
  ynames <- as.character(attr(terms(response), "variables"))[-1]
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  nobs <- nrow(mf)
  cnt <- model.weights(mf)
  if (is.null(cnt)) {
    data$cnt <- rep(1, nobs)
  } else {
    data$cnt <- cnt
    mf$"(weights)" <- NULL
  }
  ## Generate sub-basis
  if (is.null(id.basis)) {
    if (is.null(nbasis)) nbasis <- max(30, ceiling(10 * nobs^(2 / 9)))
    if (nbasis >= nobs) nbasis <- nobs
    if (!is.null(seed)) set.seed(seed)
    id.basis <- sample(nobs, nbasis, prob = cnt)
  } else {
    if (max(id.basis) > nobs | min(id.basis) < 1) {
      stop("gss error in sscden1: id.basis out of range")
    }
    nbasis <- length(id.basis)
  }
  ## Check inputs
  mt <- attr(mf, "terms")
  vars <- as.character(attr(mt, "variables"))[-1]
  if (!all(ynames %in% vars)) stop("gss error in sscden1: response missing in model")
  xnames <- vars[!(vars %in% ynames)]
  if (is.null(xnames)) stop("gss error in sscden1: missing covariate")
  ## Set type for given ydomain
  fac.list <- NULL
  for (ylab in ynames) {
    y <- mf[[ylab]]
    if (is.factor(y)) {
      fac.list <- c(fac.list, ylab)
      ydomain[[ylab]] <- NULL
    } else {
      if (!is.vector(y) & is.null(yquad)) {
        stop("gss error in sscden1: no default quadrature")
      }
      if (is.vector(y)) {
        if (is.null(ydomain[[ylab]])) {
          mn <- min(y)
          mx <- max(y)
          ydomain[[ylab]] <- c(mn, mx) + c(-1, 1) * (mx - mn) * .05
        } else {
          ydomain[[ylab]] <- c(min(ydomain[[ylab]]), max(ydomain[[ylab]]))
        }
        if (is.null(type[[ylab]])) {
          type[[ylab]] <- list("cubic", ydomain[[ylab]])
        } else {
          if (length(type[[ylab]]) == 1) {
            type[[ylab]] <- list(type[[ylab]][[1]], ydomain[[ylab]])
          }
        }
      }
    }
  }
  ## Generate terms
  term <- mkterm(mf, type)
  term$labels <- term$labels[term$labels != "1"]
  ## obtain unique covariate observations
  x <- xx <- mf[, xnames, drop = FALSE]
  xx <- apply(xx, 1, function(x) paste(x, collapse = "\r"))
  x.dup.ind <- duplicated(xx)
  if (!is.null(cnt)) xx <- rep(xx, cnt)
  xx.wt <- as.vector(table(xx)[unique(xx)])
  xx.wt <- xx.wt / sum(xx.wt)
  nx <- length(xx.wt)
  ## calculate rho
  if (is.null(rho$fun)) {
    type <- rho[[1]]
    if (type == "y") {
      yfac <- TRUE
      for (ylab in ynames) yfac <- yfac & is.factor(mf[, ylab])
      if (!yfac) {
        if (is.null(cnt)) cntt <- rep(1, dim(mf)[1])
        rho <- ssden0(response,
          data = data, weights = cnt, id.basis = id.basis,
          alpha = 2, domain = ydomain, quad = yquad
        )
        qd.pt <- rho$quad$pt
        qd.wt <- rho$quad$wt
        env <- list(ydomain = ydomain, qd.pt = qd.pt, qd.wt = qd.wt, rho = rho)
        fun <- function(x, y, env, outer.prod = FALSE) {
          if (!outer.prod) {
            dssden(env$rho, y)
          } else {
            t(matrix(dssden(env$rho, y), dim(y)[1], dim(x)[1]))
          }
        }
      } else {
        qd.pt <- data.frame(levels(mf[, ynames[1]]), stringsAsFactors = TRUE)
        if (length(ynames) > 1) {
          for (ylab in ynames[-1]) {
            wk <- expand.grid(levels(mf[, ylab]), 1:dim(qd.pt)[1])
            qd.pt <- data.frame(qd.pt[wk[, 2], ], wk[, 1], stringsAsFactors = TRUE)
          }
        }
        colnames(qd.pt) <- ynames
        qd.wt <- as.vector(table(mf[, rev(ynames)]))
        qd.wt <- qd.wt / sum(qd.wt)
        env <- list(qd.pt = qd.pt, qd.wt = qd.wt)
        fun <- function(x, y, env, outer.prod = FALSE) {
          if (!outer.prod) {
            rep(1, dim(x)[1])
          } else {
            matrix(1, dim(x)[1], dim(y)[1])
          }
        }
      }
      rho <- list(fun = fun, env = env)
    }
    if (type == "xy") {
      ydomain <- data.frame(ydomain)
      mn <- ydomain[1, ]
      mx <- ydomain[2, ]
      dm <- ncol(ydomain)
      if (dm == 1) {
        ## Gauss-Legendre quadrature
        quad <- gauss.quad(200, c(mn, mx))
        quad$pt <- data.frame(quad$pt)
        colnames(quad$pt) <- colnames(ydomain)
      } else {
        ## Smolyak cubature
        qdsz.depth <- switch(min(dm, 6) - 1,
          18,
          14,
          10,
          9,
          7
        )
        quad <- smolyak.quad(dm, qdsz.depth)
        for (i in 1:ncol(ydomain)) {
          ylab <- colnames(ydomain)[i]
          wk <- mf[[ylab]]
          jk <- ssden(~wk,
            domain = data.frame(wk = ydomain[, i]), alpha = 2,
            id.basis = id.basis, weights = cnt
          )
          quad$pt[, i] <- qssden(jk, quad$pt[, i])
          quad$wt <- quad$wt / dssden(jk, quad$pt[, i])
        }
        jk <- wk <- NULL
        quad$pt <- data.frame(quad$pt)
        colnames(quad$pt) <- colnames(ydomain)
      }
      ## Incorporate factors in quadrature
      if (!is.null(fac.list)) {
        for (i in 1:length(fac.list)) {
          wk <- expand.grid(levels(mf[[fac.list[i]]]), 1:length(quad$wt))
          quad$wt <- quad$wt[wk[, 2]]
          col.names <- c(fac.list[i], colnames(quad$pt))
          quad$pt <- data.frame(wk[, 1], quad$pt[wk[, 2], ], stringsAsFactors = TRUE)
          colnames(quad$pt) <- col.names
        }
      }
      rho <- list(NULL)
      for (ylab in ynames) {
        if (is.numeric(mf[[ylab]])) {
          form <- as.formula(paste(ylab, "~", paste(xnames, collapse = "+")))
          rho[[ylab]] <- ssanova(form, data = mf, id.basis = id.basis)
        }
        if (is.factor(mf[[ylab]])) {
          form <- as.formula(paste("~(", paste(xnames, collapse = "+"), ")*", ylab))
          resp <- as.formula(paste("~", ylab))
          rho[[ylab]] <- ssllrm(form, resp, data = mf, id.basis = id.basis)
        }
      }
      env <- list(ynames = ynames, ydomain = ydomain, qd.pt = quad$pt, qd.wt = quad$wt, rho = rho)
      fun <- function(x, y, env, outer.prod = FALSE) {
        z <- 1
        for (ylab in env$ynames) {
          yy <- y[[ylab]]
          if (is.numeric(yy)) {
            mu <- predict(env$rho[[ylab]], x)
            sigma <- sqrt(env$rho[[ylab]]$varht)
            ymn <- env$ydomain[1, ylab]
            ymx <- env$ydomain[2, ylab]
            if (!outer.prod) {
              wk <- dnorm((yy - mu) / sigma) /
                (pnorm((ymx - mu) / sigma) - pnorm((ymn - mu) / sigma))
              z <- z * wk
            } else {
              wk <- t(outer(yy, mu, dnorm, sigma)) /
                (pnorm((ymx - mu) / sigma) - pnorm((ymn - mu) / sigma))
              z <- z * wk
            }
          }
          if (is.factor(yy)) {
            wk <- predict(env$rho[[ylab]], x)
            if (!outer.prod) {
              wk1 <- NULL
              for (i in 1:length(yy)) {
                wk1 <- c(wk1, wk[i, yy[i] == env$rho[[ylab]]$qd.pt])
              }
              z <- z * wk1
            } else {
              wk1 <- NULL
              for (i in 1:length(yy)) {
                wk1 <- cbind(wk1, wk[, yy[i] == env$rho[[ylab]]$qd.pt])
              }
              z <- z * wk1
            }
          }
        }
        z
      }
      rho <- list(fun = fun, env = env)
    }
  }
  ## Generate s, r, int.s, and int.r
  rho.wk <- rho$fun(x[!x.dup.ind, , drop = FALSE], rho$env$qd.pt, rho$env, outer = TRUE)
  rho.wk <- t(t(rho.wk) * rho$env$qd.wt)
  rho.wk1 <- apply(rho.wk * xx.wt, 2, sum)
  nmesh <- length(rho$env$qd.wt)
  s <- r <- int.s <- int.r <- NULL
  r0 <- NULL
  id.s <- id.r <- NULL
  id.s.list <- id.r.list <- list(NULL)
  nu <- nq <- 0
  nq0 <- 0
  for (label in term$labels) {
    vlist <- term[[label]]$vlist
    x.list <- xnames[xnames %in% vlist]
    y.list <- ynames[ynames %in% vlist]
    xy <- mf[, vlist]
    xy.basis <- mf[id.basis, vlist]
    qd.xy <- data.frame(matrix(0, nmesh, length(vlist)))
    names(qd.xy) <- vlist
    qd.xy[, y.list] <- rho$env$qd.pt[, y.list]
    if (length(x.list)) {
      xx <- x[!x.dup.ind, x.list, drop = FALSE]
    } else {
      xx <- NULL
    }
    nphi <- term[[label]]$nphi
    nrk <- term[[label]]$nrk
    if (nphi) {
      phi <- term[[label]]$phi
      id.s.list[[label]] <- NULL
      for (i in 1:nphi) {
        nu <- nu + 1
        s.wk <- phi$fun(xy, nu = i, env = phi$env)
        s <- cbind(s, s.wk)
        if (is.null(xx)) {
          id.s <- c(id.s, nu)
          id.s.list[[label]] <- c(id.s.list[[label]], nu)
          qd.s.wk <- phi$fun(qd.xy[, , drop = TRUE], nu = i, env = phi$env)
          int.s <- c(int.s, sum(qd.s.wk * rho.wk1))
        } else {
          if (length(y.list) == 0) {
            names(xx) <- x.list
            int.s <- c(int.s, sum(phi$fun(xx[, , drop = TRUE], i, phi$env) * xx.wt))
          } else {
            id.s <- c(id.s, nu)
            id.s.list[[label]] <- c(id.s.list[[label]], nu)
            int.s.wk <- 0
            for (j in 1:nx) {
              qd.xy[, x.list] <- xx[rep(j, nmesh), ]
              qd.s.wk <- phi$fun(qd.xy, i, phi$env)
              int.s.wk <- int.s.wk + sum(qd.s.wk * rho.wk[j, ]) * xx.wt[j]
            }
            int.s <- c(int.s, int.s.wk)
          }
        }
      }
    }
    if (nrk) {
      if (nrk == 1) {
        rk <- term[[label]]$rk
        id.r.list[[label]] <- NULL
        for (i in 1:nrk) {
          nq <- nq + 1
          nq0 <- nq0 + 1
          r.wk <- rk$fun(xy, xy.basis, nu = i, env = rk$env, out = TRUE)
          r <- array(c(r, r.wk), c(nobs, nbasis, nq))
          r0 <- array(c(r0, r.wk), c(nobs, nbasis, nq0))
          if (is.null(xx)) {
            id.r <- c(id.r, nq)
            id.r.list[[label]] <- c(id.r.list[[label]], nq)
            qd.r.wk <- rk$fun(qd.xy[, , drop = TRUE], xy.basis, nu = i, env = rk$env, out = TRUE)
            int.r <- cbind(int.r, apply(rho.wk1 * qd.r.wk, 2, sum))
          } else {
            if (length(y.list) == 0) {
              names(xx) <- x.list
              qd.r.wk <- rk$fun(xx[, , drop = TRUE], xy.basis, i, rk$env, TRUE)
              int.r <- cbind(int.r, apply(xx.wt * qd.r.wk, 2, sum))
            } else {
              id.r <- c(id.r, nq)
              id.r.list[[label]] <- c(id.r.list[[label]], nq)
              int.r.wk <- 0
              for (j in 1:nx) {
                qd.xy[, x.list] <- xx[rep(j, nmesh), ]
                qd.r.wk <- rk$fun(qd.xy, xy.basis, i, rk$env, TRUE)
                int.r.wk <- int.r.wk + apply(rho.wk[j, ] * qd.r.wk, 2, sum) * xx.wt[j]
              }
              int.r <- cbind(int.r, int.r.wk)
            }
          }
        }
      }
      if (nrk > 1) {
        rtemp <- NULL
        int.r_temp <- NULL
        rk <- term[[label]]$rk
        phi <- term[[label]]$phi
        id.r.list[[label]] <- NULL
        for (i in 1:nrk) {
          nq0 <- nq0 + 1
          rtemp.wk <- rk$fun(xy, xy.basis, nu = i, env = rk$env, out = TRUE)
          rtemp <- array(c(rtemp, rtemp.wk), c(nobs, nbasis, i))
          r0 <- array(c(r0, rtemp.wk), c(nobs, nbasis, nq0))
          if (is.null(xx)) {
            qd.r.wk <- rk$fun(qd.xy[, , drop = TRUE], xy.basis, nu = i, env = rk$env, out = TRUE)
            int.r_temp <- cbind(int.r_temp, apply(rho.wk1 * qd.r.wk, 2, sum))
          } else {
            if (length(y.list) == 0) {
              names(xx) <- x.list
              qd.r.wk <- rk$fun(xx[, , drop = TRUE], xy.basis, i, rk$env, TRUE)
              int.r_temp <- cbind(int.r_temp, apply(xx.wt * qd.r.wk, 2, sum))
            } else {
              int.r.wk <- 0
              for (j in 1:nx) {
                qd.xy[, x.list] <- xx[rep(j, nmesh), ]
                qd.r.wk <- rk$fun(qd.xy, xy.basis, i, rk$env, TRUE)
                int.r.wk <- int.r.wk + apply(rho.wk[j, ] * qd.r.wk, 2, sum) * xx.wt[j]
              }
              int.r_temp <- cbind(int.r_temp, int.r.wk)
            }
          }
        }
        rk0 <- matrix(0, dim(rtemp)[1], dim(rtemp)[2])
        for (j in 1:nphi) {
          phix <- phi$fun(xy, j, phi$env)
          phiy <- phi$fun(xy.basis, j, phi$env)
          rk0 <- rk0 + outer(phix, phiy)
        }
        nq0 <- nq0 + 1
        rtemp <- array(c(rtemp, rk0), c(nobs, nbasis, nrk + 1))
        r0 <- array(c(r0, rk0), c(nobs, nbasis, nq0))
        rk1 <- matrix(0, dim(rtemp)[1], dim(rtemp)[2])
        for (i in 1:dim(rtemp)[3]) {
          rk1 <- rk1 + rtemp[, , i]
        }
        nq <- nq + 1
        r <- array(c(r, rk1), c(nobs, nbasis, nq))
        id.r <- c(id.r, nq)
        id.r.list[[label]] <- c(id.r.list[[label]], nq)
        if (is.null(xx)) {
          qd.r.wk_temp <- 0
          for (j in 1:nphi) {
            phix <- phi$fun(qd.xy[, , drop = TRUE], j, phi$env)
            phiy <- phi$fun(xy.basis, j, phi$env)
            qd.r.wk_temp <- qd.r.wk_temp + outer(phix, phiy)
          }
          int.r_temp <- cbind(int.r_temp, apply(rho.wk1 * qd.r.wk_temp, 2, sum))
        } else {
          if (length(y.list) == 0) {
            names(xx) <- x.list
            qd.r.wk_temp <- 0
            for (j in 1:nphi) {
              phix <- phi$fun(xx[, , drop = TRUE], j, phi$env)
              phiy <- phi$fun(xy.basis, j, phi$env)
              qd.r.wk_temp <- qd.r.wk_temp + outer(phix, phiy)
            }
            int.r_temp <- cbind(int.r_temp, apply(xx.wt * qd.r.wk_temp, 2, sum))
          } else {
            int.r.wk <- 0
            for (j in 1:nx) {
              qd.xy[, x.list] <- xx[rep(j, nmesh), ]
              qd.r.wk_temp <- 0
              for (k in 1:nphi) {
                phix <- phi$fun(qd.xy, k, phi$env)
                phiy <- phi$fun(xy.basis, k, phi$env)
                qd.r.wk_temp <- qd.r.wk_temp + outer(phix, phiy)
              }
              int.r.wk <- int.r.wk + apply(rho.wk[j, ] * qd.r.wk_temp, 2, sum) * xx.wt[j]
            }
            int.r_temp <- cbind(int.r_temp, int.r.wk)
          }
        }

        int.r <- cbind(int.r, apply(int.r_temp, 1, sum))
      }
    }
  }
  if (!is.null(s)) {
    s <- s[, 1:p]
    int.s <- int.s[1:p]
  }
  ## Check s rank
  if (!is.null(s)) {
    nnull <- dim(s)[2]
    if (qr(s)$rank < nnull) {
      stop("gss error in sscden1: unpenalized MLE is not unique")
    }
  }
  if (is.null(w2)) {
    w2 <- rep(1, p - 1)
  } else {
    w2 <- w2
  }
  ## Fit the model
  z <- mspcdsty1(s, r, id.basis, cnt, int.s, int.r, prec, maxiter, alpha, skip.iter, theta2, w2, p)
  R1 <- matrix(0, dim(r)[1], dim(r)[2])
  R2 <- matrix(0, dim(r)[1], dim(r)[2])
  for (i in 1:p) {
    R1 <- R1 + 10^z$theta1[i] * r[, , i]
  }
  for (i in (p + 1):(2 * p - 1)) {
    R2 <- R2 + z$theta2[i - p] * r[, , i] / w2[i - p]
  }
  R <- R1 + R2
  if (!is.null(s)) {
    estimatedg <- s %*% z$d + R %*% z$c
  } else {
    estimatedg <- R %*% z$c
  }
  loss_int <- dot(z$int_c, z$c) + dot(int.s, z$d)
  rbasis <- r[id.basis, , ]
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
    call = match.call(), mf = mf, cnt = cnt, terms = term, desc = desc, rho = rho,
    alpha = alpha, ynames = ynames, xnames = xnames,
    x.dup.ind = x.dup.ind, xx.wt = xx.wt, id.s.list = id.s.list, id.r.list = id.r.list,
    id.s = id.s, id.r = id.r, id.basis = id.basis, skip.iter = skip.iter, s = s, r = r, r0 = r0, R = R, R1 = R1, R2 = R2, rbasis = rbasis, g = estimatedg, int.s = int.s, int.r = int.r, loss_int = loss_int, w2 = w2
  ), z)
  class(obj) <- c("sscden1", "sscden")
  obj
}

## Fit (multiple smoothing parameter) conditional density model
mspcdsty1 <- function(s, r, id.basis, cnt, int.s, int.r, prec, maxiter, alpha, skip.iter, theta2, w2, p) {
  nobs <- dim(r)[1]
  nxi <- dim(r)[2]
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
  ## cv functions
  cv.s <- function(lambda) {
    fit <- .Fortran("cdennewton10",
      cd = as.double(cd), as.integer(nxis),
      as.double(10^lambda * q.wk), as.integer(nxi),
      as.double(cbind(r.wk, s)), as.integer(nobs),
      as.integer(sum(cnt)), as.integer(cnt),
      as.double(c(int.r.wk, int.s)),
      as.double(prec), as.integer(maxiter),
      as.double(.Machine$double.eps), integer(nxis),
      wk = double(2 * nobs + nxis * (nxis + 3)),
      info = integer(1), PACKAGE = "gss"
    )
    if (fit$info == 1) stop("gss error in sscden1: Newton iteration diverges")
    if (fit$info == 2) warning("gss warning in sscden1: Newton iteration fails to converge")
    aa <- fit$wk[1:nobs]
    assign("cd", fit$cd, inherits = TRUE)
    eta0 <- cbind(r.wk, s) %*% cd
    wwt <- wt * exp(-eta0)
    wwt <- wwt / sum(wwt)
    assign("scal", sum(wt * exp(-eta0)), inherits = TRUE)
    trc <- sum(wwt * exp(aa / (1 - aa))) - 1
    cv <- sum(c(int.r.wk, int.s) * cd) + log(scal) + alpha * trc
    alpha.wk <- max(0, log.la0 - lambda - 5) * (3 - alpha) + alpha
    alpha.wk <- min(alpha.wk, 3)
    adj <- ifelse(alpha.wk > alpha, (alpha.wk - alpha) * trc, 0)
    cv + adj
  }
  cv.m <- function(theta1) {
    ind.wk <- theta1 != theta1.old
    if (sum(ind.wk) == p) {
      r1.wk0 <- int.r1.wk0 <- 0
      for (i in 1:p) {
        r1.wk0 <- r1.wk0 + 10^theta1[i] * r[, , i]
        int.r1.wk0 <- int.r1.wk0 + 10^theta1[i] * int.r[, i]
      }
      assign("r1.wk", r1.wk0 + 0, inherits = TRUE)
      assign("int.r1.wk", int.r1.wk0 + 0, inherits = TRUE)
      assign("theta1.old", theta1 + 0, inherits = TRUE)
    } else {
      r1.wk0 <- r1.wk
      int.r1.wk0 <- int.r1.wk
      for (i in (1:p)[ind.wk]) {
        theta1.wk <- (10^(theta1[i] - theta1.old[i]) - 1) * 10^theta1.old[i]
        r1.wk0 <- r1.wk0 + theta1.wk * r[, , i]
        int.r1.wk0 <- int.r1.wk0 + theta1.wk * int.r[, i]
      }
    }
    r.wk0 <- r1.wk0
    int.r.wk0 <- int.r1.wk0
    for (i in (p + 1):(2 * p - 1)) {
      r.wk0 <- r.wk0 + theta2[i - p] * r[, , i] / w2[i - p]
      int.r.wk0 <- int.r.wk0 + theta2[i - p] * int.r[, i] / w2[i - p]
    }
    q.wk <- r.wk0[id.basis, ]
    fit <- .Fortran("cdennewton10",
      cd = as.double(cd), as.integer(nxis),
      as.double(10^lambda * q.wk), as.integer(nxi),
      as.double(cbind(r.wk0, s)), as.integer(nobs),
      as.integer(sum(cnt)), as.integer(cnt),
      as.double(c(int.r.wk0, int.s)),
      as.double(prec), as.integer(maxiter),
      as.double(.Machine$double.eps), integer(nxis),
      wk = double(2 * nobs + nxis * (nxis + 3)),
      info = integer(1), PACKAGE = "gss"
    )
    if (fit$info == 1) stop("gss error in sscden1: Newton iteration diverges")
    if (fit$info == 2) warning("gss warning in sscden1: Newton iteration fails to converge")
    aa <- fit$wk[1:nobs]
    assign("cd", fit$cd, inherits = TRUE)
    eta0 <- cbind(r.wk, s) %*% cd
    wwt <- wt * exp(-eta0)
    wwt <- wwt / sum(wwt)
    assign("scal", sum(wt * exp(-eta0)), inherits = TRUE)
    trc <- sum(wwt * exp(aa / (1 - aa))) - 1
    cv <- sum(c(int.r.wk0, int.s) * cd) + log(scal) + alpha * trc
    alpha.wk <- max(0, theta1 - log.th0 - 5) * (3 - alpha) + alpha
    alpha.wk <- min(alpha.wk, 3)
    adj <- ifelse(alpha.wk > alpha, (alpha.wk - alpha) * trc, 0)
    cv + adj
  }
  cv.wk <- function(theta1) cv.scale * cv.m(theta1) + cv.shift
  ## Initialization
  theta1.0 <- -log10(apply(r[id.basis, , 1:p, drop = FALSE], 3, function(x) sum(diag(x))))
  theta2.0 <- -log10(apply(r[id.basis, , (p + 1):(2 * p - 1), drop = FALSE], 3, function(x) sum(diag(x))))
  theta2_isnull <- FALSE
  theta1 <- theta1.0
  if (is.null(theta2)) {
    theta2 <- 10^theta2.0
    theta2_isnull <- TRUE
  } else {
    theta2 <- theta2
  }
  r1.wk <- int.r1.wk <- 0
  for (i in 1:p) {
    r1.wk <- r1.wk + 10^theta1[i] * r[, , i]
    int.r1.wk <- int.r1.wk + 10^theta1[i] * int.r[, i]
  }
  r.wk <- r1.wk
  int.r.wk <- int.r1.wk
  for (i in (p + 1):(2 * p - 1)) {
    r.wk <- r.wk + theta2[i - p] * r[, , i]
    int.r.wk <- int.r.wk + theta2[i - p] * int.r[, i] / w2[i - p]
  }
  if (!nnull) {
    mu.r <- apply(wt * r.wk, 2, sum)
    v.r <- apply(wt * r.wk^2, 2, sum)
    v.r <- v.r - mu.r^2
    theta.wk <- 0
  } else {
    mu.r <- apply(wt * r.wk, 2, sum)
    v.r <- apply(wt * r.wk^2, 2, sum)
    v.r <- v.r - mu.r^2
    mu.s <- apply(wt * s, 2, sum)
    v.s <- apply(wt * s^2, 2, sum)
    v.s <- v.s - mu.s^2
    theta.wk <- log10(sum(v.s) / nnull / sum(v.r) * nxi) / 2
  }
  theta1 <- theta1 + theta.wk
  theta2 <- theta2 * 10^theta.wk
  r.wk <- 10^theta.wk * r.wk
  int.r.wk <- 10^theta.wk * int.r.wk
  q.wk <- r.wk[id.basis, ]
  log.la0 <- log10(sum(v.r) / sum(diag(q.wk))) + 2 * theta.wk
  ## fixed theta iteration
  cd <- rep(0, nxi + nnull)
  scal <- NULL
  la <- log.la0
  mn0 <- log.la0 - 6
  mx0 <- log.la0 + 6
  repeat {
    mn <- max(la - 1, mn0)
    mx <- min(la + 1, mx0)
    zz <- nlm0(cv.s, c(mn, mx))
    if ((min(zz$est - mn, mx - zz$est) >= 1e-1) ||
      (min(zz$est - mn0, mx0 - zz$est) < 1e-1)) {
      break
    } else {
      la <- zz$est
    }
  }
  if (p == 1) {
    lambda <- zz$est
    c <- cd[1:nxi]
    if (nnull) {
      d <- cd[nxi + (1:nnull)]
    } else {
      d <- NULL
    }
    return(list(lambda = lambda, theta1 = theta1, theta2 = theta2, c = c, d = d, cv = zz$min, scal = scal, int_c = int.r.wk))
  }
  ## theta adjustment
  r1.wk <- int.r1.wk <- 0
  for (i in 1:p) {
    theta1[i] <- 2 * theta1[i] + log10((t(cd[1:nxi])) %*% r[id.basis, , i] %*% (cd[1:nxi]))
    r1.wk <- r1.wk + 10^theta1[i] * r[, , i]
    int.r1.wk <- int.r1.wk + 10^theta1[i] * int.r[, i]
  }
  r.wk <- r1.wk
  int.r.wk <- int.r1.wk
  for (i in (p + 1):(2 * p - 1)) {
    if (theta2_isnull) {
      theta2[i - p] <- theta2[i - p]^2 * (t(cd[1:nxi]) %*% r[id.basis, , i] %*% (cd[1:nxi]))
    }
    r.wk <- r.wk + theta2[i - p] * r[, , i] / w2[i - p]
    int.r.wk <- int.r.wk + theta2[i - p] * int.r[, i] / w2[i - p]
  }
  if (!nnull) {
    mu.r <- apply(wt * r.wk, 2, sum)
    v.r <- apply(wt * r.wk^2, 2, sum)
    v.r <- v.r - mu.r^2
    theta.wk <- 0
  } else {
    mu.r <- apply(wt * r.wk, 2, sum)
    v.r <- apply(wt * r.wk^2, 2, sum)
    v.r <- v.r - mu.r^2
    mu.s <- apply(wt * s, 2, sum)
    v.s <- apply(wt * s^2, 2, sum)
    v.s <- v.s - mu.s^2
    theta.wk <- log10(sum(v.s) / nnull / sum(v.r) * nxi) / 2
  }
  theta1 <- theta1 + theta.wk
  theta2 <- theta2 * 10^theta.wk
  r.wk <- 10^theta.wk * r.wk
  int.r.wk <- 10^theta.wk * int.r.wk
  q.wk <- r.wk[id.basis, ]
  log.la0 <- log10(sum(v.r) / sum(diag(q.wk))) + 2 * theta.wk
  log.th0 <- theta1 - log.la0
  ## fixed theta iteration
  cd <- rep(0, nxi + nnull)
  la <- log.la0
  mn0 <- log.la0 - 6
  mx0 <- log.la0 + 6
  repeat {
    mn <- max(la - 1, mn0)
    mx <- min(la + 1, mx0)
    zz <- nlm0(cv.s, c(mn, mx))
    if ((min(zz$est - mn, mx - zz$est) >= 1e-1) ||
      (min(zz$est - mn0, mx0 - zz$est) < 1e-1)) {
      break
    } else {
      la <- zz$est
    }
  }
  lambda <- zz$est
  ## early return
  if (skip.iter) {
    c <- cd[1:nxi]
    if (nnull) {
      d <- cd[nxi + (1:nnull)]
    } else {
      d <- NULL
    }
    return(list(lambda = lambda, theta1 = theta1, theta2 = theta2, c = c, d = d, cv = zz$min, scal = scal, int_c = int.r.wk))
  }
  ## theta search
  counter <- 0
  r1.wk <- int.r1.wk <- 0
  for (i in 1:p) {
    r1.wk <- r1.wk + 10^theta1[i] * r[, , i]
    int.r1.wk <- int.r1.wk + 10^theta1[i] * int.r[, i]
  }
  r.wk <- r1.wk
  int.r.wk <- int.r1.wk
  for (i in (p + 1):(2 * p - 1)) {
    r.wk <- r.wk + theta2[i - p] * r[, , i] / w2[i - p]
    int.r.wk <- int.r.wk + theta2[i - p] * int.r[, i] / w2[i - p]
  }
  theta1.old <- theta1
  tmp <- abs(cv.m(theta1))
  cv.scale <- 1
  cv.shift <- 0
  if (tmp < 1 & tmp > 10^(-4)) {
    cv.scale <- 10 / tmp
    cv.shift <- 0
  }
  if (tmp < 10^(-4)) {
    cv.scale <- 10^2
    cv.shift <- 10
  }
  repeat {
    zz <- nlm(cv.wk, theta1, stepmax = 1, ndigit = 7)
    if (zz$code <= 3) break
    theta1 <- zz$est
    counter <- counter + 1
    if (counter >= 5) {
      warning("gss warning in sscden1: CV iteration fails to converge")
      break
    }
  }
  ## return
  c <- cd[1:nxi]
  if (nnull) {
    d <- cd[nxi + (1:nnull)]
  } else {
    d <- NULL
  }
  cv <- (zz$min - cv.shift) / cv.scale
  list(lambda = lambda, theta1 = zz$est, theta2 = theta2, c = c, d = d, cv = cv, scal = scal, int_c = int.r.wk)
}
