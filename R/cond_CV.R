#' Neighborhood selection method with cross-validation to select the tuning parameter.
#' @description Apply 5-fold cross-validation to select tuning parameter.
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
cond_select <- function(data, formula, response, type = NULL, alpha, subset, na.action, rho, ydomain, yquad, prec, maxiter, skip.iter,
                        M0, M_list, maxiteration, tolerance, id.basis = NULL, theta2, w2 = NULL, n, p) {
  loop <- 1
  ltheta2 <- length(theta2)
  Theta2 <- matrix(0, ltheta2, maxiteration + 1)
  Theta2[, loop] <- theta2
  l <- rep(0, maxiteration + 1)
  l[loop] <- 0
  diff2 <- 1
  d <- rep(0, maxiteration + 1)
  d[1] <- diff2
  zerodiff <- 1

  while (diff2 >= tolerance & loop < maxiteration & (zerodiff) >= 1) {
    loop <- loop + 1
    densityfit <- sscden_selection(
      formula = formula, response = response, data = data, w2 = w2, type = type, alpha = alpha, subset = subset,
      na.action = na.action, rho = rho, ydomain = ydomain, yquad = yquad, prec = prec, maxiter = maxiter, skip.iter = skip.iter, id.basis = id.basis, p = p, theta2 = theta2
    )
    M <- tune.M(densityfit, w2, M0, M_list, k = 5, n, formula, response, data, p)[[1]]
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
    Gmatrix <- (-kronecker(diag(ltheta2), t(c1)) %*% t(U2) %*% t2 + t(B2) %*% c1 + la1 * kronecker(diag(ltheta2), t(c1)) %*% t(u2) %*% c1) / w2
    Hpart1 <- U2 %*% kronecker(diag(ltheta2), c1)
    Hpart2 <- t(t2) %*% U2 %*% kronecker(diag(ltheta2), c1)
    Hpart3 <- sqrt(t2) * Hpart1
    Hmatrix <- (t(Hpart3) %*% Hpart3 - t(Hpart2) %*% Hpart2)
    Hmatrix <- diag(1 / w2) %*% Hmatrix %*% diag(1 / w2)

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
    theta2[theta2 < 1e-8] <- 0
    Theta2[, loop] <- theta2
    if (loop <= 3) {
      zerodiff <- 1
    } else {
      zerodiff <- (sum(theta2 > 0) - sum(Theta2[, loop - 3] > 0))
    }
    diff2 <- sqrt(sum((Theta2[, loop] - Theta2[, loop - 1])^2)) / (sqrt(sum((Theta2[, loop - 1])^2)) + 1e-6)
    d[loop] <- diff2
  }

  thre_list <- c(1e-2, 1e-3, 1e-4)
  count <- 0
  thre_mat <- rep(0, length(thre_list))
  int.r <- densityfit$int.r
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
    int.r.wk <- 0
    for (i in 1:p) {
      int.r.wk <- int.r.wk + 10^theta1[i] * int.r[, i]
    }
    for (i in (p + 1):(2 * p - 1)) {
      int.r.wk <- int.r.wk + theta20[i - p] * int.r[, i] / w2[i - p]
    }
    thre_mat[count] <- log(score / n) + (dot(int.r.wk, c1) + dot(densityfit$int.s, d1)) / n
  }
  min_pos <- which.min(thre_mat)
  bound <- thre_list[min_pos]
  theta2[theta2 < bound] <- 0
  list(loop = loop, theta2 = theta2, densityfit = densityfit)
}


tune.M <- function(densityfit, w2, M0, M_list, k = 5, n, formula, response, data, p) {
  M_list1 <- M0 * M_list
  cv_mat <- rep(0, length(M_list1))
  count <- 0
  id.basis <- densityfit$id.basis
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
  la1 <- 10^densityfit$lambda / 2
  n_k <- n / k
  n_k_omit <- n - n_k
  for (M in M_list1) {
    count <- count + 1
    loglikFold <- matrix(NA, ncol = k, nrow = 1)
    intFold <- matrix(NA, ncol = k, nrow = 1)
    grps <- cut(1:n, k, labels = FALSE)
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
      fit_int <- sscden_int(formula = formula, response = response, data = data, omit = omit, id.basis = id.basis, seed = 5732, p = p)
      B2 <- fit_int$int.r[, (p + 1):(ltheta2 + p)] * n / n_k_omit
      la1 <- 10^densityfit$lambda / 2
      Gmatrix <- (-kronecker(diag(ltheta2), t(c1)) %*% t(U2) %*% t2 + t(B2) %*% c1 + la1 * kronecker(diag(ltheta2), t(c1)) %*% t(u2) %*% c1) / w2
      Hpart1 <- U2 %*% kronecker(diag(ltheta2), c1)
      Hpart2 <- t(t2) %*% U2 %*% kronecker(diag(ltheta2), c1)
      Hpart3 <- sqrt(t2) * Hpart1
      Hmatrix <- (t(Hpart3) %*% Hpart3 - t(Hpart2) %*% Hpart2)
      Hmatrix <- diag(1 / w2) %*% Hmatrix %*% diag(1 / w2)

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
      theta2[theta2 < 1e-8] <- 0
      fit_int_omit <- sscden_int_omit(formula = formula, response = response, data = data, omit = omit, id.basis = id.basis, seed = 5732, p = p)
      int.r_k <- fit_int_omit$int.r
      int.r.wk <- 0
      for (i in 1:p) {
        int.r.wk <- int.r.wk + 10^theta1.1[i] * int.r_k[, i]
      }
      for (i in (p + 1):(2 * p - 1)) {
        int.r.wk <- int.r.wk + theta2[i - p] * int.r_k[, i] / w2[i - p]
      }
      r <- densityfit$r
      R <- densityfit$R1
      for (i in (p + 1):(2 * p - 1)) {
        R <- R + theta2[i - p] * r[, , i] / w2[i - p]
      }
      cv_g <- densityfit$s %*% d1 + R %*% c1
      cv_g <- cv_g[omit]
      score <- (sum(exp(-cv_g)))
      loglikFold[, kth] <- score
      intFold[, kth] <- dot(int.r.wk, c1)
    }

    cv_mat[count] <- log((apply(loglikFold, 1, sum)) / n) + (apply(intFold, 1, sum))
  }


  min_error_pos <- which.min(cv_mat)
  M_opt <- M_list1[min_error_pos]
  list(M_opt = M_opt)
}
sscden_int <- function(formula, response, type = NULL, data = list(), omit, weights,
                       subset, na.action = na.omit, alpha = 1.4,
                       id.basis = NULL, nbasis = NULL, seed = NULL, rho = list("xy"),
                       ydomain = as.list(NULL), yquad = NULL,
                       prec = 1e-7, maxiter = 30, skip.iter = TRUE, p = 2, theta2 = NULL, w2 = NULL) {
  ## Obtain model frame and model terms
  mf <- match.call()
  mf$response <- mf$type <- mf$alpha <- NULL
  mf$id.basis <- mf$nbasis <- mf$seed <- mf$rho <- NULL
  mf$ydomain <- mf$yquad <- mf$theta2 <- mf$w2 <- mf$omit <- NULL
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
  v <- seq(1, nx, length.out = nx)
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
  id.s <- id.r <- NULL
  id.s.list <- id.r.list <- list(NULL)
  nu <- nq <- 0
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
            for (j in v[-omit]) {
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
          r.wk <- rk$fun(xy, xy.basis, nu = i, env = rk$env, out = TRUE)
          r <- array(c(r, r.wk), c(nobs, nbasis, nq))
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
              for (j in v[-omit]) {
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
          rtemp.wk <- rk$fun(xy, xy.basis, nu = i, env = rk$env, out = TRUE)
          rtemp <- array(c(rtemp, rtemp.wk), c(nobs, nbasis, i))
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
              for (j in v[-omit]) {
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
        rtemp <- array(c(rtemp, rk0), c(nobs, nbasis, nrk + 1))
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
            for (j in v[-omit]) {
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
    id.s = id.s, id.r = id.r, id.basis = id.basis, skip.iter = skip.iter, s = s, r = r, int.s = int.s, int.r = int.r
  ))
  class(obj) <- c("sscden1", "sscden")
  obj
}
sscden_int_omit <- function(formula, response, type = NULL, data = list(), omit, weights,
                            subset, na.action = na.omit, alpha = 1.4,
                            id.basis = NULL, nbasis = NULL, seed = NULL, rho = list("xy"),
                            ydomain = as.list(NULL), yquad = NULL,
                            prec = 1e-7, maxiter = 30, skip.iter = TRUE, p = 2, theta2 = NULL, w2 = NULL) {
  ## Obtain model frame and model terms
  mf <- match.call()
  mf$response <- mf$type <- mf$alpha <- NULL
  mf$id.basis <- mf$nbasis <- mf$seed <- mf$rho <- NULL
  mf$ydomain <- mf$yquad <- mf$theta2 <- mf$w2 <- mf$omit <- NULL
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
  v <- seq(1, nx, length.out = nx)
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
  id.s <- id.r <- NULL
  id.s.list <- id.r.list <- list(NULL)
  nu <- nq <- 0
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
            for (j in v[omit]) {
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
          r.wk <- rk$fun(xy, xy.basis, nu = i, env = rk$env, out = TRUE)
          r <- array(c(r, r.wk), c(nobs, nbasis, nq))
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
              for (j in v[omit]) {
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
          rtemp.wk <- rk$fun(xy, xy.basis, nu = i, env = rk$env, out = TRUE)
          rtemp <- array(c(rtemp, rtemp.wk), c(nobs, nbasis, i))
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
              for (j in v[omit]) {
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
        rtemp <- array(c(rtemp, rk0), c(nobs, nbasis, nrk + 1))
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
            for (j in v[omit]) {
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
    id.s = id.s, id.r = id.r, id.basis = id.basis, skip.iter = skip.iter, s = s, r = r, int.s = int.s, int.r = int.r
  ))
  class(obj) <- c("sscden1", "sscden")
  obj
}
