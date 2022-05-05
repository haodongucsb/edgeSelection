#' Semi-parametric method for edge selection.
#' @export
#' @import gss
#' @param data Data Frame.
#' @param ... Any options can be defined.
#' \itemize{
#' \item \code{px} Dimension of variables in the semi-parametric method to be estimated using non-parametric method.
#' \item \code{maxLambda} Number used to generate the range of the tuning parameter for selection of Lambda matrix.
#' \item \code{maxTheta} Number used to generate the range of the tuning parameter for selection of Theta matrix.
#' \item \code{iterLambda} Number of iterations to find optimal Lambda matrix.
#' \item \code{iterTheta} Number of iterations to find optimal Theta matrix.
#' \item \code{cutoff} Cutoff value for squared projection.
#' \item \code{N} Number of simulations in non-parametric part to calculate the standard error.
#' \item \code{semiMethod} Method type to select optimal matrix in semi-parametric method.
#' }
#' @usage
#' selection.semi(data, ...)
#' @examples
#' library(huge)
#' n <- 200; p <- 20; px <- 5
#' # Simulate high dimension data
#' set.seed(5732)
#' z <- huge.generator(n, d = p, graph = "random", prob = .2, verbose = FALSE, vis = FALSE, v = .65)
#' data <- data.frame(z$data)
#' # Assume we estimate the first five variables using non-parametric method.
#' edge.selection(data = data, family ="semi", px = px)
selection.semi <- function(data, ...) {
  n <- dim(data)[1]
  p <- dim(data)[2]
  px <- 0
  maxLambda <- 0.5
  maxTheta <- 0.5
  iterLambda <- 5
  iterTheta <- 5
  cutoff <- 0.03
  N <- 1e3
  semiMethod <- "quic"
  params <- list(...)
  for (name in names(params)) {
    assign(name, params[[name]])
  }
  if (px == 0) {
    stop("specify indexes of non-Gaussian variables")
  }
  x <- data[, 1:px]
  y <- data[, (px + 1):p]
  paraOut <- selectBICForAlg1(y, x, maxLambda, maxTheta, iterLambda, iterTheta, 1e-3,
    useEBIC = FALSE, method = semiMethod
  )
  Theta <- paraOut$opt_Theta
  Lambda <- paraOut$opt_Lam
  optlam_cov <- paraOut$optlam_cov
  optlam_F <- paraOut$optlam_F
  Delta <- (-1 / 2) * t(Theta) %*% Mpower(Lambda, -1) %*% Theta
  mn <- apply(x, 2, min)
  mx <- apply(x, 2, max)
  numterm <- px + px * (px - 1) / 2
  domain <- data.frame(rbind(mn, mx))
  # semi parametric estimate
  semi.ssres <- ssX(x, Delta, semiflag = TRUE, cutoff = cutoff, N = N)
  semi.CxxEst <- XandX(semi.ssres, diag(x = .5, px))
  semi.CxxEst <- semi.CxxEst + t(semi.CxxEst)
  edgeMatrix <- rbind(cbind(Lambda, Theta), cbind(t(Theta), semi.CxxEst))
  edgeMatrix
}
Mpower <- function(Sig, q) {
  a <- svd(Sig)
  d <- a$d^(q)
  m <- a$u %*% diag(d) %*% t(a$v)
  m
}
cdF_nw <- function(Sigma, Cx, Cyx, Fini, lam, thr_F) {
  diff_F <- 2 * thr_F + 1
  eps <- 1e-7
  p <- dim(Fini)[1]
  px <- dim(Fini)[2]
  F <- matrix(0, nrow = p, ncol = px)
  Fini0 <- Fini
  for (i in 1:p)
  {
    for (j in 1:px)
    {
      ac_set <- Cyx[i, j] + (Cx %*% t(Fini) %*% Sigma)[j, i]
      if (abs(ac_set) <= lam / 2 && Fini[i, j] == 0) next
      gij <- -ac_set + Cx[j, j] * Sigma[i, i] * Fini[i, j]
      if (abs(gij) > lam * 0.5) {
        preabs <- abs(gij) - lam * 0.5
        coef <- max(eps, Cx[j, j]) * max(eps, Sigma[i, i])
        absF <- preabs / coef
        F[i, j] <- sign(gij) * absF
      } else {
        F[i, j] <- 0
      }
      Fini[i, j] <- F[i, j]
    }
  }
  list(F = F)
}
Alg1 <- function(y, x, lam_cov, lam_F, EnoughIter, thr) {
  eps <- 1e-4
  prevfx <- 1e15
  n <- dim(y)[1]
  p <- dim(y)[2]
  px <- dim(x)[2]
  smnum <- 1e-7
  Cy <- 1 / n * t(y) %*% y
  Cyx <- 1 / n * t(y) %*% x
  Cx <- 1 / n * t(x) %*% x
  diff_Theta <- diff_Lambda <- 2 * thr + 1
  Ipx <- diag(px)
  S <- Cy - Cyx %*% solve(Cx) %*% t(Cyx)
  # initial values
  if (abs(det(S)) < smnum) {
    Lambda <- solve(S + smnum * diag(p))
  } else {
    Lambda <- solve(S) # theta
  }
  Lambda_init <- Lambda
  # from MLE
  Theta_init <- -Lambda %*% Cyx %*% solve(Cx) # theta_xy
  # code start
  # Lambda_init <- diag(p)
  # Theta_init <- matrix(0,p,px)
  # tC%*% Cyx %*% solve(Cx)  true theta_xy
  iter <- 0
  Sigma_old <- S
  Theta_old <- Theta_init
  Lambda_old <- Lambda_init
  sXmatrix <- Cy %*% Lambda_old + Sigma_old %*% Theta_old %*% Cx %*% t(Theta_old)
  trsX <- sum(diag(sXmatrix))
  l1normLambda <- lam_cov * (sum(Lambda_old != 0) - p)
  l1normTheta <- lam_F * sum(Theta_old != 0)
  fx <- -log(det(Lambda_old)) + trsX + l1normLambda + l1normTheta
  record1 <- record2 <- rep(0, EnoughIter)
  # while( (iter<EnoughIter) && ((diff_Theta > thr) ||
  #                              (diff_Lambda > thr))
  #        && (abs((fx - prevfx)/fx) > eps))
  while ((iter < EnoughIter) && (abs((fx - prevfx) / fx) > eps))
  # while( (iter<EnoughIter) && ((diff_Theta > thr) ||
  #                                (diff_Lambda > thr)))
  {
    iter <- iter + 1
    # print("Updating Lambda...")
    dat <- y + x %*% t(Theta_old) %*% Sigma_old
    S_temp <- 1 / n * t(dat) %*% dat
    # S_temp <- cov2cor(S_temp)
    obj_quic <- QUIC(S_temp,
      rho = lam_cov,
      path = NULL, msg = 0,
      X.init = Lambda_old, W.init = Sigma_old,
      maxIter = 1
    )
    Sigma_new <- obj_quic$W
    Lambda_new <- obj_quic$X
    obF1 <- cdF_nw(Sigma_new, Cx, Cyx, Theta_old, lam_F, 1e-3)
    Theta_new <- obF1$F
    diff_Theta <- mean(abs(Theta_new - Theta_old))
    diff_Lambda <- mean(abs(Lambda_new - Lambda_old))
    record1[iter] <- mean(abs(Theta_new - Theta_old))
    record2[iter] <- mean(abs(Lambda_new - Lambda_old))
    sXmatrix_new <- Cy %*% Lambda_new + Sigma_new %*% Theta_new %*% Cx %*% t(Theta_new)
    trsX <- sum(diag(sXmatrix_new))
    l1normLambda <- lam_cov * (sum(Lambda_new != 0) - p)
    l1normTheta <- lam_F * sum(Theta_new != 0)
    fx1 <- -log(det(Lambda_new)) + trsX + l1normLambda + l1normTheta
    prevfx <- fx
    fx <- fx1
    Theta_old <- Theta_new
    Lambda_old <- Lambda_new
    Sigma_old <- Sigma_new
    # cat("Stopped at iteration: ", iter, "\n")
  }
  list(Theta = Theta_new, Lambda = Lambda_new)
}

Alg1_alasso <- function(y, x, lam_cov, lam_F, EnoughIter, thr) {
  eps <- 2e-2
  prevfx <- 1e15
  n <- dim(y)[1]
  p <- dim(y)[2]
  px <- dim(x)[2]
  smnum <- 1e-7
  Cy <- 1 / n * t(y) %*% y
  Cyx <- 1 / n * t(y) %*% x
  Cx <- 1 / n * t(x) %*% x
  gamma <- .5
  diff_Theta <- diff_Lambda <- 2 * thr + 1
  Ipx <- diag(px)
  S <- Cy - Cyx %*% solve(Cx) %*% t(Cyx)
  # initial values
  if (abs(det(S)) < smnum) {
    Lambda <- solve(S + smnum * diag(p))
  } else {
    Lambda <- solve(S) # theta
  }
  Lambda_init <- Lambda
  # from MLE
  Theta_init <- -Lambda %*% Cyx %*% solve(Cx) # theta_xy
  # tC%*% Cyx %*% solve(Cx)  true theta_xy
  iter <- 0
  Sigma_old <- S
  Theta_old <- Theta_init
  Lambda_old <- Lambda_init
  sXmatrix <- Cy %*% Lambda_old + Sigma_old %*% Theta_old %*% Cx %*% t(Theta_old)
  trsX <- sum(diag(sXmatrix))
  l1normLambda <- lam_cov * (sum(Lambda_old != 0) - p)
  l1normTheta <- lam_F * sum(Theta_old != 0)
  fx <- -log(det(Lambda_old)) + trsX + l1normLambda + l1normTheta
  record1 <- record2 <- rep(0, EnoughIter)
  while ((iter < EnoughIter) && (abs((fx - prevfx) / fx) > eps)) {
    iter <- iter + 1
    # cat("iteration at: ", iter, "\n")
    # cat("change in objective function: ", abs((fx - prevfx)/fx), "\n")
    # print("Updating Lambda...")
    dat <- y + x %*% t(Theta_old) %*% Sigma_old
    S_temp <- 1 / n * t(dat) %*% dat
    fit.lasso <- glasso(S_temp, lam_cov, maxit = 1, penalize.diagonal = FALSE)
    wi.lasso <- fit.lasso$wi
    rhomat <- lam_cov / p / 2 * matrix(1, p, p) / (pmax(abs(wi.lasso)^gamma, 1e-05))
    fit.alasso <- glasso(S_temp, rhomat, maxit = 1, penalize.diagonal = FALSE)
    Lambda_new <- fit.alasso$wi
    Sigma_new <- fit.alasso$w
    obF1 <- cdF_nw(Sigma_new, Cx, Cyx, Theta_old, lam_F, 1e-3)
    Theta_new <- obF1$F
    sXmatrix_new <- Cy %*% Lambda_new + Sigma_new %*% Theta_new %*% Cx %*% t(Theta_new)
    trsX <- sum(diag(sXmatrix_new))
    l1normLambda <- lam_cov * (sum(Lambda_new != 0) - p)
    l1normTheta <- lam_F * sum(Theta_new != 0)
    fx1 <- -log(det(Lambda_new)) + trsX + l1normLambda + l1normTheta
    # cat("fx1 value:", fx1, "\n")
    prevfx <- fx
    fx <- fx1
    if (is.na(abs((fx - prevfx) / fx))) break
    Theta_old <- Theta_new
    Lambda_old <- Lambda_new
    Sigma_old <- Sigma_new
  }
  list(Theta = Theta_new, Lambda = Lambda_new)
}
scadrightderv <- function(lamhat, a, lam) {
  pmax(lam * ((lamhat <= lam) + pmax(a * lam - lamhat, 0) * (lamhat > lam) / (a - 1)), 1e-10)
}
Alg1_scad <- function(y, x, lam_cov, lam_F, EnoughIter, thr) {
  eps <- 1e-4
  prevfx <- 1e15
  n <- dim(y)[1]
  p <- dim(y)[2]
  px <- dim(x)[2]
  smnum <- 1e-7
  Cy <- 1 / n * t(y) %*% y
  Cyx <- 1 / n * t(y) %*% x
  Cx <- 1 / n * t(x) %*% x
  gamma <- .5
  diff_Theta <- diff_Lambda <- 2 * thr + 1
  Ipx <- diag(px)
  S <- Cy - Cyx %*% solve(Cx) %*% t(Cyx)
  # initial values
  if (abs(det(S)) < smnum) {
    Lambda <- solve(S + smnum * diag(p))
  } else {
    Lambda <- solve(S) # theta
  }
  Lambda_init <- Lambda
  # from MLE
  Theta_init <- -Lambda %*% Cyx %*% solve(Cx) # theta_xy
  # tC%*% Cyx %*% solve(Cx)  true theta_xy
  iter <- 0
  Sigma_old <- S
  Theta_old <- Theta_init
  Lambda_old <- Lambda_init
  record1 <- record2 <- rep(0, EnoughIter)
  sXmatrix <- Cy %*% Lambda_old + Sigma_old %*% Theta_old %*% Cx %*% t(Theta_old)
  trsX <- sum(diag(sXmatrix))
  l1normLambda <- lam_cov * (sum(Lambda_old != 0) - p)
  l1normTheta <- lam_F * sum(Theta_old != 0)
  fx <- -log(det(Lambda_old)) + trsX + l1normLambda + l1normTheta
  while ((iter < EnoughIter) && (abs((fx - prevfx) / fx) > eps)) {
    iter <- iter + 1
    # print("Updating Lambda...")
    dat <- y + x %*% t(Theta_old) %*% Sigma_old
    S_temp <- 1 / n * t(dat) %*% dat

    fit.lasso <- glasso(S_temp, lam_cov, maxit = 1, penalize.diagonal = FALSE)
    wi.lasso <- fit.lasso$wi
    wi.scad <- wi.lasso
    epsi <- 1
    count <- 1
    while (epsi > 1e-04) {
      if (epsi < 0.001 & count > 20) {
        break
      }
      count <- count + 1
      wi.scad.old <- wi.scad
      rhomat <- scadrightderv(abs(wi.scad.old), 3.7, lam_cov)
      fit.scad <- glasso(S_temp, rhomat, maxit = 1, penalize.diagonal = FALSE)
      wi.scad <- fit.scad$wi
      epsi <- mean(abs(wi.scad - wi.scad.old))
      if (count > 50) {
        warning("scad iteration does not converge")
        break
      }
    }
    Lambda_new <- fit.scad$wi
    Sigma_new <- fit.scad$w

    obF1 <- cdF_nw(Sigma_new, Cx, Cyx, Theta_old, lam_F, 1e-3)
    Theta_new <- obF1$F
    diff_Theta <- mean(abs(Theta_new - Theta_old))
    diff_Lambda <- mean(abs(Lambda_new - Lambda_old))
    record1[iter] <- mean(abs(Theta_new - Theta_old))
    record2[iter] <- mean(abs(Lambda_new - Lambda_old))
    Theta_old <- Theta_new
    Lambda_old <- Lambda_new
    Sigma_old <- Sigma_new
  }
  list(Theta = Theta_new, Lambda = Lambda_new)
}

CallBIC <- function(y, x, Lambda, Theta, opt) {
  n <- dim(y)[1]
  p <- dim(y)[2]
  px <- dim(x)[2]
  Sy <- 1 / n * t(y) %*% y
  Syx <- 1 / n * t(y) %*% x
  Sx <- 1 / n * t(x) %*% x
  ebic.gamma <- .5
  if (opt == 1) # for glasso
    {
      z <- cbind(y, x)
      S <- cov(z)
    } else if (opt == 2) {
    dat <- y + x %*% t(Theta) %*% solve(Lambda)
    S <- 1 / n * t(dat) %*% dat
  }
  Lambda_offd <- Lambda - diag(diag(Lambda))
  sn <- sum(Lambda_offd != 0)
  if (opt == 1) {
    neglik <- -n * log(det(Lambda)) + n * sum(diag(S %*% Lambda))
    BIC <- neglik + log(n) * sn / 2
    ebic <- BIC + 4 * ebic.gamma * log(p)
  } else if (opt == 2) {
    kn <- sum(Theta != 0)
    prod <- Sy %*% Lambda + 2 * Syx %*% t(Theta) + solve(Lambda) %*% Theta %*% Sx %*% t(Theta)
    neglik <- -n * log(det(Lambda)) + n * sum(diag(prod))
    BIC <- neglik + log(n) * (dim(Theta)[1] + sn / 2 + kn)
    # BIC <- neglik + log(n)*(dim(Theta)[1] + p*px/total_df*sn/2 + p^2/total_df*kn)
    ebic <- BIC + 4 * ebic.gamma * log(p)
  }
  list(BIC = BIC, EBIC = ebic, neglik = neglik)
}

selectBICForAlg1 <- function(trainy, trainx, Lam_cov, Lam_F, N, N_j, thr, useEBIC = TRUE, method = "quic") {
  smallnum <- 1e-4
  EnoughIter <- 20
  n1 <- dim(trainy)[1]
  p <- dim(trainy)[2]
  px <- dim(trainx)[2]
  # trainx_mean <- apply(trainx,2,mean)
  # trainx <- sweep(data.matrix(trainx), 2, trainx_mean)
  if (method == "quic") {
    mod <- Alg1(trainy, trainx, smallnum, smallnum, EnoughIter, thr)
  } else if (method == "alasso") {
    mod <- Alg1_alasso(trainy, trainx, smallnum, smallnum, EnoughIter, thr)
  } else if (method == "scad") {
    mod <- Alg1_scad(trainy, trainx, smallnum, smallnum, EnoughIter, thr)
  } else {
    stop("Inappropriate penalty type")
  }
  if (useEBIC) {
    cvob <- CallBIC(trainy, trainx, mod$Lambda, mod$Theta, 2)$EBIC
  } else {
    cvob <- CallBIC(trainy, trainx, mod$Lambda, mod$Theta, 2)$BIC
  }
  opt_cv_egec <- cvob
  optlam_egec_cov <- smallnum
  optlam_egec_F <- smallnum
  optLambda_egec <- mod$Lambda
  optTheta_egec <- mod$Theta
  optSig_egec <- solve(optLambda_egec)
  bicmatrix <- matrix(NA, N + 1, N_j + 1)
  count <- 0
  for (i in (N + 1):1)
  {
    lam_cov <- Lam_cov / N * (N + 1 - i)
    for (j in (N_j + 1):1)
    {
      lam_F <- Lam_F / N_j * (N_j + 1 - j)
      # cat("current i = ", i, ", j = ", j, "\n")
      if (method == "quic") {
        mod2 <- Alg1(trainy, trainx, lam_cov, lam_F, EnoughIter, thr)
      } else if (method == "alasso") {
        mod2 <- Alg1_alasso(trainy, trainx, lam_cov, lam_F, EnoughIter, thr)
      } else if (method == "scad") {
        mod2 <- Alg1_scad(trainy, trainx, lam_cov, lam_F, EnoughIter, thr)
      } else {
        stop("Inappropriate penalty type")
      }
      if (useEBIC) {
        cvob2 <- CallBIC(trainy, trainx, mod2$Lambda, mod2$Theta, 2)$EBIC
      } else {
        cvob2 <- CallBIC(trainy, trainx, mod2$Lambda, mod2$Theta, 2)$BIC
      }
      bicmatrix[i, j] <- cvob2
      Theta2 <- mod2$Lambda
      Theta2s <- matrix(as.integer(abs(Theta2) > smallnum), nrow = p)
      cv_egec <- cvob2
      if (cv_egec < opt_cv_egec) {
        opt_cv_egec <- cv_egec
        optlam_egec_cov <- lam_cov
        optlam_egec_F <- lam_F
        optTheta_egec <- mod2$Theta
        optadj_egec <- Theta2s - diag(diag(Theta2s))
        optLambda_egec <- Theta2
      }
    }
  }
  list(
    optBIC = opt_cv_egec, optlam_cov = optlam_egec_cov, optlam_F = optlam_egec_F,
    opt_Lam = optLambda_egec, opt_Theta = optTheta_egec,
    opt_adj = optadj_egec,
    bicmatrix = bicmatrix
  )
}
gen.data.sp <- function(fit, N = 3e3) {
  fit.dm <- fit$domain
  sample.list <- NULL
  for (label in names(fit$mf)) {
    fit.pho <- fit$rho[[label]]
    xx <- seq(fit.dm[1, label], fit.dm[2, label], by = 0.01)
    tmp <- my.dssden(fit.pho, xx)
    xx.wt <- tmp / sum(tmp)
    xx.sp <- sample(xx, prob = xx.wt, replace = TRUE, size = N)
    # xx <- sample(x[, label], N, replace = TRUE)
    sample.list[[label]] <- xx.sp
  }
  return(sample.list)
}

cal.semi.mse <- function(object, delta, N = 3e3, semiflag = TRUE) {
  cat("calculating semi-mse. \n")
  fit <- object
  rho1 <- sum(object$rho.int)
  rho2 <- rho1^2 - sum(object$rho.int^2) + sum(object$rho.int2)
  s <- object$int$s
  r <- object$int$r
  s.rho <- object$int$s.rho - s * rho1
  r.rho <- object$int$r.rho - r * rho1
  int2 <- mkint2(
    object$mf, object$int$var.type,
    object$id.basis, object$quad, object$terms
  )
  ss <- int2$ss
  sr <- int2$sr
  rr <- int2$rr
  d <- object$d
  c <- object$c
  theta <- object$theta
  nq <- length(theta)
  s.eta <- ss %*% d
  r.eta <- tmp <- NULL
  r.wk <- r.rho.wk <- sr.wk <- rr.wk <- 0
  for (i in 1:nq) {
    tmp <- c(tmp, 10^(2 * theta[i]) * sum(diag(rr[, , i, i])))
    s.eta <- s.eta + 10^theta[i] * sr[, , i] %*% c
    if (length(d) == 1) {
      r.eta.wk <- sr[, , i] * d
    } else {
      r.eta.wk <- t(sr[, , i]) %*% d
    }
    r.wk <- r.wk + 10^theta[i] * r[, i]
    r.rho.wk <- r.rho.wk + 10^theta[i] * r.rho[, i]
    sr.wk <- sr.wk + 10^theta[i] * sr[, , i]
    for (j in 1:nq) {
      r.eta.wk <- r.eta.wk + 10^theta[j] * rr[, , i, j] %*% c
      rr.wk <- rr.wk + 10^(theta[i] + theta[j]) * rr[, , i, j]
    }
    r.eta <- cbind(r.eta, r.eta.wk)
  }
  s.eta <- s.eta - s * (sum(s * d) + sum(r.wk * c))
  r.eta <- r.eta - r * (sum(s * d) + sum(r.wk * c))
  ss <- ss - outer(s, s, "*")
  sr.wk <- sr.wk - outer(s, r.wk, "*")
  rr.wk <- rr.wk - outer(r.wk, r.wk, "*")
  rho.eta <- sum(s.rho * d) + sum(r.rho.wk * c)
  eta2 <- sum(c * (rr.wk %*% c)) + sum(d * (ss %*% d)) + 2 * sum(d * (sr.wk %*% c))
  old.mse <- eta2 + rho2 - rho1^2 + 2 * rho.eta
  # cat("old mse is", old.mse, "\n")


  if (semiflag) {
    # generate data frame for monte-carlo simulation
    sample.list <- gen.data.sp(fit, N = N)
    # generate components needed to calculate mse and se
    delta.x <- delta.x2 <- B_21 <- B_22 <- A_31 <- A_32 <- eta.hat <- eta.hat2 <- rep(NaN, N)
    eta.mu.record <- rep(NaN, N)
    s.x.full <- r.x.full <- NULL
    s.eta.mc.full <- r.eta.mc.full <- NULL

    for (i in 1:N) {
      data.pt <- NULL
      eta.mu.tmp <- eta.mu <- NULL
      for (label in names(sample.list)) {
        fit.rho <- fit$rho[[label]]
        data.pt <- c(data.pt, sample.list[[label]][i])
        eta.mu.tmp <- c(eta.mu.tmp, my.dssden(fit.rho, sample.list[[label]][i]))
        eta.mu <- sum(log(eta.mu.tmp))
        eta.mu.record[[i]] <- eta.mu
      }
      delta.x[i] <- t(data.pt) %*% delta %*% data.pt
      delta.x2[i] <- delta.x[i]^2
      data.df <- as.data.frame(t(data.pt))
      colnames(data.df) <- names(sample.list)
      multi_density.info <- my.dssden(fit, data.df)
      eta.hat[i] <- multi_density.info$pdf.unscaled
      eta.hat2[i] <- eta.hat[i]^2
      # s.x.full[[i]] <- multi_density.info$s
      # r.x.full[[i]] <- multi_density.info$r
      # s.eta.mc.full[[i]] <- s.x.full[[i]] * eta.hat[[i]]
      # r.eta.mc.full[[i]] <- r.x.full[[i]] * eta.hat[[i]]
      A_31[i] <- (eta.hat[i] + delta.x[i])
      A_32[i] <- (eta.hat[i] + delta.x[i])^2
      B_21[i] <- (eta.hat[i] + eta.mu) * delta.x[i]
      B_22[i] <- (eta.hat[i] + eta.mu)
    }
    # A-B verifies the correctness of this mc approximation. it's close to old.mse
    A <- 1 / N * (sum(eta.hat2) + sum(eta.mu.record^2) + 2 * sum(eta.mu.record * eta.hat))
    B <- (1 / N * sum(eta.hat + eta.mu.record))^2
    A - B
    A_2 <- 1 / N * sum(delta.x2) - (1 / N * sum(delta.x))^2
    B_2 <- 1 / N * sum(B_21) - (1 / N * sum(B_22)) * (1 / N * sum(delta.x))
    semi.mse <- old.mse + A_2 + 2 * B_2
    zeta2.int <- 1 / N * sum(A_32) - (1 / N * sum(A_31))^2
    list(
      old.mse = old.mse, semi.mse = semi.mse, zeta2.int = zeta2.int, sample.list = sample.list,
      eta.hat = eta.hat, delta.x = delta.x
    )
  }
}

my.project.mc.v2 <- function(object, delta, semiinfo, include, drop1 = FALSE, N = 3e3, semiflag = TRUE) {
  old.mse <- semiinfo$old.mse
  semi.mse <- semiinfo$semi.mse
  zeta2.int <- semiinfo$zeta2.int
  sample.list <- semiinfo$sample.list
  eta.hat <- semiinfo$eta.hat
  delta.x <- semiinfo$delta.x
  s.delta <- r.delta <- NULL
  fit <- object
  rho1 <- sum(object$rho.int)
  rho2 <- rho1^2 - sum(object$rho.int^2) + sum(object$rho.int2)
  s <- object$int$s
  r <- object$int$r
  s.rho <- object$int$s.rho - s * rho1
  r.rho <- object$int$r.rho - r * rho1
  int2 <- mkint2(
    object$mf, object$int$var.type,
    object$id.basis, object$quad, object$terms
  )
  ss <- int2$ss
  sr <- int2$sr
  rr <- int2$rr
  d <- object$d
  c <- object$c
  theta <- object$theta
  nq <- length(theta)
  s.eta <- ss %*% d
  r.eta <- tmp <- NULL
  r.wk <- r.rho.wk <- sr.wk <- rr.wk <- 0
  for (i in 1:nq) {
    tmp <- c(tmp, 10^(2 * theta[i]) * sum(diag(rr[, , i, i])))
    s.eta <- s.eta + 10^theta[i] * sr[, , i] %*% c
    if (length(d) == 1) {
      r.eta.wk <- sr[, , i] * d
    } else {
      r.eta.wk <- t(sr[, , i]) %*% d
    }
    r.wk <- r.wk + 10^theta[i] * r[, i]
    r.rho.wk <- r.rho.wk + 10^theta[i] * r.rho[, i]
    sr.wk <- sr.wk + 10^theta[i] * sr[, , i]
    for (j in 1:nq) {
      r.eta.wk <- r.eta.wk + 10^theta[j] * rr[, , i, j] %*% c
      rr.wk <- rr.wk + 10^(theta[i] + theta[j]) * rr[, , i, j]
    }
    r.eta <- cbind(r.eta, r.eta.wk)
  }
  s.eta <- s.eta - s * (sum(s * d) + sum(r.wk * c))
  r.eta <- r.eta - r * (sum(s * d) + sum(r.wk * c))
  ss <- ss - outer(s, s, "*")
  sr.wk <- sr.wk - outer(s, r.wk, "*")
  rr.wk <- rr.wk - outer(r.wk, r.wk, "*")
  rho.eta <- sum(s.rho * d) + sum(r.rho.wk * c)
  eta2 <- sum(c * (rr.wk %*% c)) + sum(d * (ss %*% d)) + 2 * sum(d * (sr.wk %*% c))

  # se
  ## calculate projection
  rkl <- function(include) {
    cat("include: ", include, "\n")
    cat("calculating zeta \n")
    # verify correctness of mu using mc approximation
    inc.wk <- union(names(object$mf), include)
    id.s <- id.q <- NULL
    for (label in inc.wk) {
      if (!any(label == object$terms$labels)) next
      term <- object$terms[[label]]
      if (term$nphi > 0) id.s <- c(id.s, term$iphi + (1:term$nphi) - 2)
      if (term$nrk > 0) id.q <- c(id.q, term$irk + (1:term$nrk) - 1)
    }
    ss.wk <- ss[id.s, id.s]
    r.eta.wk <- r.wk <- sr.wk <- rr.wk <- 0
    for (i in id.q) {
      r.eta.wk <- r.eta.wk + 10^theta[i] * r.eta[, i]
      r.wk <- r.wk + 10^theta[i] * r[, i]
      sr.wk <- sr.wk + 10^theta[i] * sr[id.s, , i]
      for (j in id.q) {
        rr.wk <- rr.wk + 10^(theta[i] + theta[j]) * rr[, , i, j]
      }
    }
    sr.wk <- sr.wk - outer(s[id.s], r.wk, "*")
    rr.wk <- rr.wk - outer(r.wk, r.wk, "*")
    v <- cbind(rbind(ss.wk, t(sr.wk)), rbind(sr.wk, rr.wk))
    if (semiflag) {
      form <- as.formula(paste("~", paste(inc.wk, collapse = " + ")))
      fit.reduced <- ssden1(form,
        domain = object$domain,
        id.basis = object$id.basis, data = object$mf
      )
      fit.reduced$theta <- object$theta[id.q]
      s.eta.mc <- r.eta.mc <- s.x <- r.x <- NULL
      for (i in 1:N) {
        data.pt <- NULL
        for (label in names(sample.list)) {
          data.pt <- c(data.pt, sample.list[[label]][i])
        }
        data.df <- as.data.frame(t(data.pt))
        colnames(data.df) <- names(sample.list)
        reduced.info <- my.dssden(fit.reduced, data.df)
        s.x[[i]] <- reduced.info$s
        r.x[[i]] <- reduced.info$r
        s.eta.mc[[i]] <- s.x[[i]] * eta.hat[[i]]
        r.eta.mc[[i]] <- r.x[[i]] * eta.hat[[i]]
        s.delta[[i]] <- s.x[[i]] * delta.x[i]
        r.delta[[i]] <- r.x[[i]] * delta.x[i]
      }
      mu_1 <- 1 / N * Reduce("+", s.eta.mc) - (1 / N * sum(eta.hat)) * (1 / N * Reduce("+", s.x))
      mu_2 <- 1 / N * Reduce("+", r.eta.mc) - (1 / N * sum(eta.hat)) * (1 / N * Reduce("+", r.x))
      s.delta.int <- 1 / N * Reduce("+", s.delta) - (1 / N * sum(delta.x)) * (1 / N * Reduce("+", s.x))
      r.delta.int <- 1 / N * Reduce("+", r.delta) - (1 / N * sum(delta.x)) * (1 / N * Reduce("+", r.x))
      mu_mc <- c(mu_1 + s.delta.int, mu_2 + r.delta.int)
      mu <- mu_mc
      nn <- length(mu)
      # z <- chol(v,pivot=TRUE)
      # v.ori <-v
      # v <- z
      # rkv <- attr(z,"rank")
      # m.eps <- .Machine$double.eps
      # while (v[rkv,rkv]<2*sqrt(m.eps)*v[1,1]) rkv <- rkv - 1
      # if (rkv<nn) v[(1:nn)>rkv,(1:nn)>rkv] <- diag(v[1,1],nn-rkv)
      # mu <- backsolve(v,mu[attr(z,"pivot")],transpose=TRUE)
      # se <- eta2 - sum(mu[1:rkv]^2)
      dc.tilda <- solve(v) %*% mu
      d.tilda <- dc.tilda[1:length(fit.reduced$d)]
      c.tilda <- dc.tilda[(length(fit.reduced$d) + 1):length(dc.tilda)]
      zeta.tilda <- sum(c.tilda * (rr.wk %*% c.tilda)) + sum(d.tilda * (ss.wk %*% d.tilda)) + 2 * sum(d.tilda * (sr.wk %*% c.tilda))
      se <- zeta2.int - zeta.tilda
    }
    # A_3 <- 1/N * sum(A_31) - (1/N * sum(eta.hat + delta.x))^2
    # A_3 <- 1/N * sum(A_31) - (1/N * sum(eta.hat + delta.x))^2
    else {
      mu <- c(s.eta[id.s], r.eta.wk)
      nn <- length(mu)
      z <- chol(v, pivot = TRUE)
      v.ori <- v
      v <- z
      rkv <- attr(z, "rank")
      m.eps <- .Machine$double.eps
      # cat("rkv value is: ", rkv, "\n")
      while (v[rkv, rkv] < 2 * sqrt(m.eps) * v[1, 1]) rkv <- rkv - 1
      # cat("now rkv value is: ", rkv, "\n")
      if (rkv < nn) v[(1:nn) > rkv, (1:nn) > rkv] <- diag(v[1, 1], nn - rkv)
      mu <- backsolve(v, mu[attr(z, "pivot")], transpose = TRUE)
      # print("print here")
      se <- eta2 - sum(mu[1:rkv]^2)
    }
  }

  ## projection
  # cat("returning ratios")
  if (drop1) {
    se <- NULL
    for (i in 1:length(include)) se <- c(se, rkl(include[-i]))
    if (semiflag) {
      ratio <- se / semi.mse
    } else {
      ratio <- se / old.mse
      cat("ratio: ", ratio)
    }
    names(se) <- names(ratio) <- include
  } else {
    se <- rkl(include)
  }
  if (semiflag) {
    ratio <- se / semi.mse
  } else {
    ratio <- ratio <- se / old.mse
  }
  list(ratio = ratio, se = se)
}

# Non-parametric
ssX <- function(x, delta, semiflag = TRUE, cutoff = 0.05, N) {
  px <- dim(x)[2]
  n <- dim(x)[1]
  mn <- apply(x, 2, min)
  mx <- apply(x, 2, max)
  numterm <- px + px * (px - 1) / 2
  domain <- data.frame(rbind(mn, mx))
  x <- data.frame(x)
  fit.x <- ssden1(~ (.)^2, domain = domain, data = x)
  # print("Full ssden1 model fitted")
  label <- fit.x$terms$labels[(px + 1):numterm]


  ## use semi-parametric projection
  semi.info <- cal.semi.mse(fit.x, delta, N = N, semiflag = semiflag)
  ratio <- my.project.mc.v2(fit.x, delta, semi.info, include = label, drop1 = TRUE, N = N, semiflag = semiflag)$ratio
  # ratio <- project(fit.x,include=label,drop1=TRUE)$ratio
  cat("Projection compared to only main effects model calculated \n")
  ord <- rev(order(ratio))
  po <- length(ord)
  if (semiflag) {
    temp <- my.project.mc.v2(fit.x, delta, semi.info, include = 0, drop1 = FALSE, N = N, semiflag = semiflag)$ratio
    if (temp < cutoff) {
      cat("ratio with only main effects: ", temp, "\n")
      return("No interaction")
    }
    for (i in 1:po) {
      cat("Doing projection at i = ", i, "\n")
      lis <- ord[1:i]
      ratio <- my.project.mc.v2(fit.x, delta, semi.info, include = label[lis], drop1 = FALSE, N = N, semiflag = semiflag)$ratio
      cat("current term: ", label[lis], ", current ratio: ", ratio, "\n")
      if (ratio < cutoff) {
        return(label[lis])
      }
    }
  } else {
    if (project(fit.x, include = 0)$ratio < cutoff) {
      return("No interaction")
    }
    for (i in 1:po) {
      print(paste("Doing projection at i = ", i, "\n"))
      lis <- ord[1:i]
      ratio <- project(fit.x, include = label[lis])$ratio
      if (ratio < cutoff) {
        return(label[lis])
      }
    }
  }
}

# translate ssden1 result into precision matrix of x
XandX <- function(Xr, CXX) {
  px <- dim(CXX)[2]
  CXXr <- CXX
  l <- length(Xr)
  if (Xr[1] != "No interaction") {
    for (k in 1:l) {
      for (i in 1:px) {
        for (j in 1:px) {
          quo <- paste("X", i, ":X", j, sep = "")
          t <- as.character(Xr[k])
          quo <- as.character(quo)
          if (quo == t) {
            CXXr[i, j] <- 1
          }
          quo2 <- paste("X", j, ":X", i, sep = "")
          quo2 <- as.character(quo2)
          if (quo2 == t) {
            CXXr[j, i] <- 1
          }
        }
      }
    }
  }
  CXXr
}
mkint2 <- function(mf, type, id.basis, quad, term) {
  ## Obtain model terms
  mt <- attr(mf, "terms")
  xvars <- as.character(attr(mt, "variables"))[-1]
  xfacs <- attr(mt, "factors")
  term.labels <- labels(mt)
  vlist <- xvars[as.logical(apply(xfacs, 1, sum))]
  ## Create phi and rk
  nbasis <- length(id.basis)
  phi.term <- rk.term <- list(NULL)
  nvar <- length(names(mf))
  ns <- nq <- 0
  for (label in term.labels) {
    ns <- ns + term[[label]]$nphi
    nq <- nq + term[[label]]$nrk
    phi.term[[label]] <- rk.term[[label]] <- list(NULL)
    vlist <- xvars[as.logical(xfacs[, label])]
    x <- mf[, vlist]
    dm <- length(vlist)
    phi <- rk <- NULL
    if (dm == 1) {
      type.wk <- type[[vlist]][[1]]
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
        range <- type[[vlist]][[2]]
        ## phi
        phi.fun <- mkphi.cubic(range)
        phi <- phi.fun$fun(xmesh, 1, phi.fun$env)
        ## rk
        rk.fun <- mkrk.cubic(range)
        rk <- rk.fun$fun(xmesh, xx, rk.fun$env, TRUE)
      }
      if (type.wk %in% c("cubic.per", "linear", "linear.per", "sphere")) {
        ## cubic periodic, linear, and linear periodic splines
        range <- type[[vlist]][[2]]
        ## rk
        if (type.wk == "cubic.per") rk.fun <- mkrk.cubic.per(range)
        if (type.wk == "linear") rk.fun <- mkrk.linear(range)
        if (type.wk == "linear.per") rk.fun <- mkrk.linear.per(range)
        if (type.wk == "sphere") rk.fun <- mkrk.sphere(range)
        rk <- rk.fun$fun(xmesh, xx, rk.fun$env, TRUE)
      }
      if (type.wk == "tp") {
        ## thin-plate splines
        par <- type[[vlist]][[2]]
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
        par <- type[[vlist]][[2]]
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
      phi.term[[label]][[vlist]] <- phi
      if (is.null(rk)) {
        rk.term[[label]][[vlist]] <- rk
      } else {
        nmesh <- length(quad[[vlist]]$wt)
        rk.term[[label]][[vlist]] <- array(rk, c(nmesh, nbasis, 1))
      }
    } else {
      bin.fac <- n.phi <- phi.list <- rk.list <- NULL
      for (i in 1:dm) {
        type.wk <- type[[vlist[i]]][[1]]
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
          range <- type[[vlist[i]]][[2]]
          ## phi
          phi.wk <- mkphi.cubic(range)
          n.phi <- c(n.phi, 1)
          ## rk
          rk.wk <- mkrk.cubic(range)
          bin.fac <- c(bin.fac, 0)
        }
        if (type.wk %in% c("cubic.per", "linear", "linear.per", "sphere")) {
          ## cubic periodic, linear, or linear periodic splines
          range <- type[[vlist[i]]][[2]]
          n.phi <- c(n.phi, 0)
          phi.wk <- NULL
          if (type.wk == "cubic.per") rk.wk <- mkrk.cubic.per(range)
          if (type.wk == "linear") rk.wk <- mkrk.linear(range)
          if (type.wk == "linear.per") rk.wk <- mkrk.linear.per(range)
          if (type.wk == "sphere") rk.wk <- mkrk.sphere(range)
          bin.fac <- c(bin.fac, 0)
        }
        if (type.wk == "tp") {
          ## thin-plate splines
          par <- type[[vlist[i]]][[2]]
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
          par <- type[[vlist[i]]][[2]]
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
            phi.term[[label]][[vlist[i]]] <-
              cbind(phi.term[[label]][[vlist[i]]], phi)
          }
        }
      }
      ## rk
      n.rk <- ifelse(n.phi, 2, 1)
      nrk <- prod(n.rk) - as.logical(nphi)
      if (nrk > 0) {
        for (nu in 1:nrk) {
          ind <- nu - !nphi
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
            nmesh <- length(quad[[vlist[i]]]$wt)
            rk.term[[label]][[vlist[i]]] <-
              array(
                c(rk.term[[label]][[vlist[i]]], rk),
                c(nmesh, nbasis, nu)
              )
          }
        }
      }
    }
  }
  ## create arrays
  ss <- matrix(1, ns, ns)
  sr <- array(1, c(ns, nbasis, nq))
  rr <- array(1, c(nbasis, nbasis, nq, nq))
  for (label1 in term.labels) {
    if (!term[[label1]]$nphi) {
      id.s1 <- NULL
    } else {
      id.s1 <- term[[label1]]$iphi + (1:term[[label1]]$nphi) - 2
    }
    if (!term[[label1]]$nrk) {
      id.r1 <- NULL
    } else {
      id.r1 <- term[[label1]]$irk + (1:term[[label1]]$nrk) - 1
    }
    irk1 <- term[[label1]]$irk
    for (label2 in term.labels) {
      if (!term[[label2]]$nphi) {
        id.s2 <- NULL
      } else {
        id.s2 <- term[[label2]]$iphi + (1:term[[label2]]$nphi) - 2
      }
      if (!term[[label2]]$nrk) {
        id.r2 <- NULL
      } else {
        id.r2 <- term[[label2]]$irk + (1:term[[label2]]$nrk) - 1
      }
      irk2 <- term[[label2]]$irk
      for (xlab in names(mf)) {
        wmesh <- quad[[xlab]]$wt
        phi1 <- phi.term[[label1]][[xlab]]
        phi2 <- phi.term[[label2]][[xlab]]
        rk1 <- rk.term[[label1]][[xlab]]
        rk2 <- rk.term[[label2]][[xlab]]
        ## ss
        if (!is.null(id.s1) & !is.null(id.s2)) {
          if ((!is.null(phi1)) & (!is.null(phi2))) {
            ss[id.s1, id.s2] <- ss[id.s1, id.s2] * (t(wmesh * phi1) %*% phi2)
          } else {
            if (!is.null(phi1)) {
              ss[id.s1, id.s2] <- ss[id.s1, id.s2] * apply(wmesh * matrix(phi1), 2, sum)
            } else {
              if (!is.null(phi2)) {
                ss[id.s1, id.s2] <- t(t(ss[id.s1, id.s2]) *
                  apply(wmesh * matrix(phi2), 2, sum))
              }
            }
          }
        }
        ## sr
        if (!is.null(id.s1) & !is.null(id.r2)) {
          if ((!is.null(phi1)) & (!is.null(rk2))) {
            for (i in id.r2) {
              sr[id.s1, , i] <- sr[id.s1, , i] * (t(wmesh * phi1) %*% rk2[, , i - irk2 + 1])
            }
          } else {
            if (!is.null(phi1)) {
              sr[id.s1, , id.r2] <- sr[id.s1, , id.r2] * apply(wmesh * matrix(phi1), 2, sum)
            } else {
              if (!is.null(rk2)) {
                for (i in id.r2) {
                  sr[id.s1, , i] <- t(t(sr[id.s1, , i]) *
                    apply(wmesh * rk2[, , i - irk2 + 1], 2, sum))
                }
              }
            }
          }
        }
        ## rr
        if (!is.null(id.r1) & !is.null(id.r2)) {
          if ((!is.null(rk1)) & (!is.null(rk2))) {
            for (i in id.r1) {
              for (j in id.r2) {
                rr[, , i, j] <- rr[, , i, j] * (t(wmesh * rk1[, , i - irk1 + 1]) %*% rk2[, , j - irk2 + 1])
              }
            }
          } else {
            if (!is.null(rk1)) {
              for (i in id.r1) {
                rr[, , i, id.r2] <- rr[, , i, id.r2] * apply(wmesh * rk1[, , i - irk1 + 1], 2, sum)
              }
            } else {
              if (!is.null(rk2)) {
                for (j in id.r2) {
                  rr[, , id.r1, j] <-
                    aperm(aperm(rr[, , id.r1, j, drop = FALSE], c(2, 1, 3, 4)) *
                      apply(wmesh * rk2[, , j - irk2 + 1], 2, sum), c(2, 1, 3, 4))
                }
              }
            }
          }
        }
      }
    }
  }
  list(ss = ss, sr = sr, rr = rr)
}

my.d.ssden1 <- ## Evaluate density estimate
  function(object, x) {
    if (!("ssden1" %in% class(object))) stop("gss error in d.ssden1: not a ssden1 object")
    ## rho
    rho <- 1
    for (xlab in names(object$mf)) {
      xx <- x[[xlab]]
      rho.wk <- object$rho[[xlab]]
      if (is.factor(xx)) rho <- rho * rho.wk[xx]
      if (is.vector(xx) & !is.factor(xx)) rho <- rho * dssden(rho.wk, xx)
      if (is.matrix(xx)) rho <- rho * dssden(rho.wk, xx)
    }
    ## exp(eta)
    s <- NULL
    r <- matrix(0, dim(x)[1], length(object$id.basis))
    nq <- 0
    for (label in object$terms$labels) {
      xx <- object$mf[object$id.basis, object$terms[[label]]$vlist]
      x.new <- x[, object$terms[[label]]$vlist]
      nphi <- object$terms[[label]]$nphi
      nrk <- object$terms[[label]]$nrk
      if (nphi) {
        phi <- object$terms[[label]]$phi
        for (i in 1:nphi) {
          s <- cbind(s, phi$fun(x.new, nu = i, env = phi$env))
        }
      }
      if (nrk) {
        rk <- object$terms[[label]]$rk
        for (i in 1:nrk) {
          nq <- nq + 1
          r <- r + 10^object$theta[nq] * rk$fun(x.new, xx, nu = i, env = rk$env, out = TRUE)
        }
      }
    }
    # unscaled version
    pdf.unscaled <- as.vector(cbind(s, r) %*% c(object$d, object$c))
    list(s = s, r = r, pdf.unscaled = pdf.unscaled)
  }
d.ssden <- ## Evaluate density estimate
  function(object, x) {
    if (class(object) != "ssden") stop("gss error in d.ssden: not a ssden object")
    if (dim(object$mf)[2] == 1 & is.vector(x)) {
      x <- data.frame(x)
      colnames(x) <- colnames(object$mf)
    }
    s <- NULL
    r <- matrix(0, dim(x)[1], length(object$id.basis))
    nq <- 0
    for (label in object$terms$labels) {
      xx <- object$mf[object$id.basis, object$terms[[label]]$vlist]
      x.new <- x[, object$terms[[label]]$vlist]
      nphi <- object$terms[[label]]$nphi
      nrk <- object$terms[[label]]$nrk
      if (nphi) {
        phi <- object$terms[[label]]$phi
        for (i in 1:nphi) {
          s <- cbind(s, phi$fun(x.new, nu = i, env = phi$env))
        }
      }
      if (nrk) {
        rk <- object$terms[[label]]$rk
        for (i in 1:nrk) {
          nq <- nq + 1
          r <- r + 10^object$theta[nq] * rk$fun(x.new, xx, nu = i, env = rk$env, out = TRUE)
        }
      }
    }
    as.vector(exp(cbind(s, r) %*% c(object$d, object$c)) / object$int)
  }

my.dssden <- ## Evaluate density estimate
  function(object, x) {
    ## check input
    if (!("ssden" %in% class(object))) stop("gss error in dssden: not a ssden object")
    if ("ssden1" %in% class(object)) {
      return(my.d.ssden1(object, x))
    } else {
      return(d.ssden(object, x))
    }
  }
