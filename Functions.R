library(tensor)
library(LaplacesDemon)
library(foreach)
library(doSNOW)
library(tidyr)
library(data.table)
library(mclust)
library(zipfR)
library(withr)

Pars_init.paral <- function(X, k, nstartR = 100, nThreads = 1, density = "MVN") {
  r_Pars_init <- function(X, k, nstartR = 100, nThreads = 1, density) {
    comb <- function(x, ...) {
      lapply(
        seq_along(x),
        function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
      )
    }
    dMVnorm <- function(X, M, U, V) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)

      return(pdf)
    }
    tr <- function(x) {
      return(sum(diag(x)))
    }
    dMVtin <- function(X, M, U, V, theta) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      intf <- function(w, del) {
        w^((p * r) / 2) * exp((-w / 2) * del)
      }
      pdfinteg <- sapply(1:num, function(i) stats::integrate(intf, lower = (1 - theta), upper = 1, del = delta[i])$value)
      pdfconst <- (2 * pi)^(-(p * r) / 2) * (1 / theta) * det(U)^(-r / 2) * det(V)^(-p / 2)
      PDF <- pdfconst * pdfinteg

      return(PDF)
    }

    # Dimensions

    p <- dim(X)[1] # rows of each matrix;
    r <- dim(X)[2] # columns of each matrix;
    num <- dim(X)[3] # sample size;

    # Create some objects

    prior <- numeric(k)
    M <- array(0, dim = c(p, r, k))
    sigmaU <- array(0, dim = c(p, p, k))
    sigmaV <- array(0, dim = c(r, r, k))
    nu <- numeric(k)

    WR <- array(0, dim = c(p, p, k))
    WC <- array(0, dim = c(r, r, k))
    tempWR <- array(0, dim = c(p, p, num))
    tempWC <- array(0, dim = c(r, r, num))
    tempM <- array(0, dim = c(p, r, num))

    post <- dens <- array(0, c(num, k), dimnames = list(1:(num), paste("comp.", 1:k, sep = "")))

    ## Random initialization ##

    eu <- matrix(0, nrow = num, ncol = k)
    classy <- numeric(num)
    rand.start <- matrix(0, nstartR, k)
    nu.start <- matrix(0, nstartR, k)

    withr::with_seed(1, for (i in 1:nstartR) {
      rand.start[i, ] <- sample(c(1:num), k)
      if (density == "MVTIN") {
        nu.start[i, ] <- runif(k, 0.6, 0.95)
      }
    })

    cluster <- makeCluster(nThreads, type = "SOCK")
    registerDoSNOW(cluster)

    pb <- txtProgressBar(max = nstartR, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    oper0 <- foreach(l = 1:nstartR, .combine = "comb", .packages = c("tensor"), .multicombine = TRUE, .init = list(list(), list(), list(), list(), list(), list(), list()), .options.snow = opts) %dopar% {

      ### part 0 ###

      sec <- rand.start[l, ]
      nu <- nu.start[l, ]

      for (j in 1:k) {
        M[, , j] <- X[, , sec[j]]
      }

      for (j in 1:k) {
        for (i in 1:(num)) {
          eu[i, j] <- norm((X[, , i] - M[, , j]), type = "F")
        }
      }

      for (i in 1:(num)) {
        classy[i] <- which.min(eu[i, ])
      }

      z <- mclust::unmap(classy)

      ### part 1 ###

      for (j in 1:k) {
        M[, , j] <- rowSums(X * z[, j][slice.index(X, 3)], dims = 2) / sum(z[, j])

        pt1 <- sweep(X, 1:2, M[, , j]) * z[, j][slice.index(sweep(X, 1:2, M[, , j]), 3)]

        pt2 <- aperm(sweep(X, 1:2, M[, , j]), c(2, 1, 3))

        WR[, , j] <- tensor(pt1, pt2, c(2, 3), c(1, 3))
      }

      phi <- tr(rowSums(WR, dims = 2)) / ((num) * p * r)

      for (j in 1:k) {
        sigmaU[, , j] <- phi * diag(1, p, p)
        sigmaV[, , j] <- diag(1, r, r)

        if (density == "MVN") {
          dens[, j] <- dMVnorm(X = X, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j])
        } else if (density == "MVTIN") {
          dens[, j] <- dMVtin(X = X, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j], theta = nu[j])
        }
      }

      if (k == 1) {
        prior <- 1
      } else {
        prior <- colMeans(z)
      }

      ### part 2 ###

      # mixture density

      numerator <- matrix(rep(prior, num), num, k, byrow = TRUE) * dens
      mixt.dens <- rowSums(numerator)
      llk <- sum(log(mixt.dens))

      # return

      list(sigmaU, sigmaV, M, prior, llk, classy, nu)
    }

    stopCluster(cluster)
    registerDoSEQ()
    close(pb)

    df <- data.frame(llk = unlist(oper0[[5]]), pos = c(1:nstartR))
    df <- df %>% drop_na()
    df <- df[!is.infinite(rowSums(df)), ]

    bestR <- head(setorderv(df, cols = "llk", order = -1), n = 1)$pos

    res <- vector(mode = "list", length = 7)

    for (i in 1:7) {
      res[[i]] <- oper0[[i]][bestR]
    }

    return(res)
  }

  results <- vector(mode = "list", length = length(k))

  for (g in 1:length(k)) {
    print(paste(paste("Initializing Parsimonious Matrix", density), paste("mixtures with k =", k[g])))

    results[[g]] <- r_Pars_init(X = X, k = k[g], nstartR = nstartR, nThreads = nThreads, density = density)
  }

  return(results)
}

Pars_mixt.paral <- function(X, k, init.par = NULL, mod.row = "all", mod.col = "all", mod.theta = "all", nThreads = 1, density = "MVN") {
  if (all(mod.row == "all")) {
    model.row <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV")
  } else {
    model.row <- mod.row
  }
  if (all(mod.col == "all")) {
    model.col <- c("II", "EI", "VI", "EE", "VE", "EV", "VV")
  } else {
    model.col <- mod.col
  }

  if (density == "MVN") {
    pt.mod <- length(model.row) * length(model.col)
    list.mod <- expand.grid(model.row, model.col)
  } else {
    if (all(mod.theta == "all")) {
      model.theta <- c("E", "V")
    } else {
      model.theta <- mod.theta
    }

    pt.mod <- length(model.row) * length(model.col) * length(model.theta)
    list.mod <- expand.grid(model.row, model.col, model.theta)
  }

  tol <- 0.001
  tol2 <- 0.001
  maxit = 500
  maxit2 = 100
  
  comb <- function(x, ...) {
    lapply(
      seq_along(x),
      function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }
  oper <- vector(mode = "list", length = length(k))

  Pars_mixt <- function(X, k, init.par = NULL, mod.row = NULL, mod.col = NULL, mod.theta = NULL, tol = 0.001, tol2 = 0.001, maxit = 500, maxit2 = 100, density = "MVN") {
    ptm <- proc.time()

    # Fuctions

    tr <- function(x) {
      return(sum(diag(x)))
    }
    dMVnorm <- function(X, M, U, V) {
      tr <- function(x) {
        return(sum(diag(x)))
      }

      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)

      return(pdf)
    }
    dMVtin <- function(X, M, U, V, theta) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      intf <- function(w, del) {
        w^((p * r) / 2) * exp((-w / 2) * del)
      }
      pdfinteg <- sapply(1:num, function(i) stats::integrate(intf, lower = (1 - theta), upper = 1, del = delta[i])$value)
      pdfconst <- (2 * pi)^(-(p * r) / 2) * (1 / theta) * det(U)^(-r / 2) * det(V)^(-p / 2)
      PDF <- pdfconst * pdfinteg

      return(PDF)
    }
    Mstep_AECM <- function(X, k, weights = NULL, M, sigmaU, sigmaV, mod.theta) {
      n <- dim(X)[3]
      if (is.null(weights)) {
        weights <- rep(1, n)
      }

      if (mod.theta == "V") {
        f1 <- function(par, weights, X, M, sigmaU, sigmaV) {
          theta <- par
          pll <- sum(weights * log(dMVtin(X = X, M = M, U = sigmaU, V = sigmaV, theta = theta)))
          return(pll)
        }
        res <- stats::optimize(
          f = f1, interval = c(0, 1), weights = weights,
          X = X, M = M, sigmaU = sigmaU, sigmaV = sigmaV, maximum = TRUE
        )
        theta <- res$maximum
      }

      if (mod.theta == "E") {
        f2 <- function(par, weights, X, M, k, sigmaU, sigmaV) {
          theta <- par
          pll <- 0

          for (j in 1:k) {
            pll <- pll + sum(weights[, j] * log(dMVtin(X = X, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j], theta = theta)))
          }

          return(pll)
        }

        res <- stats::optimize(
          f = f2, interval = c(0, 1), weights = weights,
          X = X, k = k, M = M, sigmaU = sigmaU, sigmaV = sigmaV, maximum = TRUE
        )

        theta <- res$maximum
      }

      return(theta)
    }

    # Dimensions

    p <- dim(X)[1] # rows of each matrix;
    r <- dim(X)[2] # columns of each matrix;
    num <- dim(X)[3] # sample size;

    # Create objects related to M and prior

    prior <- numeric(k)
    M <- array(0, dim = c(p, r, k))
    tempM <- array(0, dim = c(p, r, num))

    ## Objects related to the two covariance matrices

    sigmaU <- array(0, dim = c(p, p, k))
    sigmaV <- array(0, dim = c(r, r, k))

    WR <- array(0, dim = c(p, p, k))
    WC <- array(0, dim = c(r, r, k))
    tempWR <- array(0, dim = c(p, p, num))
    tempWC <- array(0, dim = c(r, r, num))

    phi.k <- numeric(k)
    temp.phi <- vector("list", k) # for VEI, VEE, VEV
    temp.phi2 <- numeric(k) # for EVI
    temp.numdeltaR <- array(0, dim = c(p, p, k)) # for VEI, VEE, VEV
    deltaU.k <- array(0, dim = c(p, p, k))
    deltaV.k <- array(0, dim = c(r, r, k))
    ftemp.r <- array(0, dim = c(p, p, k)) # MM object
    ftemp.c <- array(0, dim = c(r, r, k)) # MM object
    numphi <- numeric(k) # for EVE, EVV
    gammaU.k <- array(0, dim = c(p, p, k))
    gammaV.k <- array(0, dim = c(r, r, k))
    tempW_EEV <- vector("list", k) # for EEV, VEV
    tempW_EV <- vector("list", k) # for EV
    tempomegaR <- array(0, dim = c(p, p, k)) # for EEV, VEV
    tempomegaC <- array(0, dim = c(r, r, k)) # for EEV, VEV
    V_EVV.U.K <- array(0, dim = c(p, p, k)) # for EVV

    ## Other objects

    theta <- numeric(k)
    post <- dens <- array(0, c(num, k), dimnames = list(1:(num), paste("comp.", 1:k, sep = "")))
    w <- matrix(0, nrow = num, ncol = k)

    # Preliminary definition of convergence criterions for EM/MM algorithms

    check <- 0
    check2 <- 0
    loglik.old <- -Inf
    loglik.new <- NULL
    ll <- NULL
    mark <- 1  
    MM.r.old <- -Inf
    MM.c.old <- -Inf
    m.iter <- 0
    m.iter2 <- 0
    cnt <- 0

    ### Algorithm ###

    oper0 <- init.par
    sigmaU <- oper0[[1]][[1]]
    sigmaV <- oper0[[2]][[1]]
    M <- oper0[[3]][[1]]
    prior <- oper0[[4]][[1]]
    classy <- oper0[[6]][[1]]
    theta <- oper0[[7]][[1]]

    if (mod.row == "VEI" | mod.row == "VEE" | mod.row == "VEV") {
      for (j in 1:k) {
        temp.phi[[j]] <- eigen(sigmaU[, , j])$values

        phi.k[j] <- (prod(temp.phi[[j]]))^(1 / p)
      }
    }
    if (mod.row == "EVE" | mod.row == "VVE") {
      TempW2R <- array(0, dim = c(p, p, k))
      Tempz <- mclust::unmap(classy)

      for (j in 1:k) {
        deltaU.k[, , j] <- diag(eigen(sigmaU[, , j])$values, p, p)
        TempW2R[, , j] <- sigmaU[, , j] * sum(Tempz[, j])
      }

      gammaU <- eigen(rowSums(TempW2R, dims = 2) / ((det(rowSums(TempW2R, dims = 2)))^(1 / p)))$vectors
    } else {
      gammaU <- matrix(0, p, p)
    }
    if (mod.row == "EII" | mod.row == "VII") {
      for (j in 1:k) {
        sigmaU[, , j] <- diag(1, p, p)
      }
    }

    if (mod.col == "VE") {
      deltaV.k <- array(0, dim = c(r, r, k))
      TempW2C <- array(0, dim = c(r, r, k))
      Tempz <- mclust::unmap(classy)

      for (j in 1:k) {
        deltaV.k[, , j] <- diag(eigen(sigmaV[, , j])$values, r, r)
        TempW2C[, , j] <- sigmaV[, , j] * sum(Tempz[, j])
      }

      gammaV <- eigen(rowSums(TempW2C, dims = 2) / ((det(rowSums(TempW2C, dims = 2)))^(1 / r)))$vectors
    } else {
      deltaV.k <- array(0, dim = c(r, r, k))
      gammaV <- matrix(0, r, r)
    }
    if (mod.col == "II") {
      for (j in 1:k) {
        sigmaV[, , j] <- diag(1, r, r)
      }
    }

    ### Estimation ###

    tryCatch(while (check < 1) {
      m.iter <- m.iter + 1

      ### E - STEP ###

      if (density == "MVN") {
        for (j in 1:k) {
          dens[, j] <- dMVnorm(X = X, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j])
        }
      }

      if (density == "MVTIN") {
        for (j in 1:k) {
          dens[, j] <- dMVtin(X = X, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j], theta = theta[j])
        }

        for (j in 1:k) {
          delta <- sapply(1:num, function(i) tr(solve(sigmaU[, , j]) %*% (X[, , i] - M[, , j]) %*% solve(sigmaV[, , j]) %*% t(X[, , i] - M[, , j])))

          numer <- 2 * (zipfR::Igamma(a = ((p * r) / 2 + 2), x = (1 - theta[j]) * delta / 2, lower = FALSE) - zipfR::Igamma(a = ((p * r) / 2 + 2), x = delta / 2, lower = FALSE))
          den <- delta * (zipfR::Igamma(a = ((p * r) / 2 + 1), x = (1 - theta[j]) * delta / 2, lower = FALSE) - zipfR::Igamma(a = ((p * r) / 2 + 1), x = delta / 2, lower = FALSE))

          numer[numer < .Machine$double.xmin] <- .Machine$double.xmin
          den[den < .Machine$double.xmin] <- .Machine$double.xmin

          wtt <- numer / den
          wtt[wtt > 1] <- 0.999
          w[, j] <- wtt
        }
      }

      numerator <- matrix(rep(prior, num), num, k, byrow = TRUE) * dens
      mixt.dens <- rowSums(numerator)
      post <- numerator / mixt.dens

      ### M - STEP ###

      if (density == "MVN") {
        for (j in 1:k) {
          M[, , j] <- rowSums(X * post[, j][slice.index(X, 3)], dims = 2) / sum(post[, j])
        }
      } else {
        for (j in 1:k) {
          M[, , j] <- rowSums(X * (w[, j] * post[, j])[slice.index(X, 3)], dims = 2) / sum(w[, j] * post[, j])
        }
      }

      # ROWS COVARIANCE MATRIX

      if (density == "MVN") {
        for (j in 1:k) {
          pt1 <- aperm(tensor(sweep(X, 1:2, M[, , j]) * post[, j][slice.index(sweep(X, 1:2, M[, , j]), 3)], solve(sigmaV[, , j]), 2, 1), c(1, 3, 2))

          pt2 <- aperm(sweep(X, 1:2, M[, , j]), c(2, 1, 3))

          WR[, , j] <- tensor(pt1, pt2, c(2, 3), c(1, 3))
        } # Scatter Matrix
      } else {
        for (j in 1:k) {
          pt1 <- aperm(tensor(sweep(X, 1:2, M[, , j]) * (w[, j] * post[, j])[slice.index(sweep(X, 1:2, M[, , j]), 3)], solve(sigmaV[, , j]), 2, 1), c(1, 3, 2))

          pt2 <- aperm(sweep(X, 1:2, M[, , j]), c(2, 1, 3))

          WR[, , j] <- tensor(pt1, pt2, c(2, 3), c(1, 3))
        } # Scatter Matrix
      }

      if (mod.row == "EII") {
        phi <- tr(rowSums(WR, dims = 2)) / ((num) * p * r)

        for (j in 1:k) {
          sigmaU[, , j] <- phi * diag(1, p, p)
        }
      }

      if (mod.row == "VII") {
        for (j in 1:k) {
          phi.k[j] <- tr(WR[, , j]) / (p * r * sum(post[, j]))
          sigmaU[, , j] <- phi.k[j] * diag(1, p, p)
        }
      }

      if (mod.row == "EEI") {
        deltaU <- diag(diag(rowSums(WR, dims = 2)), p, p) / (det(diag(diag(rowSums(WR, dims = 2)), p, p)))^(1 / p)

        phi <- (det(diag(diag(rowSums(WR, dims = 2)), p, p)))^(1 / p) / ((num) * r)

        for (j in 1:k) {
          sigmaU[, , j] <- phi * deltaU
        }
      }

      if (mod.row == "VEI") {
        for (j in 1:k) {
          temp.numdeltaR[, , j] <- (1 / phi.k[j]) * WR[, , j]
        }

        deltaU <- diag(diag(rowSums(temp.numdeltaR, dims = 2)), p, p) / (det(diag(diag(rowSums(temp.numdeltaR, dims = 2)), p, p)))^(1 / p)

        for (j in 1:k) {
          phi.k[j] <- (tr(solve(deltaU) %*% WR[, , j])) / (p * r * sum(post[, j]))

          sigmaU[, , j] <- phi.k[j] * deltaU
        }
      }

      if (mod.row == "EVI") {
        for (j in 1:k) {
          deltaU.k[, , j] <- diag(diag(WR[, , j]), p, p) / (det(diag(diag(WR[, , j]), p, p)))^(1 / p)

          temp.phi2[j] <- det(diag(diag(WR[, , j]), p, p))^(1 / p)
        }

        phi <- sum(temp.phi2) / ((num) * r)

        for (j in 1:k) {
          sigmaU[, , j] <- phi * deltaU.k[, , j]
        }
      }

      if (mod.row == "VVI") {
        for (j in 1:k) {
          deltaU.k[, , j] <- diag(diag(WR[, , j]), p, p) / (det(diag(diag(WR[, , j]), p, p)))^(1 / p)

          phi.k[j] <- det(diag(diag(WR[, , j]), p, p))^(1 / p) / (r * sum(post[, j]))

          sigmaU[, , j] <- phi.k[j] * deltaU.k[, , j]
        }
      }

      if (mod.row == "EEE") {
        for (j in 1:k) {
          sigmaU[, , j] <- rowSums(WR, dims = 2) / ((num) * r)
        }
      }

      if (mod.row == "VEE") {
        for (j in 1:k) {
          temp.numdeltaR[, , j] <- (1 / phi.k[j]) * WR[, , j]
        }

        deltaU <- rowSums(temp.numdeltaR, dims = 2) / ((det(rowSums(temp.numdeltaR, dims = 2)))^(1 / p))

        for (j in 1:k) {
          phi.k[j] <- tr(solve(deltaU) %*% WR[, , j]) / (p * r * sum(post[, j]))

          sigmaU[, , j] <- phi.k[j] * deltaU
        }
      }

      if (mod.row == "EVE") {
        while (check2 < 1) {
          m.iter2 <- m.iter2 + 1

          for (j in 1:k) {
            ftemp.r[, , j] <- tcrossprod(solve(deltaU.k[, , j]), gammaU) %*% WR[, , j] - max(eigen(WR[, , j])$values) * tcrossprod(solve(deltaU.k[, , j]), gammaU)
          }

          f <- rowSums(ftemp.r, dims = 2)

          MM.r.new <- tr(f %*% gammaU)

          if ((abs(MM.r.new - MM.r.old)) < tol2 | m.iter2 == maxit2) {
            check2 <- 1
            res.svd <- svd(f)
            gammaU <- tcrossprod(res.svd$v, res.svd$u)
          } else {
            res.svd <- svd(f)
            gammaU <- tcrossprod(res.svd$v, res.svd$u)
          }

          MM.r.old <- MM.r.new
        }

        m.iter2 <- 0
        check2 <- 0
        MM.r.old <- -Inf

        for (j in 1:k) {
          deltaU.k[, , j] <- diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p) / (det(diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p)))^(1 / p)

          numphi[j] <- tr(gammaU %*% tcrossprod(solve(deltaU.k[, , j]), gammaU) %*% WR[, , j])
        }

        phi <- sum(numphi) / ((num) * p * r)

        for (j in 1:k) {
          sigmaU[, , j] <- phi * gammaU %*% tcrossprod(deltaU.k[, , j], gammaU)
        }
      }

      if (mod.row == "VVE") {
        while (check2 < 1) {
          m.iter2 <- m.iter2 + 1

          for (j in 1:k) {
            ftemp.r[, , j] <- tcrossprod(solve(deltaU.k[, , j]), gammaU) %*% WR[, , j] - max(eigen(WR[, , j])$values) * tcrossprod(solve(deltaU.k[, , j]), gammaU)
          }

          f <- rowSums(ftemp.r, dims = 2)

          MM.r.new <- tr(f %*% gammaU)

          if ((abs(MM.r.new - MM.r.old)) < tol2 | m.iter2 == maxit2) {
            check2 <- 1
            res.svd <- svd(f)
            gammaU <- tcrossprod(res.svd$v, res.svd$u)
          } else {
            res.svd <- svd(f)
            gammaU <- tcrossprod(res.svd$v, res.svd$u)
          }

          MM.r.old <- MM.r.new
        }

        m.iter2 <- 0
        check2 <- 0
        MM.r.old <- -Inf

        for (j in 1:k) {
          deltaU.k[, , j] <- diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p) / (det(diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p)))^(1 / p)
          phi.k[j] <- (det(diag(diag(crossprod(gammaU, WR[, , j]) %*% gammaU), p, p))^(1 / p)) / (r * sum(post[, j]))
          sigmaU[, , j] <- phi.k[j] * gammaU %*% tcrossprod(deltaU.k[, , j], gammaU)
        }
      }

      if (mod.row == "EEV") {
        for (j in 1:k) {
          tempW_EEV[[j]] <- eigen(WR[, , j])

          gammaU.k[, , j] <- tempW_EEV[[j]][["vectors"]]

          tempomegaR[, , j] <- diag(tempW_EEV[[j]][["values"]], p, p)
        }

        deltaU <- rowSums(tempomegaR, dims = 2) / ((det(rowSums(tempomegaR, dims = 2)))^(1 / p))

        phi <- ((det(rowSums(tempomegaR, dims = 2)))^(1 / p)) / ((num) * r)

        for (j in 1:k) {
          sigmaU[, , j] <- phi * gammaU.k[, , j] %*% tcrossprod(deltaU, gammaU.k[, , j])
        }
      }

      if (mod.row == "VEV") {
        for (j in 1:k) {
          tempW_EEV[[j]] <- eigen(WR[, , j])

          gammaU.k[, , j] <- tempW_EEV[[j]][["vectors"]]

          tempomegaR[, , j] <- diag(tempW_EEV[[j]][["values"]], p, p)

          temp.numdeltaR[, , j] <- (1 / phi.k[j]) * tempomegaR[, , j]
        }

        deltaU <- rowSums(temp.numdeltaR, dims = 2) / ((det(rowSums(temp.numdeltaR, dims = 2)))^(1 / p))

        for (j in 1:k) {
          phi.k[j] <- tr(tempomegaR[, , j] %*% solve(deltaU)) / (p * r * sum(post[, j]))

          sigmaU[, , j] <- phi.k[j] * gammaU.k[, , j] %*% tcrossprod(deltaU, gammaU.k[, , j])
        }
      }

      if (mod.row == "EVV") {
        for (j in 1:k) {
          V_EVV.U.K[, , j] <- WR[, , j] / ((det(WR[, , j]))^(1 / p))

          numphi[j] <- det(WR[, , j])^(1 / p)
        }

        phi <- sum(numphi) / ((num) * r)

        for (j in 1:k) {
          sigmaU[, , j] <- phi * V_EVV.U.K[, , j]
        }
      }

      if (mod.row == "VVV") {
        for (j in 1:k) {
          sigmaU[, , j] <- WR[, , j] / (r * sum(post[, j]))
        }
      }

      # COLUMNS COVARIANCE MATRIX

      if (density == "MVN") {
        for (j in 1:k) {
          pt1 <- aperm(tensor(aperm(sweep(X, 1:2, M[, , j]), c(2, 1, 3)) * post[, j][slice.index(aperm(sweep(X, 1:2, M[, , j]), c(2, 1, 3)), 3)], solve(sigmaU[, , j]), 2, 1), c(1, 3, 2))

          pt2 <- sweep(X, 1:2, M[, , j])

          WC[, , j] <- tensor(pt1, pt2, c(2, 3), c(1, 3))
        } # Scatter Matrix
      } else {
        for (j in 1:k) {
          pt1 <- aperm(tensor(aperm(sweep(X, 1:2, M[, , j]), c(2, 1, 3)) * (w[, j] * post[, j])[slice.index(aperm(sweep(X, 1:2, M[, , j]), c(2, 1, 3)), 3)], solve(sigmaU[, , j]), 2, 1), c(1, 3, 2))

          pt2 <- sweep(X, 1:2, M[, , j])

          WC[, , j] <- tensor(pt1, pt2, c(2, 3), c(1, 3))
        } # Scatter Matrix
      }

      if (mod.col == "II") {
        for (j in 1:k) {
          sigmaV[, , j] <- diag(1, r, r)
        }
      }

      if (mod.col == "EI") {
        deltaV <- diag(diag(rowSums(WC, dims = 2)), r, r) / (det(diag(diag(rowSums(WC, dims = 2)), r, r)))^(1 / r)

        for (j in 1:k) {
          sigmaV[, , j] <- deltaV
        }
      }

      if (mod.col == "VI") {
        for (j in 1:k) {
          sigmaV[, , j] <- diag(diag(WC[, , j]), r, r) / (det(diag(diag(WC[, , j]), r, r)))^(1 / r)
        }
      }

      if (mod.col == "EE") {
        for (j in 1:k) {
          sigmaV[, , j] <- rowSums(WC, dims = 2) / ((det(rowSums(WC, dims = 2)))^(1 / r))
        }
      }

      if (mod.col == "VE") {
        while (check2 < 1) {
          m.iter2 <- m.iter2 + 1

          for (j in 1:k) {
            ftemp.c[, , j] <- tcrossprod(solve(deltaV.k[, , j]), gammaV) %*% WC[, , j] - max(eigen(WC[, , j])$values) * tcrossprod(solve(deltaV.k[, , j]), gammaV)
          }

          f.C <- rowSums(ftemp.c, dims = 2)

          MM.c.new <- tr(f.C %*% gammaV)

          if ((abs(MM.c.new - MM.c.old)) < tol2 | m.iter2 == maxit2) {
            check2 <- 1
            res.svd.C <- svd(f.C)
            gammaV <- tcrossprod(res.svd.C$v, res.svd.C$u)
          } else {
            res.svd.C <- svd(f.C)
            gammaV <- tcrossprod(res.svd.C$v, res.svd.C$u)
          }

          MM.c.old <- MM.c.new
        }

        m.iter2 <- 0
        check2 <- 0
        MM.c.old <- -Inf

        for (j in 1:k) {
          deltaV.k[, , j] <- diag(diag(crossprod(gammaV, WC[, , j]) %*% gammaV), r, r) / (det(diag(diag(crossprod(gammaV, WC[, , j]) %*% gammaV), r, r)))^(1 / r)
        }

        for (j in 1:k) {
          sigmaV[, , j] <- gammaV %*% tcrossprod(deltaV.k[, , j], gammaV)
        }
      }

      if (mod.col == "EV") {
        for (j in 1:k) {
          tempW_EV[[j]] <- eigen(WC[, , j])

          gammaV.k[, , j] <- tempW_EV[[j]][["vectors"]]

          tempomegaC[, , j] <- diag(tempW_EV[[j]][["values"]], r, r)
        }

        deltaV <- rowSums(tempomegaC, dims = 2) / ((det(rowSums(tempomegaC, dims = 2)))^(1 / r))

        for (j in 1:k) {
          sigmaV[, , j] <- gammaV.k[, , j] %*% tcrossprod(deltaV, gammaV.k[, , j])
        }
      }

      if (mod.col == "VV") {
        for (j in 1:k) {
          sigmaV[, , j] <- WC[, , j] / ((det(WC[, , j]))^(1 / r))
        }
      }

      if (k == 1) {
        prior <- 1
      } else {
        prior <- colMeans(post)
      }

      if (density == "MVTIN") {
        if (mod.theta == "E") {
          theta <- rep(Mstep_AECM(X = X, k = k, weights = post, M = M, sigmaU = sigmaU, sigmaV = sigmaV, mod.theta = mod.theta), k)
        } else {
          for (j in 1:k) {
            theta[j] <- Mstep_AECM(X = X, weights = post[, j], M = M[, , j], sigmaU = sigmaU[, , j], sigmaV = sigmaV[, , j], mod.theta = mod.theta)
          }
        }

        for (j in 1:k) {
          dens[, j] <- dMVtin(X = X, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j], theta = theta[j])
        }
      }

      if (density == "MVN") {
        for (j in 1:k) {
          dens[, j] <- dMVnorm(X = X, M = M[, , j], U = sigmaU[, , j], V = sigmaV[, , j])
        }
      }

      # mixture density

      numerator <- matrix(rep(prior, num), num, k, byrow = TRUE) * dens
      mixt.dens <- rowSums(numerator)
      loglik.new <- sum(log(mixt.dens))
      ll <- c(ll, loglik.new)

      # stopping rule

      if ((floor(100 * loglik.new) / 100 - floor(100 * loglik.old) / 100) < tol) {
        check <- 1
      }

      if (floor(100 * loglik.new) / 100 < floor(100 * loglik.old) / 100) {
        mark <- 1
      } else {
        mark <- 0
      }

      if (m.iter == maxit) {
        check <- 1
      }

      loglik.old <- loglik.new
    }, error = function(e) {
      cnt <<- 1
    })

    #### Output ####

    if (cnt == 0 & mark == 0) {
      plot(ll, type = "l", xlab = "Iterations", ylab = "Log-Likelihoods")

      # --------------------- #
      #     Classification    #
      # --------------------- #

      if (k == 1) {
        classification <- rep(1, num)
      } else {
        colnames(post) <- c(1:k)
        classification <- as.numeric(colnames(post)[max.col(post, ties.method = "first")])
      }

      # -------------------- #
      # Information criteria #
      # -------------------- #

      # Number of parameters

      meanpar <- (p * r) * k

      if (mod.row == "EII") {
        rowpar <- 1
      }
      if (mod.row == "VII") {
        rowpar <- k
      }
      if (mod.row == "EEI") {
        rowpar <- p
      }
      if (mod.row == "VEI") {
        rowpar <- k + (p - 1)
      }
      if (mod.row == "EVI") {
        rowpar <- 1 + k * (p - 1)
      }
      if (mod.row == "VVI") {
        rowpar <- k * p
      }
      if (mod.row == "EEE") {
        rowpar <- p * (p + 1) / 2
      }
      if (mod.row == "VEE") {
        rowpar <- k - 1 + p * (p + 1) / 2
      }
      if (mod.row == "EVE") {
        rowpar <- 1 + k * (p - 1) + p * (p - 1) / 2
      }
      if (mod.row == "VVE") {
        rowpar <- k * p + p * (p - 1) / 2
      }
      if (mod.row == "EEV") {
        rowpar <- p + k * p * (p - 1) / 2
      }
      if (mod.row == "VEV") {
        rowpar <- k + (p - 1) + (k * p * (p - 1) / 2)
      }
      if (mod.row == "EVV") {
        rowpar <- 1 + k * (p * ((p + 1) / 2) - 1)
      }
      if (mod.row == "VVV") {
        rowpar <- k * p * (p + 1) / 2
      }

      if (mod.col == "II") {
        colpar <- 0
      }
      if (mod.col == "EI") {
        colpar <- r - 1
      }
      if (mod.col == "VI") {
        colpar <- k * (r - 1)
      }
      if (mod.col == "EE") {
        colpar <- r * ((r + 1) / 2) - 1
      }
      if (mod.col == "VE") {
        colpar <- k * (r - 1) + r * (r - 1) / 2
      }
      if (mod.col == "EV") {
        colpar <- (r - 1) + k * r * (r - 1) / 2
      }
      if (mod.col == "VV") {
        colpar <- k * (r * ((r + 1) / 2) - 1)
      }

      weights <- k - 1

      npar <- meanpar + rowpar + colpar + weights

      if (density == "MVTIN") {
        if (mod.theta == "V") {
          npar <- npar + k
        } else {
          npar <- npar + 1
        }
      }

      name <- c(mod.row, mod.col, mod.theta)

      # to be minimized

      BIC <- -2 * loglik.new + npar * log(num)

      ptm2 <- proc.time() - ptm
      time <- ptm2[3]

      return(list(
        den = density, name = name, prior = prior, M = M, sigmaU = sigmaU, sigmaV = sigmaV, theta = theta,
        loglik = loglik.new, mark = mark, check = check, npar = npar, iter = m.iter, time = time, BIC = BIC, group = classification
      ))
    } else {
      return(NA)
    }
  }

  for (g in 1:length(k)) {
    if (density == "MVN") {
      print(paste(paste(paste(paste("Fitting parsimonious matrix", density), paste("mixtures with k =", k[g])), paste("and", mod.row)), paste("-", mod.col), paste("parsimonious structure")))
    } else {
      print(paste(paste(paste(paste(paste("Fitting parsimonious matrix", density), paste("mixtures with k =", k[g])), paste("and", mod.row)), paste("-", mod.col), paste("-", mod.theta)), paste("parsimonious structure")))
    }

    cluster <- makeCluster(nThreads, type = "SOCK")
    registerDoSNOW(cluster)

    pb <- txtProgressBar(max = pt.mod, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    oper[[g]] <- foreach(l = 1:pt.mod, .combine = "comb", .packages = c("tensor", "zipfR", "expint"), .multicombine = TRUE, .init = list(list()), .options.snow = opts) %dopar% {
      res <- tryCatch(Pars_mixt(X = X, k = k[g], init.par = init.par[[g]], mod.row = as.character(list.mod[l, 1]), mod.col = as.character(list.mod[l, 2]), mod.theta = as.character(list.mod[l, 3]), tol = tol, tol2 = tol2, maxit = maxit, maxit2 = maxit2, density = density), error = function(e) {
        NA
      })

      list(res)
    }

    stopCluster(cluster)
    registerDoSEQ()
    close(pb)
  }

  return(oper)
}

r_data.gen <- function(num, p, r, k, pi, M, U, V, theta, density) {
  X <- array(0, dim = c(p, r, num))
  temp <- sample(1:k, size = num, replace = TRUE, prob = pi)

  if (density == "MVN") {
    for (i in 1:num) {
      X[, , i] <- rmatrixnorm(M = M[, , temp[i]], U = U[, , temp[i]], V = V[, , temp[i]])
    }
  }
  if (density == "MVTIN") {
    for (i in 1:num) {
      w <- stats::runif(n = 1, min = 1 - theta[temp[i]], 1)

      X[, , i] <- rmatrixnorm(M = M[, , temp[i]], U = U[, , temp[i]] / w, V = V[, , temp[i]])
    }
  }

  return(list(X = X, trueclass = temp))
}

extract.bestM <- function(results, mod.row = "all", mod.col = "all", mod.theta = "all", density) {
  k <- length(results)
  num.mod <- length(results[[1]][[1]])
  if (all(mod.row == "all")) {
    model.row <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV")
  } else {
    model.row <- mod.row
  }
  if (all(mod.col == "all")) {
    model.col <- c("II", "EI", "VI", "EE", "VE", "EV", "VV")
  } else {
    model.col <- mod.col
  }

  if (density == "MVN") {
    list.mod <- expand.grid(model.row, model.col)
  } else {
    if (all(mod.theta == "all")) {
      model.theta <- c("E", "V")
    } else {
      model.theta <- mod.theta
    }
    list.mod <- expand.grid(model.row, model.col, model.theta)
  }

  list.mod2 <- do.call("rbind", replicate(k, list.mod, simplify = FALSE))
  count.k <- sort(rep(1:k, num.mod))
  count.mod <- rep(1:num.mod, k)
  list.mod3 <- data.frame(list.mod2, count.k, count.mod)

  allBIC <- numeric(k * num.mod)
  cont <- 0

  for (j in 1:k) {
    for (i in 1:num.mod) {
      if (!all(is.na(results[[j]][[1]][[i]]))) {
        if (min(results[[j]][[1]][[i]][["prior"]]) > 0.05) {
          cont <- cont + 1
          allBIC[cont] <- results[[j]][[1]][[i]][["BIC"]]
        } else {
          cont <- cont + 1
          allBIC[cont] <- NA
        }
      } else {
        cont <- cont + 1
        allBIC[cont] <- NA
      }
    }
  }

  allBIC[allBIC == -Inf] <- Inf

  winBIC <- which.min(allBIC)

  tempBIC <- list.mod3[winBIC, ]

  if (density == "MVN") {
    bestBIC <- results[[as.numeric(tempBIC[3])]][[1]][[as.numeric(tempBIC[4])]]
  } else {
    bestBIC <- results[[as.numeric(tempBIC[4])]][[1]][[as.numeric(tempBIC[5])]]
  }

  return(list(bestBIC = bestBIC))
}
