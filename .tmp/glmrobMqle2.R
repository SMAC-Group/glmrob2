#### "Inconsistent" Mallows quasi-likelihood estimator of E. Cantoni and E. Ronchetti (2001)
#### based on M. Maechler's code in R package "robustbase"
#### https://github.com/cran/robustbase/blob/master/R/glmrobMqle.R

glmrobMqle2 <-
  function(X, y, weights = NULL, start = NULL, offset = NULL,
           family, weights.on.x = "none",
           control = glmrobMqle.control(), intercept = TRUE,
           trace = FALSE)
  {
    ## To DO:
    ## o weights are not really implemented as *extra* user weights; rather as "glm-weights"
    ## o offset is not fully implemented (really? -- should have test case!)

    if(!is.matrix(X)) X <- as.matrix(X)
    ## never used:
    ##     xnames <- dimnames(X)[[2]]
    ##     ynames <- if (is.matrix(y)) rownames(y) else names(y)
    nobs <- NROW(y)
    stopifnot(nobs == nrow(X))
    if (is.null(weights))
      weights <- rep.int(1, nobs)
    else if(any(weights <= 0))
      stop("All weights must be positive")
    if (is.null(offset))
      offset <- rep.int(0, nobs) else if(!all(offset==0))
        warning("'offset' not fully implemented")
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
      stop("illegal 'family' argument")
    mu.eta <- family$mu.eta
    if (is.null(valideta <- family$valideta)) valideta <- function(eta) TRUE
    if (is.null(validmu	 <- family$validmu))  validmu <-  function(mu) TRUE

    ncoef <- ncol(X)
    w.x <- robXweights(weights.on.x, X=X, intercept=intercept)

    ### Initializations
    stopifnot(control$maxit >= 1, (tcc <- control$tcc) >= 0)

    ## note that etastart and mustart are used to make 'family$initialize' run
    etastart <- NULL;  mustart <- NULL
    ## note that 'weights' are used and set by binomial()$initialize !
    eval(family$initialize) ## --> n, mustart, y and weights (=ni)
    ni <- as.vector(weights)# dropping attributes for computation
    ##
    if(is.null(start))
      start <- glm.fit(x = X, y = y, weights = weights, offset = offset,
                       family = family)$coefficients
    if(any(ina <- is.na(start))) {
      cat("initial start 'theta' has NA's; eliminating columns X[, j];",
          "j = ", pasteK(which(ina)),"\n")
      theta.na <- start
      X <- X[, !ina, drop = FALSE]
      start <- glm.fit(x = X, y = y, weights = weights, offset = offset,
                       family = family)$coefficients
      if(any(is.na(start)))
        stop("start 'theta' has still NA's .. badly singular x\n")
      ## FIXME
      ncoef <- length(start)
    }

    thetaOld <- theta <- as.vector(start) # as.v*(): dropping attributes
    eta <- as.vector(X %*% theta)
    mu <- linkinv(eta) # mu estimates pi (in [0,1]) at the binomial model
    if (!(validmu(mu) && valideta(eta)))
      stop("Cannot find valid starting values: You need help")
    ##
    switch(family$family,
           "binomial" = {
             Epsi.init <- EpsiBin.init
             Epsi <- EpsiBin
             EpsiS <- EpsiSBin
             Epsi2 <- Epsi2Bin
             phiEst <- phiEst.cl <- 1
           },
           "poisson" = {
             Epsi.init <- EpsiPois.init
             Epsi <- EpsiPois
             EpsiS <- EpsiSPois
             Epsi2 <- Epsi2Pois
             phiEst <- phiEst.cl <- expression({1})
           },
           "gaussian" = {
             Epsi.init <- EpsiGaussian.init
             Epsi <- EpsiGaussian
             EpsiS <- EpsiSGaussian
             Epsi2 <- Epsi2Gaussian
             phiEst.cl <- phiGaussianEst.cl
             phiEst <- phiGaussianEst
           },
           "Gamma" = { ## added by ARu
             Epsi.init <- EpsiGamma.init
             Epsi <- EpsiGamma
             EpsiS <- EpsiSGamma
             Epsi2 <- Epsi2Gamma
             phiEst.cl <- phiGammaEst.cl
             phiEst <- phiGammaEst
           },
           ## else
           stop(gettextf("family '%s' not yet implemented", family$family),
                domain=NA)
    )

    sV <- NULL # FIXME workaround for codetools

    comp.V.resid <- expression({
      Vmu <- variance(mu)
      if (any(is.na(Vmu)))  stop("NAs in V(mu)")
      if (any(Vmu == 0))    stop("0s in V(mu)")
      sVF <- sqrt(Vmu)   # square root of variance function
      residP <- (y - mu)* sni/sVF  # Pearson residuals
    })

    comp.scaling <- expression({
      sV <- sVF * sqrt(phi)
      residPS <- residP/sqrt(phi) # scaled Pearson residuals
    })

    comp.Epsi.init <- expression({
      ## d mu / d eta :
      dmu.deta <- mu.eta(eta)
      if (any(is.na(dmu.deta))) stop("NAs in d(mu)/d(eta)")
      ## "Epsi init" :
      H <- floor(mu*ni - tcc* sni*sV)
      K <- floor(mu*ni + tcc* sni*sV)
      eval(Epsi.init)
    })


    ### Iterations

    if(trace && ncoef) {
      cat("Initial theta: \n")
      local({names(theta) <- names(start); print(theta) })

      digits <- max(1, getOption("digits") - 5)
      w.th.1 <- 6+digits # width of one number; need 8 for 2 digits: "-4.8e-11"
      width.th <- ncoef*(w.th.1 + 1) - 1
      cat(sprintf("%3s | %*s | %12s\n",
                  "it", width.th, "d{theta}", "rel.change"))
      mFormat <- function(x, wid) {
        r <- formatC(x, digits=digits, width=wid)
        sprintf("%*s", wid, sub("e([+-])0","e\\1", r))
      }
    }

    sni <- sqrt(ni)
    eval(comp.V.resid) #-> (Vmu, sVF, residP)
    phi <- eval(phiEst.cl)
    ## Determine the range of phi values based on the distribution of |residP|
    Rphi <- c(1e-12, 3*median(abs(residP)))^2
    conv <- FALSE
    if(ncoef) for (nit in 1:control$maxit) {
      eval(comp.scaling) #-> (sV, residPS)
      eval(comp.Epsi.init)
      ## Computation of alpha and (7) using matrix column means:

      # ------ Modif here --------
      # previous code:
      # cpsi <- pmax.int(-tcc, pmin.int(residPS,tcc)) - eval(Epsi)
      # new code:
      # cpsi <-
      # ------ End modif --------

      EEq <- colMeans(cpsi * w.x * sni/sV * dmu.deta * X)
      ##
      ## Solve  1/n (t(X) %*% B %*% X) %*% delta.coef	  = EEq
      DiagB <- eval(EpsiS) /(sni*sV) * w.x * (ni*dmu.deta)^2
      if(any(n0 <- ni == 0)) DiagB[n0] <- 0 # instead of NaN
      Dtheta <- solve(crossprod(X, DiagB*X)/nobs, EEq)
      if (any(!is.finite(Dtheta))) {
        warning("Non-finite coefficients at iteration ", nit)
        break
      }
      theta <- thetaOld + Dtheta
      eta <- as.vector(X %*% theta) + offset
      mu <- linkinv(eta)

      ## estimation of the dispersion parameter
      eval(comp.V.resid)
      phi <- eval(phiEst)

      ## Check convergence: relative error < tolerance
      relE <- sqrt(sum(Dtheta^2)/max(1e-20, sum(thetaOld^2)))
      conv <- relE <= control$acc
      if(trace) {
        cat(sprintf("%3d | %*s | %12g\n", nit, width.th,
                    paste(mFormat(Dtheta, w.th.1),
                          collapse=" "), relE))
      }
      if(conv)
        break
      thetaOld <- theta
    } ## end of iteration
    else { ## ncoef == 0
      conv <- TRUE
      nit <- 0
    }
    if (!conv)
      warning("Algorithm did not converge")

    eps <- 10 * .Machine$double.eps
    switch(family$family,
           "binomial" = {
             if (any(mu/weights > 1 - eps) || any(mu/weights < eps))
               warning("fitted probabilities numerically 0 or 1 occurred")
           },
           "poisson" = {
             if (any(mu < eps))
               warning("fitted rates numerically 0 occurred")
           })

    eval(comp.V.resid) #-> (Vmu, sVF, residP)
    eval(comp.scaling) #-> (sV, residPS)

    ## Estimated asymptotic covariance of the robust estimator
    if(ncoef) {
      eval(comp.Epsi.init)
      alpha <- colMeans(eval(Epsi) * w.x * sni/sV * dmu.deta * X)
      DiagA <- eval(Epsi2) / (ni*sV^2)* w.x^2* (ni*dmu.deta)^2
      matQ  <- crossprod(X, DiagA*X)/nobs - tcrossprod(alpha, alpha)

      DiagB <- eval(EpsiS) / (sni*sV)* w.x * (ni*dmu.deta)^2
      if(any(n0 <- ni == 0)) DiagB[n0] <- 0 # instead of NaN
      matM <- crossprod(X, DiagB*X)/nobs
      matMinv <- solve(matM)
      asCov <-  matMinv %*% matQ %*% matMinv / nobs
    } else { ## ncoef == 0
      matM <- matQ <- asCov <- matrix(NA_real_, 0,0)
    }

    if(any(ina)) {# put NA's back, extending theta[] to "original length"
      ok <- !ina
      theta.na[ok] <- theta ; theta <- theta.na
      ## also extend the "p x p" matrices with NA's --
      ##No : lm() and glm() also do *not* do this
      ##No  p <- length(theta)
      ##No  nm <- names(theta)
      ##No  M <- matrix(NA_real_, p, p, dimnames = list(nm,nm))
      ##No  Mn <- M; Mn[ok, ok] <- asCov ; asCov <- Mn
      ##No  Mn <- M; Mn[ok, ok] <- matM  ; matM  <- Mn
      ##No  Mn <- M; Mn[ok, ok] <- matQ  ; matQ  <- Mn
    }

    w.r <- pmin(1, tcc/abs(residPS))
    names(mu) <- names(eta) <- names(residPS) # re-add after computation
    list(coefficients = theta, residuals = residP, # s.resid = residPS,
         fitted.values = mu,
         w.r = w.r, w.x = w.x, ni = ni, dispersion = phi, cov = asCov,
         matM = matM, matQ = matQ, tcc = tcc, family = family,
         linear.predictors = eta, deviance = NULL, iter = nit, y = y,
         converged = conv)
  }
