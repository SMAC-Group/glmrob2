#### The code is largely inspired, even sometimes copied, from:
#### o glmrob/glmrobMqle functions of R package "robustbase";
#### o glm/glm.fit functions of R package "stats".

#' @title Robust Fitting of Generalized Linear Models (under development)
#' \code{glmrob2} is used to fit generalized linear models.
#' The iterative re-weighted least squared is used for fitting, similarly
#' to \code{\link[stats]{glm}}. Currently only \code{\link[stats]{binomial}} with \code{logit}
#' link is supported. Be aware estimators are not Fisher consistent.
#' @param formula an object of class formula (see \code{\link[stats]{formula}}).
#' @param family see
#' @param data see
#' @param weights see
#' @param subset see
#' @param na.action see
#' @param start see
#' @param etastart see
#' @param mustart see
#' @param offset see
#' @param control see
#' @param model see
#' @param x see
#' @param y see
#' @param contrasts see
#' @param ... see
#' @importFrom Rdpack reprompt
#' @importFrom stats .getXlevels family gaussian binomial Gamma is.empty.model model.extract model.matrix model.offset model.response model.weights
#' @importFrom utils getFromNamespace
#' @importFrom robustbase glmrob
#' @author Samuel Orso
#' @example /inst/examples/eg1.R
#' @export
glmrob2 <- function(formula, family = gaussian, data, weights,
                subset, na.action, start = NULL,
                etastart, mustart, offset,
                control = list(...),
                model = TRUE, #method = "glm.fit",
                x = FALSE, y = TRUE,
                contrasts = NULL, ...)
{
  call <- match.call()
  ## family
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  ## extract x, y, etc from the model formula and frame
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  # if(identical(method, "model.frame")) return(mf)

  # if (!is.character(method) && !is.function(method))
  #   stop("invalid 'method' argument")
  ## for back-compatibility in return result
  # if (identical(method, "glm.fit"))
  control <- do.call("glmrob2.control", control)

  mt <- attr(mf, "terms") # allow model.frame to have updated it

  Y <- model.response(mf, "any") # e.g. factors are allowed
  ## avoid problems with 1D arrays, but keep names
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(,NROW(Y), 0L)
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  weights <- as.vector(model.weights(mf))
  if(!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  ## check weights and offset
  if( !is.null(weights) && any(weights < 0) )
    stop("negative weights not allowed")

  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
  }
  ## these allow starting values to be expressed in terms of other vars.
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")

  ## We want to set the name on this call and the one below for the
  ## sake of messages from the fitter function
  fit <- eval(call("glmrob2.fit", #if(is.function(method)) "method" else method,
                   x = X, y = Y, weights = weights, start = start,
                   etastart = etastart, mustart = mustart,
                   offset = offset, family = family, control = control,
                   intercept = attr(mt, "intercept") > 0L))

  ## This calculated the null deviance from the intercept-only model
  ## if there is one, otherwise from the offset-only model.
  ## We need to recalculate by a proper fit if there is intercept and
  ## offset.
  ##
  ## The glm.fit calculation could be wrong if the link depends on the
  ## observations, so we allow the null deviance to be forced to be
  ## re-calculated by setting an offset (provided there is an intercept).
  ## Prior to 2.4.0 this was only done for non-zero offsets.
  if(length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <-
      eval(call("glmrob2.fit", #if(is.function(method)) "method" else method,
                x = X[, "(Intercept)", drop=FALSE], y = Y,
                weights = weights, offset = offset, family = family,
                control = control, intercept = TRUE))
    ## That fit might not have converged ....
    if(!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }
  if(model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if(x) fit$x <- X
  if(!y) fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula,
                     terms = mt, data = data,
                     offset = offset, control = control, method = "glmrob2.fit",
                     contrasts = attr(X, "contrasts"),
                     xlevels = .getXlevels(mt, mf)))
  class(fit) <- c(fit$class, c("glmrob2", "glm"))
  fit
}

#' @rdname glmrob2
#' @param intercept see
#' @export
glmrob2.fit <- function (x, y, weights = rep(1, nobs), start = NULL,
                         etastart = NULL, mustart = NULL, offset = rep(0, nobs),
                         family = gaussian(), control = list(), intercept = TRUE) {
  control <- do.call("glmrob2.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if(is.matrix(y)) rownames(y) else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  ## define weights and offset if needed
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)

  ## get family functions:
  # TODO: for the moment only poisson and binomial families are supported (due to phi estimation)
  if(all(family$family != "binomial", family$family != "poisson"))
    stop(gettextf("family '%s' not yet implemented", family$family),
       domain=NA)

  variance <- family$variance
  linkinv  <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv) )
    stop("'family' argument seems not to be a valid family object", call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if(is.null(x)) if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu  <- unless.null(family$validmu,  function(mu) TRUE)
  if(is.null(mustart)) {
    ## calculates mustart and may change y and weights and set n (!)
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  #### we need two additional functions:
  #### o second derivative of link function and first
  #### o first derivative of variance
  #### FIXME: consider creating a "family"
  switch(family$family,
         "binomial" = {
           if(family$link!="logit") stop("Only implemented for 'logit' link")
           g2 <- function(mu) (2 * mu - 1) / mu^2 / (1 - mu)^2
           variancep <- function(mu) 1 - 2 * mu
         },
         ## else
         stop(gettextf("family '%s' not yet implemented", family$family),
              domain=NA)
  )

  ## get robust functions (from robustbase):
  .psi.conv.cc <- getFromNamespace(".psi.conv.cc", ns = "robustbase")
  .psi2ipsi <- getFromNamespace(".psi2ipsi", ns = "robustbase")
  C_psifun <- getNativeSymbolInfo("R_psifun", PACKAGE = getLoadedDLLs()$robustbase)
  c.psi <- .psi.conv.cc(control$psi, control$tuning.psi) # control tuning constant
  ipsi <- .psi2ipsi(control$psi) # control name and get its "number"

  ## additional functions (from stats)
  C_Cdqrls <- getNativeSymbolInfo("Cdqrls", PACKAGE = getLoadedDLLs()$stats)

  ##
  if(EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
      stop("invalid linear predictor values in empty model", call. = FALSE)
    mu <- linkinv(eta)
    ## calculate initial deviance and coefficient
    if (!validmu(mu))
      stop("invalid fitted means in empty model", call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep_len(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  } else {
    coefold <- NULL
    eta <-
      if(!is.null(etastart)) etastart
    else if(!is.null(start))
      if (length(start) != nvars)
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", nvars, paste(deparse(xnames), collapse=", ")),
             domain = NA)
    else {
      coefold <- start
      offset + as.vector(if (NCOL(x) == 1L) x * start else x %*% start)
    }
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
      stop("cannot find valid starting values: please specify some", call. = FALSE)
    ## calculate initial deviance and coefficient
    devold <- sum(dev.resids(y, mu, weights))
    boundary <- conv <- FALSE

    ## from glmrobMqle
    if(control$trace && nvars) {
      # cat("Initial coefficients: \n")
      # local({names(coef) <- names(start); print(coef) })

      digits <- max(1, getOption("digits") - 5)
      w.th.1 <- 6+digits # width of one number; need 8 for 2 digits: "-4.8e-11"
      width.th <- nvars*(w.th.1 + 1) - 1
      cat(sprintf("%3s | %*s | %12s\n",
                  "it", width.th, "d{theta}", "rel.change"))
      mFormat <- function(x, wid) {
        rr <- formatC(x, digits=digits, width=wid)
        sprintf("%*s", wid, sub("e([+-])0","e\\1", rr))
      }
    }

    ##------------- THE Iteratively Reweighting L.S. iteration -----------
    for (iter in 1L:control$maxit) {
      good <- weights > 0
      varmu <- variance(mu)[good]
      if (anyNA(varmu))
        stop("NAs in V(mu)")
      if (any(varmu == 0))
        stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta)
      if (any(is.na(mu.eta.val[good])))
        stop("NAs in d(mu)/d(eta)")
      ## drop observations for which w will be zero
      good <- (weights > 0) & (mu.eta.val != 0)

      if (all(!good)) {
        conv <- FALSE
        warning(gettextf("no observations informative at iteration %d",
                         iter), domain = NA)
        break
      }
      #### Modify z and w for robust glm
      # Classical:
      # z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
      # w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])

      # Robust:
      phi <- 1 # TODO: estimate phi (for more families)
      ## see for e.g., section 5.3.3 of Robust Methods in Biostatistics (2009).
      varmu <- variance(mu)[good]
      mu.eta.val <- mu.eta(eta)[good]
      r <- (y - mu)[good] / sqrt(varmu) # Pearson's residuals
      psi <- .Call(C_psifun, r, c.psi, ipsi, 0L) # psi function evaluated on residuals
      psip <- .Call(C_psifun, r, c.psi, ipsi, 1L) # 1st derivative of psi function evaluated on residuals
      g2.val <- g2(mu)[good]
      varmup <- variancep(mu)[good]
      # alpha == 1 in classical (when considering Fisher scoring)
      alpha <- psi * (varmup / sqrt(phi * varmu) / 2 + sqrt(varmu) * g2.val * mu.eta.val) + psip * (1 / sqrt(phi) + r * varmup / sqrt(varmu) / 2)
      z <- psi * sqrt(varmu) / sqrt(phi) / mu.eta.val / alpha + eta[good]
      w <- mu.eta.val^2 * alpha / varmu

      #### End modification
      ngoodobs <- as.integer(nobs - sum(!good))
      ## call Fortran code via C wrapper
      fit <- .Call(C_Cdqrls, x[good, , drop = FALSE] * w, z * w,
                   min(1e-7, control$epsilon/1000), check=FALSE)
      if (any(!is.finite(fit$coefficients))) {
        conv <- FALSE
        warning(gettextf("non-finite coefficients at iteration %d", iter), domain = NA)
        break
      }
      ## stop if not enough parameters
      if (nobs < fit$rank)
        stop(sprintf(ngettext(nobs,
                              "X matrix has rank %d, but only %d observation",
                              "X matrix has rank %d, but only %d observations"),
                     fit$rank, nobs), domain = NA)
      ## calculate updated values of eta and mu with the new coef:
      start[fit$pivot] <- fit$coefficients
      eta <- drop(x %*% start)
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights))
      if (control$trace)
        cat("Deviance = ", dev, " Iterations - ", iter, "\n", sep = "")
      ## check for divergence
      boundary <- FALSE
      if (!is.finite(dev)) {
        if(is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
        warning("step size truncated due to divergence", call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit)
            stop("inner loop 1; cannot correct step size", call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        boundary <- TRUE
        if (control$trace)
          cat("Step halved: new deviance = ", dev, "\n", sep = "")
      }
      ## check for fitted values outside domain.
      if (!(valideta(eta) && validmu(mu))) {
        if(is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
        warning("step size truncated: out of bounds", call. = FALSE)
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit)
            stop("inner loop 2; cannot correct step size", call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
        }
        boundary <- TRUE
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace)
          cat("Step halved: new deviance = ", dev, "\n", sep = "")
      }
      ## check for convergence
      ## deviance might not be the best criterion for robust glm
      ## original from glm:
      # if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
      #   conv <- TRUE
      #   coef <- start
      #   break
      # } else {
      #   devold <- dev
      #   coef <- coefold <- start
      # }

      ## code modified from glmrobMqle:
      ## Check convergence: relative error < tolerance
      Dtheta <- start - coefold
      relE <- sqrt(sum(Dtheta^2)/max(1e-20, sum(coefold^2)))
      conv <- relE <= control$acc
      if(control$trace) {
        cat(sprintf("%3d | %*s | %12g\n", iter, width.th,
                    paste(mFormat(Dtheta, w.th.1),
                          collapse=" "), relE))
      }
      if(conv){
        coef <- start
        break
      } else {
        devold <- dev
        coef <- coefold <- start
      }

    } ##-------------- end IRLS iteration -------------------------------

    if (!conv)
      warning("glm.fit: algorithm did not converge", call. = FALSE)
    if (boundary)
      warning("glm.fit: algorithm stopped at boundary value", call. = FALSE)
    eps <- 10*.Machine$double.eps
    if (family$family == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps))
        warning("glm.fit: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
    }
    if (family$family == "poisson") {
      if (any(mu < eps))
        warning("glm.fit: fitted rates numerically 0 occurred", call. = FALSE)
    }
    ## If X matrix was not full rank then columns were pivoted,
    ## hence we need to re-label the names ...
    ## Original code changed as suggested by BDR---give NA rather
    ## than 0 for non-estimable parameters
    if (fit$rank < nvars) coef[fit$pivot][seq.int(fit$rank+1, nvars)] <- NA
    xxnames <- xnames[fit$pivot]
    ## update by accurate calculation, including 0-weight cases.
    residuals <-  (y - mu)/mu.eta(eta)
    ##        residuals <- rep.int(NA, nobs)
    ##        residuals[good] <- z - (eta - offset)[good] # z does not have offset in.
    fit$qr <- as.matrix(fit$qr)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
    }
    else Rmat <- fit$qr[1L:nvars, 1L:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  # for compatibility with lm, which has a full-length weights vector
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if(!EMPTY)
    names(fit$effects) <-
    c(xxnames[seq_len(fit$rank)], rep.int("", sum(good) - fit$rank))
  ## calculate null deviance -- corrected in glm() if offset and intercept
  wtdmu <-
    if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  ## calculate df
  n.ok <- nobs - sum(weights==0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if(EMPTY) 0 else fit$rank
  resdf  <- n.ok - rank
  ## calculate AIC
  aic.model <- aic(y, nobs, mu, weights, dev) + 2*rank
  ##     ^^ is only initialize()d for "binomial" [yuck!]
  list(coefficients = coef, residuals = residuals, fitted.values = mu,
       effects = if(!EMPTY) fit$effects, R = if(!EMPTY) Rmat, rank = rank,
       qr = if(!EMPTY) structure(fit[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
       family = family,
       linear.predictors = eta, deviance = dev, aic = aic.model,
       null.deviance = nulldev, iter = iter, weights = wt,
       prior.weights = weights, df.residual = resdf, df.null = nulldf,
       y = y, converged = conv, boundary = boundary)

}

#' @title Auxiliary for Controlling robust GLM Fitting
#' @param acc see
#' @param test.acc see
#' @param maxit see
#' @param tuning.psi see
#' @param psi see
#' @param epsilon see
#' @param trace see
#' @export
glmrob2.control <-
  function(acc = 1e-04, test.acc = "coef", maxit = 50, tuning.psi = NULL, psi = "bisquare", epsilon = 1e-08, trace = FALSE)
  {
    if (!is.numeric(acc) || acc <= 0)
      stop("value of acc must be > 0")
    if (test.acc != "coef")
      stop("Only 'test.acc = \"coef\"' is currently implemented")
    if (!is.numeric(maxit) || maxit <= 0)
      stop("maximum number of iterations must be > 0")
    if (!is.numeric(epsilon) || epsilon <= 0)
      stop("value of 'epsilon' must be > 0")

    # robust functions
    .Mpsi.tuning.default <- getFromNamespace(".Mpsi.tuning.default", ns = "robustbase")
    .regularize.Mpsi <- getFromNamespace(".regularize.Mpsi", ns = "robustbase")
    psi <- .regularize.Mpsi(psi)
    if(is.null(tuning.psi))
      tuning.psi <- .Mpsi.tuning.default(psi)
    list(acc = acc, test.acc = test.acc, maxit = maxit, tuning.psi = tuning.psi, psi = psi, epsilon = epsilon, trace = trace)
}
