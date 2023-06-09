% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmrob2.R
\name{glmrob2}
\alias{glmrob2}
\alias{glmrob2.fit}
\title{Robust Fitting of Generalized Linear Models (under development)
\code{glmrob2} is used to fit generalized linear models.
The iterative re-weighted least squared is used for fitting, similarly
to \code{\link[stats]{glm}}. Currently only \code{\link[stats]{binomial}} with \code{logit}
link is supported. Be aware estimators are not Fisher consistent.}
\usage{
glmrob2(
  formula,
  family = gaussian,
  data,
  weights,
  subset,
  na.action,
  start = NULL,
  etastart,
  mustart,
  offset,
  control = list(...),
  model = TRUE,
  x = FALSE,
  y = TRUE,
  contrasts = NULL,
  ...
)

glmrob2.fit(
  x,
  y,
  weights = rep(1, nobs),
  start = NULL,
  etastart = NULL,
  mustart = NULL,
  offset = rep(0, nobs),
  family = gaussian(),
  control = list(),
  intercept = TRUE
)
}
\arguments{
\item{formula}{an object of class formula (see \code{\link[stats]{formula}}).}

\item{family}{see}

\item{data}{see}

\item{weights}{see}

\item{subset}{see}

\item{na.action}{see}

\item{start}{see}

\item{etastart}{see}

\item{mustart}{see}

\item{offset}{see}

\item{control}{see}

\item{model}{see}

\item{x}{see}

\item{y}{see}

\item{contrasts}{see}

\item{...}{see}

\item{intercept}{see}
}
\description{
Robust Fitting of Generalized Linear Models (under development)
\code{glmrob2} is used to fit generalized linear models.
The iterative re-weighted least squared is used for fitting, similarly
to \code{\link[stats]{glm}}. Currently only \code{\link[stats]{binomial}} with \code{logit}
link is supported. Be aware estimators are not Fisher consistent.
}
\examples{
## Binary response --------------
data(vaso)

## classical fit
Vfit1 <- glm(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso)
coef(Vfit1)

## robust fit (robustbase)
Vfit2 <- glmrob(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso,
                method="Mqle", control=glmrobMqle.control(tcc=3.5))
coef(Vfit2)

## robust fit 2 (new method, allows for different weight functions)
Vfit3 <- glmrob2(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso)
coef(Vfit3)
}
\author{
Samuel Orso
}
