\name{singlelevel}
\alias{singlelevel}
\title{Using a GLM with a single participant's data}
\description{
  Called by mlmsdt when only a single unique subject number is in the data file.
  It is assumed this is not called directly by users (use glm instead).
}
\usage{
singlelevel(isold, sold, covs = {
}, lk = "logit", int = FALSE, modify = TRUE)
}
\arguments{
  \item{isold}{whether item is old or new}
  \item{sold}{whether person says old or new}
  \item{covs}{list of covariates}
  \item{lk}{link function, "probit" or "logit" usually}
  \item{int}{interactions between covariates}
  \item{modify}{centering covariates and -5, +.5 coding for isold}
}
\value{
  a glm.object is returned and values for d' or lnOR printed to screen
}
\author{Daniel. B. Wright}
\examples{
x <- runif(100,0,1)
singlelevel(rbinom(100,1,x),rbinom(100,1,x))
}
