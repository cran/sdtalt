\name{mlmsdt}
\alias{mlmsdt}
\title{Using multilevel modeling for repeated 2x2 tables}
\description{
  Using Bates' lme4 package (2007) to run multilevel generalized linear models
  as an alternative to the traditional approaches to sdt. Details in Wright,
  Horry and Skagerberg (in press, Behavior Research Methods).
}
\usage{
mlmsdt(subno, isold, sold, covs = NULL, lk = "logit", vardiff = TRUE, modify = TRUE, int = FALSE, item = NULL, ...)
}
\arguments{
  \item{subno}{unique subject number for each individual}
  \item{isold}{whether the person says old (1) or new (0)}
  \item{sold}{whether the item is old (1) or new (0)}
  \item{covs}{a list of covariates}
  \item{lk}{link function, default = "logit" but any that can be used in lmer}
  \item{vardiff}{whether to allow different variances for old and new items}
  \item{modify}{used for centering covariates}
  \item{int}{to look at interactions among covariates}
  \item{item}{if item numbers available, a random intercept for them included}
  \item{\dots}{other parameters can be passed}
}
\details{
 subno -
    if only one unique subno is found, the function singlelevel is run.
    This is because lmer requires more than one subject. \cr
 isold and sold -
    the assumption is that these have the old values higher than the low values,
    but as long as they are in the same direction the d' and lnOR output will
    be fine. The 0/1 dummy variable coding is for ease of interpreting the
    coefficients. modify also affects these. \cr
 lk -
    it is assumed only logit and probit are used, so the output
    is designed for this but others can be used.\cr
 item -
    this allows you to run item response models. Baayen et al. (2008) is worth
    consulting if using this.\cr
 More complex models can be run accessing the lmer function directly.
}
\value{
  Output to the screen includes values for diagnosticity (lnOR or d'). \cr An
  S4 object is created.
}
\references{
Baayen, R.H., Davidson, D.J., \& Bates, D.M. (2008). Mixed-effects modeling
with crossed random effects for subjects and items.
\emph{Journal of Memory and Language}.\cr
Bates, D. (2007). lme4: Linear mixed-effects models using S4 classes.
Version 0.99875-9. http://cran.r-project.org/.\cr
Wright, D.B., Horry, R., \& Skagerberg, E.M. (in press). Functions for
traditional and multilevel approaches to signal detection theory.
\emph{Behavior Research Methods}.
}
\author{Daniel B. Wright}
\seealso{sdt}
\examples{
set.seed(47)
id <- rep(1:10,each=10)
isold <- rbinom(100,1,.5)
covariate <- runif(100,0,.3)
item <- rep(1:10,10)
sayold <- rbinom(100,1,isold/4+covariate+item/50+id/50)
mlmsdt(id,isold,sayold,covariate,item=item)
}