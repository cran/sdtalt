\name{sdt}
\alias{sdt}
\title{Measures from SDT}
\description{
  Several statistics from signal detection theory are calculated for each
  individual. The user can define their own, too. Summary statistics (mean,
  trimmed mean, standard deviation) for the sample can be calculated and
  BCa confidence intervals produced.
}
\usage{
sdt(hits, fas, misses, cr, flat = 0, pmeans = FALSE, meas = "all", wk = 0.5, runboot = FALSE, trim = 0, R = 2000, confl = 0.95, newst = NULL, bound = FALSE, ...)
}
\arguments{
  \item{hits}{vector of hits}
  \item{fas}{vector of false alarms}
  \item{misses}{vector of misses}
  \item{cr}{vector of correct rejections}
  \item{flat}{value added to all cells, default = 0, .5 is common}
  \item{pmeans}{print the means and standard deviations}
  \item{meas}{which measures. See below. Default is all 15.}
  \item{wk}{if kappa used, the weighting value}
  \item{runboot}{whether to run bootstrap for BCa confidence intervals}
  \item{trim}{how much to trim the means. Default = 0, .2 is common}
  \item{R}{how many replicates for the bootstrap. Default = 2000.}
  \item{confl}{confidence level for intervals, between 0 and .999}
  \item{newst}{name of any user-defined statistic is to be used. See below.}
  \item{bound}{whether to bound infinite value to nearest finite value.}
  \item{\dots}{other parameters passed to the function}
}
\details{
 meas can take a list of up to 15 statistics.
 \tabular{ll}{
  HR     \tab  Hit rate \cr
  FAR    \tab  False alarm rate \cr
  d      \tab  d' \cr
  csdt   \tab  C \cr
  A      \tab  A' \cr
  B      \tab  B'' \cr
  lnbeta \tab  lnbeta \cr
  beta   \tab	 beta  \cr
  OR     \tab  Odds ratio \cr
  lnOR	 \tab  lnOR  \cr
  kappa  \tab  Weighted kappa \cr
  phi    \tab  phi  \cr
  Q      \tab  Yule's Q  \cr
  eta    \tab  Choice-theory measure eta \cr
  PC     \tab  Proportion correct }
  newst is the name of a user-defined statistic that must have 4 arguments,
 for hits, false alarms, misses, and correct rejections, in this order. \cr If
 any user wants a statistics added to this function, please contact Dan Wright.}
\value{
 If pmeans=TRUE then sample statistics are printed to screen. A dataframe of
 the same length as input is created with statistics for each of those listed
 in the meas option (15 by default).
}
\references{
 Wright, D.B., Horry, R., \& Skagerberg, E.M. (in press). Functions for
 traditional and multilevel approaches to signal detection theory.
 \emph{Behavior Research Methods}.}
\author{Daniel B. Wright}
\examples{
hits <- rbinom(100,25,.6)
fa <- rbinom(100,25,.2)
miss <- rbinom(100,25,.4)
cr <- rbinom(100,25,.7)
sdtout <- sdt(hits,fa,miss,cr)
sdtout[1:3,]
sdt(hits,fa,miss,cr,meas=c("d","A"))[1:3,]
sdt(hits,fa,miss,cr,meas=c("d","A"),flat=.5)[1:3,]
HC <- function(hits,fas,misses,cr)
   return((hits+cr-fas-misses)/(hits+cr+fas+misses))
HCsqrt <- function(hits,fas,misses,cr)
   return(sqrt((hits+cr-fas-misses)/(hits+cr+fas+misses)))
HCstats <- sdt(hits,fa,miss,cr, meas={},newst=c(HC,HCsqrt))
HCstats[1:3,]
sdt(hits,fa,miss,cr, meas=c("d","A"), pmeans=TRUE)[1:3,]
sdt(hits,fa,miss,cr,meas=c("d","A"), pmeans=TRUE,trim=.2)[1:3,]
sdt(hits,fa,miss,cr, meas=c("A"),pmeans=TRUE,trim=.2,runboot=TRUE,confl=.90)[1:3,]
}


