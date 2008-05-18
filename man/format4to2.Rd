\name{format4to2}
\alias{format4to2}
\title{Changes data for sdt to data for mlmsdt}
\description{
  Changes data in the form of hits, false alarms, misses, and correct
  rejections (which is suitable for the sdt function) for each individual
  into a file with one trial per line suitable for mlmsdt. See Wright,
  Horry and Skagerberg (in press, Behavior Research Methods).
}
\usage{
format4to2(x, cnames = c("subno", "isold", "saysold"), code = 0.5, ...)
}
\arguments{
  \item{x}{dataframe with subject number, hits, false alarms, misses, and CR}
  \item{cnames}{names to output}
  \item{code}{code for isold and saysold}
  \item{\dots}{other parameters passed}}
\details{
  code -
    default is -.5 for new and +.5 for old for isold. Put 1 in for 0 and 1.
    saysold is 0 and 1
}
\value{
  A three column dataframe with length equal the number of trials is
  produced.
}
\references{
Wright, D.B., Horry, R., \& Skagerberg, E.M. (in press). Functions for
traditional and multilevel approaches to signal detection theory.
\emph{Behavior Research Methods}. }
\author{Daniel B. Wright}
\examples{
format4 <- rbinom(100,25,.5)
dim(format4) <- c(25,4)
format4 <- cbind(1:25,format4)
dim(format4)
format4[1:5,]
format2 <- format4to2(format4)
dim(format2)
format2[1:5,]
}