\name{format2to4}
\alias{format2to4}
\title{Changes data for sdt to mlmsdt form}
\description{
  Takes three variables (subject number, whether an item is old, and what
  the person says) and changes them into a dataframe with hits, false alarms,
  misses and correct rejections for each subject.
}
\usage{
format2to4(subno, isold, sold, cnames = c("subno", "hits", "fas", "misses", "crs"), ...)
}
\arguments{
  \item{subno}{unique subject number for each individual}
  \item{isold}{whether the item is old - higher number old}
  \item{sold}{whether the person saysold - higher number old}
  \item{cnames}{names for the five variables in the output dataframe}
  \item{\dots}{other parameters passed}
}
\value{
  An nx5 dataframe
}
\references{
Wright, D.B., Horry, R., \& Skagerberg, E.M. (in press). Functions for
traditional and multilevel approaches to signal detection theory.
\emph{Behavior Research Methods}.
}
\author{Daniel B. Wright}
\examples{
madeup <- data.frame(rep(1:10,each=10),rep(0:1,50),rbinom(100,1,.5))
colnames(madeup) <- c("subno","isold","sayold")
madeup[1:5,]
madeup4 <- format2to4(madeup[,1],madeup[,2],madeup[,3])
madeup4
}