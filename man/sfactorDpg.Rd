% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sfactorDpg.R
\name{sfactorDpg}
\alias{sfactorDpg}
\title{Estimate Astrophysical S-factor}
\format{\describe{
\item{x}{
The function has two  arguments:  ecm, a.scale}
}}
\usage{
sfactorDpg(ecm = ecm, a.scale = a.scale)
}
\arguments{
\item{ecm}{ecm}

\item{a.scale}{a.scale}
}
\value{
S-factor
}
\description{
Provides a confusion matrix of classification statistics following logistic regression.
}
\examples{
library(nuclear)

N <- 300
obsx1 <- exp(seq(log(1e-3), log(1),length.out=N))
plot(obsx1,sfactorDpg(obsx1),
col="red",cex=1.25,type="l",ylab="S-factor",xlab="E",log="x")

}
\author{
Rafael de Souza, UNC,  and Christian Illiadis, UNC
}
\keyword{S-factor}
