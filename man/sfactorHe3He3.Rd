% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sfactorHe3He3.R
\name{sfactorHe3He3}
\alias{sfactorHe3He3}
\title{Estimate Astrophysical S-factor}
\format{\describe{
\item{x}{
The function has 5  arguments: ecm,S0, S1, S2,i.screen}
}}
\usage{
sfactorHe3He3(ecm = ecm, S0 = S0, S1 = S1,S2 = S2, i.screen = i.screen)
}
\arguments{
\item{ecm}{ecm}

\item{S0}{S0}

\item{S1}{S1}

\item{S2}{S2}

\item{i.screen}{i.screen}
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
obsx1 <- exp(seq(log(1e-2), log(1),length.out=N))
plot(obsx1,sfactorHe3He3(obsx1,5.14,-2.69,2.14,325*1e-6),
col="red",cex=1.25,type="l",ylab="S-factor",xlab="E",log="x")

}
\author{
Rafael de Souza, UNC,  and Christian Illiadis, UNC
}
\keyword{S-factor}
