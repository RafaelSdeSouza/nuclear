% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NumRateHe3dp.R
\name{NumRateHe3dp}
\alias{NumRateHe3dp}
\title{Estimate reaction rates}
\format{\describe{
\item{x}{
The function has six arguments: e0, gi, gf, ri, rf, T9}
}}
\usage{
NumRateHe3dp(e0 = e0, gi = gi, gf = gf, ri = ri, rf = rf, T9 = T9)
}
\arguments{
\item{e0}{e0}

\item{gi}{gi}

\item{gf}{gf}

\item{ri}{ri}

\item{rf}{rf}

\item{T9}{T9}
}
\value{
nuclear_rate
}
\description{
Calculate numerical reaction rates
}
\examples{
library(nuclear)

 NumRateHe3dp(0.35,1.0085,0.025425,5,5,10)

}
\author{
Rafael de Souza, UNC,  and Christian Illiadis, UNC
}
\keyword{Nuclear}
\keyword{rate}
