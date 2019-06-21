\name{ddn}
\alias{ddn}
\docType{data}
\title{
Astrophysical S-factor versus center-of-mass energy for the d(d,n)reaction.
}
\description{
Data set compiled from  Leonard et al. (2006), Greife et al. (1995), Brown et al. (1990), Krauss et al. (1987) (B)  and Krauss et al. (1987) (M).
}
\usage{data("ddn")}
\format{
  A data frame with 52 observations on the following 3 variables.
  \describe{
    \item{\code{E}}{Energy}
    \item{\code{S}}{S-factor}
    \item{\code{Stat}}{S-factor uncertainties}
  }
}
\details{
Colum  Syst displayes the systhematic uncertainty of each experiment given by column dat.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ddn)
head(ddn)
}
\keyword{datasets}
