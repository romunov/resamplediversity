\name{bear.diversity}
\alias{bear.diversity}
\docType{data}
\title{
Table of brown bear diversity data from a number of studies around the world.
}

\description{A \code{data.frame} of results of pertinent studies. See vignette for full citations.}

\usage{data(bear.diversity)}
\format{
  A data frame with 31 observations on the following 7 variables.
  \describe{
    \item{\code{Population}}{Studied population.}
    \item{\code{N}}{Number of samples analyzed.}
    \item{\code{Study}}{Number of studied - look at the dataframe "included.studies", included.studies$ID}
    \item{\code{A}}{Allelic diversity}
    \item{\code{SEA}}{Standard error of allelic diversity}
    \item{\code{He}}{Expected heterozygosity.}
    \item{\code{SEHe}}{Standard error of expected heterozygosity.}
  }
}

\source{
Skrbinšek T, Jelenčič M, Waits LP, Potočnik H, Kos I, Trontelj P(2012) Using a reference population yardstick to calibrate and compare genetic diversity reported in different studies: an example from the brown bear. Heredity, In press.
}

\examples{
data(bear.diversity)
bear.diversity
}
\keyword{datasets}
