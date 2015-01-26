\name{included.studies}
\alias{included.studies}
\docType{data}
\title{
Table of genetic studies of brown bears included in the manuscript.
}
\description{A list of pertinent studies used in our paper. See vignette for full citations.}
\usage{data(included.studies)}
\format{
  A data frame with 10 observations on the following 7 variables.
  \describe{
    \item{\code{ID}}{Study ID... the same as in "Study" column of bear.diversity dataset.}
    \item{\code{Reference}}{Literature reference. See manuscript references for details.}
    \item{\code{GeoArea}}{Geographic area studied.}
    \item{\code{Aim}}{Why was the study conducted.}
    \item{\code{NP}}{Number of populations studied.}
    \item{\code{LocUsed}}{Number of microsatellite loci used.}
    \item{\code{LocCommon}}{Number of loci in common with the reference population in Northern Dinaric Mountains.}
  }
}
\source{
Skrbinšek T, Jelenčič M, Waits LP, Potočnik H, Kos I, Trontelj P(2012) Using a reference population yardstick to calibrate and compare genetic diversity reported in different studies: an example from the brown bear. Heredity, In press.
}
\examples{
data(included.studies)
included.studies
}
\keyword{datasets}
