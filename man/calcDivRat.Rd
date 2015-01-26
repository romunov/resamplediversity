\name{calcDivRat}
\alias{calcDivRat}
\title{
Calculate genetic diversity ratio and its standard error.
}
\description{
Calculates diversity ratio index (Ar or Hert) between the reference population and studied population, calculates its standard error as a standard error of division (see the manuscript for details).
}
\usage{
calcDivRat(ref, SEref, obs, SEobs, type = "U")
}
\arguments{
  \item{ref}{
Sample-size corrected diversity index in reference population (Heterozygosity of Allelic Richness)
}
  \item{SEref}{
Standard error of ref.
}
  \item{obs}{
Diversity index in the studied population (Heterozygosity of Allelic Richness)
}
  \item{SEobs}{
Standard error of obs.
}
  \item{type}{
TYPE of the diversity index - A for allelic richness, He for expected heterozygosity. Changes the names of columns of the returned data frame. "U" for unknow is set as default if the user doesn't define this parameter.
}
}
\value{
Returns a dataframe with two columns with the ratio and its standard error. Column names depend on the "type" defined in the function.
}
\references{
Skrbinšek T, Jelenčič M, Waits LP, Potočnik H, Kos I, Trontelj P(2012) Using a reference population yardstick to calibrate and compare genetic diversity reported in different studies: an example from the brown bear. Heredity, In press.
}
\author{
Tomaž Skrbinšek \email{tomaz.skrbinsek@gmail.com}
}
\examples{
# For usage, see vignette
vignette("resamplediversity")
}