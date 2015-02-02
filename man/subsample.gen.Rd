\name{subsample.gen}
\alias{subsample.gen}

\title{
Correct for unequal sampling / calculate allelic richness
}
\description{
Resampling routine to correct for unequal sampling / calculate allelic richness in comparison of allelic richness between two populations~
}
\usage{
subsample.gen(genotypes, nboots = 1000, nsamps, loci, verbose = FALSE)
}
\arguments{
  \item{genotypes}{
Genind object of genotypes (\link[adegenet]{adegenet} package) for the reference popoulation.
}
  \item{nboots}{
Number of bootstrap resamples. Default=1000.
}
  \item{nsamps}{
Number of samples to use in each subsample - equal to the sample size in the compared population.
}
  \item{loci}{
Subset of loci. Vector of generic locus names in the genotypes object c("L01","L03",...) common in both the reference and the compared populations.
}
  \item{verbose}{
If \code{TRUE}, print summary for each iteration. Can be very verbose.
}
}

\value{
Data frame with parameter means over all subsamples.
A: allelic diversity
SEA: standard error of A
He: expected deterozygosity
SEHe: standard error of expected heterozygosity
Ho: observed deterozygosity
SEHo: standard error of observed heterozygosity
}
\references{
Skrbinšek T, Jelenčič M, Waits LP, Potočnik H, Kos I, Trontelj P(2012) Using a reference population yardstick to calibrate and compare genetic diversity reported in different studies: an example from the brown bear. Heredity, In press.

Leberg PL (2002) Estimating allelic richness: Effects of sample size and bottlenecks. Molecular Ecology, 11, 2445-2449.
}
\author{
Tomaž Skrbinšek \email{tomaz.skrbinsek@gmail.com}
}
\note{
SLOW! Be patient. Produces a lot of text as I don't know how to stop the summary function of adegenet package to print out every iteration.
}

\seealso{
\code{\link[resamplediversity:runall]{runall}}, \link[adegenet:genind]{genind}
}
\examples{
# For examples, see vignette
vignette("resamplediversity")
}
