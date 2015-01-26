\name{runall}
\alias{runall}

\title{
Batch run subsampling corrections for allelic richness w/ function subsample.gen 
}
\description{
Batch wrapper for subsample.gen function. Batch runs corrections for unequal sampling for a subset of loci and a vector of sample sizes. Streamlines usage of subsample.gen function to process a number of populations of different sample size that were analyzed using the same panel of markers.
}
\usage{
runall(N, genotypes, loci, nboots = 1000)
}

\arguments{
  \item{N}{
Vector of sample sizes. Sample sizes are the numbers of samples to use in each subsample - equal to the sample sizes in the set of compared populations.
}
  \item{genotypes}{
Genind object of genotypes (\link[adegenet]{adegenet} package) for the reference popoulation.
}
  \item{loci}{
Subset of loci. Vector of generic locus names in the genotypes object c("L01","L03",...) common in both the reference and the compared populations.
}
  \item{nboots}{
Number of bootstrap resamples. Default=1000.
}
}

\value{
Dataframe with summarized data (parameter means over all subsamples) for each sample size.
N: sample size
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
\code{\link[resamplediversity:subsample.gen]{subsample.gen}}, \link[adegenet:genind]{genind}
}
\examples{
# For examples, see vignette
vignette("resamplediversity")
}

