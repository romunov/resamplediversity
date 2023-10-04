\name{dinaric.genotypes}
\alias{dinaric.genotypes}
\docType{data}
\title{
Genotypes of brown bears (\emph{Ursus arctos}) from Northern Dinaric Mountains
}
\description{
The dataset contains genotypes of 513 brown bears genotyped on 20 microsatellite loci. They were obtained by genotyping tissue samples of legally killed bears between 2003 and 2008. The genotypes are mostly complete, there are missing data on a single locus in two samples. Species: \emph{Ursus arctos}. Geographic location: Slovenia, 45°29'N - 46°32'N; 13°44'E - 15°25'E
}
\usage{data(dinaric.genotypes)}

\details{
This dataset mostly corresponds to the dataset used in the manuscript:
Skrbinšek T, Jelenčič M, Waits L, Kos I, Jerina K, Trontelj P (2012). Monitoring the effective population size of a brown bear (\emph{Ursus arctos}) population using new single-sample approaches. Molecular Ecology 21(4): 862-875. Data: doi: 10.5061/dryad.22rm1728.
This present dataset doesn't have age data for the animals, but includes genotypes of three additional individuals for which we had no age data available.

Due to changes in newever versions of the \code{adegenet} package (>=2.0), the \code{genind} object had to be ported. During the process, some (empty) alleles were not transferred. This shouldn't influence the results from the vignette since these faulty microsatellites (Mu05, Mu09) were not included in the subset of analyzed loci. \emph{See original release of the package on GitHub for the data and code that was bundled that accompanied the original paper.}
}
\source{
Skrbinšek T, Jelenčič M, Waits LP, Potočnik H, Kos I, Trontelj P (2012) Using a reference population yardstick to calibrate and compare genetic diversity reported in different studies: an example from the brown bear. Heredity, In press.
}
\seealso{
\link[adegenet:genind]{genind}
}
\examples{
data(dinaric.genotypes)
# summary(dinaric.genotypes)
}
\keyword{datasets}
