# FUNCTIONS

# Function: subsample.gen

# subsample.gen(genotypes,nboots=1000,nsamps,loci)

# Description: Resampling routine to correct for unequal sampling / calculate allelic richness
# in comparison of allelic richness between two populations
# see: Leberg PL (2002) Estimating allelic richness: Effects of sample size and
#     bottlenecks. Molecular Ecology, 11, 2445-2449.
# SLOW! Be patient. Produces a lot of text as I don't know how to stop the summary
# function of adegenet to print out every iteration.

# Parameters:
# genotypes = adegenet genotypes object of the larger sample
# nboots = number of resamples (about 1000 should be ok, more doesn't hurt)
# nsamps = number of samples in the smaller sample
# loci = vector of generic (i.e. "L02") names of loci used in the bootstrap

# Example: Smaller dataset of 17 samples, loci from "loci_usa" (see below), 1000 resamples,
# 		  Dinaric bears used as a reference population:
# test=subsample.gen(dinaric.genotypes,1000,17,loci_usa)

#--------------------------------------------------------

# Function: runall
# runall(N,genotypes,loci,nboots=1000)

# Description: Runs corrections for unequal sampling for a subset of loci and
# a vector of sample sizes.
# N = vector of sample sizes
# genotypes = genind object of genotypes (adegenet package)
# loci = vector of generic locus names
# nboots = number of bootstraps
# verbose = default FALSE, if TRUE, will print summary for each iteration
#' @export
#' @importFrom adegenet nInd locNames "indNames<-" genind2df df2genind
#' @importFrom stats sd
subsample.gen <- function(genotypes, nboots = 1000, nsamps, loci, verbose = FALSE) {
  # resampling routine to correct for unequal sampling / calculate allelic richness
  # in comparison of allelic richness between two populations
  # see: Leberg PL (2002) Estimating allelic richness: Effects of sample size and
  #     bottlenecks. Molecular Ecology, 11, 2445-2449.
  # SLOW! Be patient. Produces a lot of text as I don't know how to stop the summary
  # function of adegenet to print out every iteration.

  # genotypes = adegenet genotypes object of the larger sample
  # nboots = number of resamples (about 1000 should be ok, more doesn't hurt)
  # nsamps = number of samples in the smaller sample
  # loci = vector of generic (i.e. "L02") names of loci used in the bootstrap

  # Example: Smaller dataset of 17 samples, loci from "loci_usa" (see below), 1000 resamples,
  #     Dinaric bears used as a reference population:
  # test=subsample.gen(dinaric.genotypes,1000,17,loci_usa)
  ngens <- nInd(genotypes)
  Al <- NULL
  SEAl <- NULL
  Hex <- NULL
  SEHex <- NULL
  Hobs <- NULL
  SEHobs <- NULL
  genotypes <- genotypes[, loc = loci]

  for (i in 1:nboots) {
    samp <- sample(1:ngens, nsamps, replace = TRUE)
    gen.sample <- genotypes[samp]

    indNames(gen.sample) <- 1:nInd(gen.sample)
    gen.sample <- suppressWarnings(genind2df(gen.sample, sep = " "))

    gen.sample <- suppressWarnings(df2genind(gen.sample, sep = " "))
    # summary as of adegenet (=1.4.2) is not exported from its
    # namespace. we need to specify exactly where summary comes
    # from by using ::
    summ <- adegenet::summary(gen.sample, verbose = FALSE)

    Al <- rbind(Al, mean(summ$loc.n.all))
    SEAl <- rbind(SEAl, sd(summ$loc.n.all) / sqrt(length(summ$loc.n.all)))
    Hex <- rbind(Hex, mean(summ$Hexp))
    SEHex <- rbind(SEHex, sd(summ$Hexp) / sqrt(length(summ$Hexp)))
    Hobs <- rbind(Hobs, mean(summ$Hobs))
    SEHobs <- rbind(SEHobs, sd(summ$Hobs) / sqrt(length(summ$Hobs)))
  }
  # for
  #*** OUTPUTS
  A <- mean(Al) # mean number of alleles
  SEA <- mean(SEAl)
  # sdA = sd(Al) #standard deviation of the number of alleles
  He <- mean(Hex) # mean expected heterozygosity
  SEHe <- mean(SEHex)
  # sdHe = sd(Hex) #standard deviation of expected heterozygosity
  Ho <- mean(Hobs)
  SEHo <- mean(SEHobs)
  # sdHo = sd (Hobs)
  out <- data.frame(A, SEA, He, SEHe, Ho, SEHo)

  if (verbose == TRUE) {
    print(locNames(genotypes))
  }

  out
}

#' @export
runall <- function(N, genotypes, loci, nboots = 1000) {
  # Runs corrections for unequal sampling for a subset of loci and a vector of sample sizes
  # N = vector of sample sizes
  # genotypes = genind object of genotypes (adegenet package)
  # loci = vector of generic locus names
  # nboots = number of bootstraps
  out <- NULL
  for (i in 1:length(N)) {
    resamp <- subsample.gen(genotypes, nboots, N[i], loci)
    out <- rbind(out, data.frame(Nsamp = N[i], resamp))
  }
  return(as.data.frame(out))
}

#' @export
calcDivRat <- function(ref, SEref, obs, SEobs, type = "U") {
  # ref=resampled reference population parameter
  # SEref = standard error of ref
  # obs = studied population parameter
  # SEobs = standard error of obs
  # type = type of parameter, "A" for allelic diversity, "He" for
  # expected heterozygosity. Default = "U" (unknown)
  ratio <- obs / ref
  SE.ratio <- sqrt((ratio^2) * (((SEobs / obs)^2) + ((SEref / ref)^2)))
  out <- data.frame(ratio, SE.ratio)
  names(out) <- c(paste(type, "r", sep = ""), paste("SE", type, "r", sep = ""))
  return(out)
}
