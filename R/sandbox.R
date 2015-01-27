install.packages("resamplediversity_1.0.tar.gz", type = "source")
library(resamplediversity)

loci_na <- c("L02", "L03", "L04", "L06", "L07",
        "L08", "L09", "L10")

data(dinaric.genotypes)
#load("./resamplediversity/data/dinaric.genotypes.rda")

resampled.ar <- subsample.gen(genotypes = dinaric.genotypes,
        nboots = 1000,
        nsamps = 50,
        loci = loci_na)