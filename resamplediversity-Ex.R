pkgname <- "resamplediversity"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('resamplediversity')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("bear.diversity.rda")
### * bear.diversity.rda

flush(stderr()); flush(stdout())

### Name: bear.diversity
### Title: Table of brown bear diversity data from a number of studies
###   around the world.
### Aliases: bear.diversity
### Keywords: datasets

### ** Examples

data(bear.diversity)
bear.diversity



cleanEx()
nameEx("calcDivRat")
### * calcDivRat

flush(stderr()); flush(stdout())

### Name: calcDivRat
### Title: Calculate genetic diversity ratio and its standard error.
### Aliases: calcDivRat

### ** Examples

# For usage, see vignette
vignette("resamplediversity")



cleanEx()
nameEx("dinaric.genotypes.rda")
### * dinaric.genotypes.rda

flush(stderr()); flush(stdout())

### Name: dinaric.genotypes
### Title: Genotypes of brown bears (_Ursus arctos_) from Northern Dinaric
###   Mountains
### Aliases: dinaric.genotypes
### Keywords: datasets

### ** Examples

data(dinaric.genotypes)
#summary(dinaric.genotypes)



cleanEx()
nameEx("included.studies.rda")
### * included.studies.rda

flush(stderr()); flush(stdout())

### Name: included.studies
### Title: Table of genetic studies of brown bears included in the
###   manuscript.
### Aliases: included.studies
### Keywords: datasets

### ** Examples

data(included.studies)
included.studies



cleanEx()
nameEx("runall")
### * runall

flush(stderr()); flush(stdout())

### Name: runall
### Title: Batch run subsampling corrections for allelic richness w/
###   function subsample.gen
### Aliases: runall

### ** Examples

# For examples, see vignette
vignette("resamplediversity")



cleanEx()
nameEx("subsample.gen")
### * subsample.gen

flush(stderr()); flush(stdout())

### Name: subsample.gen
### Title: Correct for unequal sampling / calculate allelic richness
### Aliases: subsample.gen

### ** Examples

# For examples, see vignette
vignette("resamplediversity")



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
