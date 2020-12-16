## This script
 ## 1. Add the scRNASeq matrices to BiocFileCache
## 2. Reads the count matrices
## 3. Adds the counts matrices to BiocFileCache as a hdf5 resource
## 4. Creates the SingleCellExperiment objects
## Steps 1 to 4 need to be run only once.

## See also https://bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/UsingHDF5Array.html

library(MultiAssayExperiment)
library(HDF5Array)
library(SingleCellExperiment)
library(BiocFileCache)

## -------------------------------------- ##
## 1. Add web resources to BiocFileCache
## -------------------------------------- ##
## bfc <- BiocFileCache()
## url1 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4226nnn/GSM4226877/suppl/GSM4226877_rna_data_Bio_Replicate_1.csv.gz"
## url2 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4226nnn/GSM4226878/suppl/GSM4226878_rna_data_Bio_Replicate_2.csv.gz"
## bfcrpath(bfc, url1)
## bfcrpath(bfc, url2)


## ---------------------- ##
## 2. Read count matrices
## ---------------------- ##
## m1 <- data.table::fread(file = bfcquery(bfc, "GSM4226877")$rpath,
##                         sep = ",", header = TRUE)
## rn <- m1[[1]]
## m1 <- as.matrix(m1[, -1])
## rownames(m1) <- rn
## m2 <- data.table::fread(file = bfcquery(bfc, "GSM4226878")$rpath,
##                         sep = ",", header = TRUE)
## rn <- m2[[1]]
## m2 <- as.matrix(m2[, -1])
## rownames(m2) <- rn

## ------------------------------------------------------- ##
## 3. Write count matrices to a new BiocFileCache resource
## ------------------------------------------------------- ##
## hdf5file <- bfcnew(bfc, "GSE142392", ext=".h5")
## writeHDF5Array(m1,
##                filepath = hdf5file,
##                name = "rna1",
##                with.dimnames = TRUE)
## writeHDF5Array(m2,
##                filepath = hdf5file,
##                name = "rna2",
##                with.dimnames = TRUE)

## ------------------------------------ ##
## Create SingleCellExperiment objects
## ------------------------------------ ##
hdf5file <- bfcquery(bfc, "GSE142392")$rpath

rna1 <- HDF5ArraySeed(file = hdf5file, name = "rna1")
sc1 <- SingleCellExperiment(assays = DelayedArray(rna1))
sc1

rna2 <- HDF5ArraySeed(file = hdf5file, name = "rna2")
sc2 <- SingleCellExperiment(assays = DelayedArray(rna2))
sc2
