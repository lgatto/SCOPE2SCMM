
## This script is divided in 3 parts
##  - Get the scRNA-Seq data
##  - Get the SCP data 
##  - Combine the 2 modalities in MAE

## Required packages
library(MultiAssayExperiment)
library(HDF5Array)
library(SingleCellExperiment)
library(BiocFileCache)


####---- Get the scRNA-Seq data ----####


## Steps:
## 1. Add the scRNASeq matrices to BiocFileCache
## 2. Reads the count matrices
## 3. Adds the counts matrices to BiocFileCache as a hdf5 resource
## 4. Creates the SingleCellExperiment objects
## 5. 
## Steps 1 to 3 need to be run only once.

## See also https://bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/UsingHDF5Array.html

bfc <- BiocFileCache()

## -------------------------------------- ##
## 1. Add web resources to BiocFileCache
## -------------------------------------- ##
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
## 4. Create SingleCellExperiment objects
## ------------------------------------ ##
hdf5file <- bfcquery(bfc, "GSE142392")$rpath

rna1 <- HDF5ArraySeed(file = hdf5file, name = "rna1")
sc1 <- SingleCellExperiment(assays = DelayedArray(rna1))

rna2 <- HDF5ArraySeed(file = hdf5file, name = "rna2")
sc2 <- SingleCellExperiment(assays = DelayedArray(rna2))

## The 2 scRNA-Seq datasets
sc1
sc2


####---- Get the SCP data ----####

## Download file from GoogleDrive (this will very soon be accessible
## from an EH package): specht2019v3.Rda
## https://drive.google.com/file/d/1OLVVRxvWVttsQF142dvgp4ZcHJEDdGRI/view?usp=sharing
load("./specht2019v3.Rda")
## The data contains different levels of proteomics information, you 
## should be only interested in the protein level (final processed 
## data).
scp <- specht2019v3[["proteins"]]

## The SCP datasets
scp


####---- Combine the 2 modalities in MAE ----####


mae <- MultiAssayExperiment(
    experiments = ExperimentList(scRNAseq1 = sc1,
                                 scRNAseq2 = sc2, 
                                 scp = scp)
)