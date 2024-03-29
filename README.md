# CamCNV
CamCNV pipeline for calling rare CNVs from Illumina data.
For further details see the paper in Genetic Epidemiology -  "Detecting rare copy number variants from Illumina genotyping arrays with the CamCNV pipeline: Segmentation of z-scores improves detection and reliability" (https://pubmed.ncbi.nlm.nih.gov/33020983/) and the analysis of a large breast cancer dataset - "Rare germline copy number variants (CNVs) and breast cancer risk" (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8766486/)

The basic steps of the pipeline as outlined in the example code are:
1. Load Log R Ratio intensities (LRR) and B Allele Frequencies (BAF), sorted by chromosome and position, from Illumina gentotyping into a NetCDF data store
2. Run a principal component adjustment on the LRR
3. Calculate the mean and standard deviation of the LRR for each probe across all samples and convert the LRR into z-scores
4. For each sample segment the z-scores using DNACopy (https://bioconductor.org/packages/release/bioc/html/DNAcopy.html) and identify potential CNVs from the mean z-scores of the segments
5. Generate additional QC scores for each CNVs

The pipeline uses a NetCDF data store and the R ncdf4 library (https://cran.r-project.org/web/packages/ncdf4/index.html) but could be easily adapted to use simple R data objects. For loading large datasets into NetCDF files I have found the Perl NetCDF libraries (https://metacpan.org/source/DHUNT/PDL-NetCDF-4.20/netcdf.pd) faster.
Note the bigPCA package (http://cran.nexr.com/web/packages/bigpca/index.html) has not been updated and requires an older version of R e.g. R-3.4.2 - other R PCA packages could be used. 

