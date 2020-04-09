# CamCNV
CamCNV pipeline for calling rare CNVs from Illumina data.

The basic steps of the pipeline as outlined in the example code are:
1. Load Log R Ratio intensities (LRR) and B Allele Frequencies (BAF) from Illumina gentotyping into a NetCDF data store
2. Run a principal component adjustment on the LRR
3. Calculate the mean and standard deviation of the LRR for each probe across all samples and convert the LRR into z-scores
4. For each sample segment the z-scores using DNACopy (https://bioconductor.org/packages/release/bioc/html/DNAcopy.html) and identify potential CNVs from the mean z-scores of the segments
5. Generate additional QC scores for each CNVs

The pipeline uses a NetCDF data store and the R ncdf4 library (https://cran.r-project.org/web/packages/ncdf4/index.html) but could be easily adapted to use simple R data objects. For loading large datasets into NetCDF files I have found the Perl NetCDF libraries (https://metacpan.org/source/DHUNT/PDL-NetCDF-4.20/netcdf.pd) faster.
