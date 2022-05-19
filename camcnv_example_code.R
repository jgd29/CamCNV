### Code Snippets for CamCNV calling pipeline
# 1. Create NetCDF datastore for a chromosome
# 2. Load Log R Ratio (LRR) intensities and B Allele Frequencies (BAF) 
# 3. Generate PCA-adjusted LRR
# 4. Calculate z-scores based on mean and sd of adjusted LRR at each probe
# 5. Segment z-scores and select segments that are potential CNVs; Count number of segments per sample
# 6. Generate additional QC scores for each CNV


library(ncdf4)
library(foreach)
library(DBI)
library(bigpca)
library(matrixStats)
library(DNAcopy)

##### STEP 1 #### 
# 1. Create a NetCDF datastore to hold the data 
#####

#Define the dimensions of the matrix - (number of samples X number of SNPs)
snpcount <- 533631
sample_count <- 24 
snpDim <- ncdim_def("snpDim", "SNP IDs", 1:snpcount) 
sampleDim <- ncdim_def("sampleDim", "Sample IDs",1:sample_count) 

#Sample variables
varOncID <- ncvar_def("sampleID", "person", sampleDim, -1, prec = "integer")
#Add optional variables for samples e.g. study, case control status etc.
varStudyID <- ncvar_def("studyID", "person", sampleDim, -1, prec = "integer")

#Variable for SNPs - SNP ID = the index position of the probes sorted by chromosome, position 
varSNPID <- ncvar_def("SNP_ID", "SNP_Index", snpDim, -1, longname="Unique_identifier_for_SNP", prec="integer")
varPos <- ncvar_def("SNP_Pos", "bases", snpDim, -1, longname="Build37_positions", prec="integer")
varChr <- ncvar_def("Chrom", "chrom", snpDim, -1, longname="chromosome", prec="integer")

#Variables at each SNP-Sample coordinate 
varGT <- ncvar_def("GT", "genotype", dim=list(snpDim,sampleDim), -1, longname="Genotype", prec="integer",compression=7)
varBAF <- ncvar_def("BAF", "B_allele", dim=list(snpDim,sampleDim), -1, longname="B_Allele_Freq", prec="single",compression=7)
varLRR <- ncvar_def("LRR_orig", "log_R", dim=list(snpDim,sampleDim), -1, longname="LogR_Ratio", prec="single",compression=7)
varLRR_adj <- ncvar_def("LRR_adj", "log_R", dim=list(snpDim,sampleDim), -1, longname="LogR_Ratio_Adjusted", prec="single",compression=7)
varZ <- ncvar_def("Z", "z_score", dim=list(snpDim,sampleDim), -1, longname="Z_score", prec="single",compression=7)
varCN <- ncvar_def("CN", "copy_number", dim=list(snpDim,sampleDim), -1, longname="Assigned_CN", prec="integer",compression=7)

#Create the empty NetCDF file with these variables
netcdf_file = "camcnv_onco.nc"
nc1 <- nc_create(netcdf_file,list(varOncID,varStudyID,varSNPID,varChr,varPos,varGT,varBAF,varLRR,varLRR_adj,varZ,varCN), force_v4=TRUE, verbose=TRUE )

#Load a file with positions and indices of the SNPs sorted by chromosome, position
snp_info <- read.csv("cnv_probes_ordered_by_pos.csv",header=TRUE)
ncvar_put( nc1, varSNPID, snp_info$index_by_pos, start=1, count=snpcount )
ncvar_put( nc1, varChr, snp_info$Chr_numeric, start=1, count=snpcount )
ncvar_put( nc1, varPos, snp_info$position_b37, start=1, count=snpcount )

nc_close(nc1)


##### STEP 2 #### 
# 2. Load into datastore LRR and BAF for each sample at each probe  
#####

nc <- nc_open("camcnv_onco.nc",write=TRUE)
start_snp <- 1
start_sample <- 1
snp_length <- 533631
sample_length <- 24

#Save our sample IDs in the NetCDF file
sample_id_input <- read.csv("sample_name_lookup.csv")
ncvar_put(nc,"sampleID", sample_id_input$Sample_ID)

#Create theB Allele Frequency and LRR input files by e.g. by exporting from Genome Studio Full Data Table
#Ensure that the data is sorted by chromosome and position
#Ensure that the columns in the files match the Sample ID in the order we have just saved into the NetCDF file by e.g. sorting by the Genome Studio Index column
baf_input <- read.csv("baf_data_table.txt",sep="\t")
#Take out header columns
baf_input <- as.matrix(baf_input[,4:27])
#Save into the datastore
ncvar_put(nc,"BAF", baf_input, c(start_snp,start_sample), c(snp_length,sample_length) )
#In this case we are saving all samples and all SNPs so could just say:
#ncvar_put(nc,"BAF", baf_input, c(1,1), c(-1,-1) )
rm(baf_input)

#Repeat process for the LRR
lrr_input <- read.csv("lrr_data_table.txt",sep="\t")
lrr_input <- as.matrix(lrr_input[,4:27])
ncvar_put(nc,"LRR_orig", lrr_input, c(1,1), c(-1,-1) )
rm(lrr_input)

#Repeat process for the genotypes
#Convert AA=0 AB=1 BB=2 NC=9
gt_input <- read.csv("genotype_data_table.txt",sep="\t")
gt_input <- as.matrix(gt_input[,4:27])
gt_input <- gsub("AA","0",gt_input)
gt_input <- gsub("AB","1",gt_input)
gt_input <- gsub("BB","2",gt_input)
gt_input <- gsub("NC","9",gt_input)
ncvar_put(nc,"GT", gt_input, c(1,1), c(-1,-1) )
rm(gt_input)


#### Optional step: Calculate DLRS to measure noise of samples
#Calculates the probe-to-probe log ratio difference of an array. This is a noise estimation which is robust to outliers. 
#See https://rdrr.io/cran/ADM3/man/dLRs.html
#Package is no longer available on CRAN so download from https://cran.r-project.org/src/contrib/Archive/ADM3/ and install from source:
#install.packages("ADM3_1.3.tar.gz", repos = NULL, type="source")
library(ADM3)

lrr_orig <- ncvar_get(nc,"LRR_orig",c(1,1),c(-1,-1))
sample_ids <- ncvar_get(nc,"sampleID")
#Bind IDs as column names
colnames(lrr_orig) <- sample_ids

dlrs <- apply(lrr_orig,2, dLRs)
write.table(dlrs,file="dlrs_before_pca.txt",quote=FALSE,col.names=FALSE,sep="\t")
#Plot DLRS to check outliers 
 filename <-  "dlrs_before_pca.pdf"
 pdf(filename)
 plot(dlrs,main="DLRS before PCA")
 dev.off()
 #Select noisy samples with high DLRS above some cutoff e.g 3 SD above mean
cut_off <-  mean(dlrs) + (3 * sd(dlrs))
exclude <- subset(dlrs,dlrs > cut_off)
exclude_dlrs_ids <- names(exclude)

##### STEP 3 #### 
##### Run Principal components adjustment (PCA) on LRR
#####

#Before PCA exclude list of SNPs that fail genotyping QC, telomeric, in MHC etc.
snp_ids <- ncvar_get(nc,"SNP_ID")
rownames(lrr_orig) <- snp_ids
exclude_snps <- read.csv("snps_to_exclude_before_pca.csv")
lrr_orig <- lrr_orig[!(rownames(lrr_orig) %in% exclude_snps$index_by_pos),]

#load into big matrix
bm_lrr <- as.big.matrix(lrr_orig)
rm(lrr_orig)

#Calculate PCs - thin to 15% keeping most informative
pc_lrr  <- big.PCA(bm_lrr,verbose=TRUE,pcs.to.keep=30,thin(bm_lrr,.15,"PCA",))

 #Print scree plot to working directory
 filename <- "pc_scree.pdf"
 pdf(filename)
 pca.scree.plot(pc_lrr$Evalues,min.dim=20,n.xax=30,verbose=TRUE)
dev.off()

#From scree plot select how many PCs are above the 'elbow' 
pc_number <- 3
rm(bm_lrr)

#Now reload complete set of LRR and adjust for n PCs
lrr_to_correct <- ncvar_get(nc,"LRR_orig",c(1,1),c(-1,-1))
lrr_corrected2 <- PC.correct(pc_lrr,lrr_to_correct,num.pcs=pc_number,pref="camcnv_example")
lrr_corr <- get.big.matrix(lrr_corrected2,verbose=TRUE)
lrr_new <- as.matrix(lrr_corr)
#Save into the NetCDF file
ncvar_put(nc,"LRR_adj", lrr_new, c(1,1), c(-1,-1) )

#Re-run DLRS to check noise has reduced and select samples who are still noisy to exclude
dlrs_adj <- apply(lrr_corr,2, dLRs)
write.table(dlrs_adj,file="dlrs_after_pca.txt",quote=FALSE,col.names=FALSE,sep="\t")

filename <- "dlrs_after_pca.pdf"
pdf(filename)
plot(dlrs_adj,main="DLRS after PCA")
dev.off()
cut_off <-  mean(dlrs_adj) + (3 * sd(dlrs_adj))
exclude <- subset(dlrs_adj,dlrs_adj > cut_off)
exclude_dlrs_ids2 <- names(exclude)

#Tidy up - .bck files can also be deleted 
rm(lrr_corr)
rm(lrr_to_correct)
rm(lrr_new)
rm(pc_lrr)

##### STEP 4 #### 
##### Calculate z-scores for the adjusted LRR
#####

#Get mean and standard deviation for each SNP
lrr_adj <- ncvar_get(nc,"LRR_adj",c(start_snp,start_sample), c(snp_length,sample_length))
snp_ids <- ncvar_get(nc,"SNP_ID")
rownames(lrr_adj) <- snp_ids
lrr_means <- rowMeans(lrr_adj, na.rm = TRUE)
lrr_sds <- rowSds(lrr_adj, na.rm = TRUE)
#Inspect the SDs for all probes and select some cut-off to exclude the noisiest probes from the segmentation process 

#For the example z-score calculation we will load the SNP means and sds calculated from a large set of samples
means_file <- read.csv("means_sds_onco.csv")
lrr_means <- means_file$lrr_means
lrr_sds <- means_file$lrr_sds

#Save z-scores
z1 <- (lrr_adj - lrr_means)/lrr_sds
ncvar_put( nc,"Z", z1,  c(start_snp,start_sample), c(snp_length,sample_length))

##### STEP 5 #### 
##### Circular binary segmentation on z-scores with DNACopy
#####

snp_id <- ncvar_get(nc,"SNP_ID")
chr <- ncvar_get(nc,"Chrom")
snp_pos <- ncvar_get(nc,"SNP_Pos")
z1 <- ncvar_get( nc,"Z",c(start_snp,start_sample), c(snp_length,sample_length))
sample_ids <- ncvar_get(nc,"sampleID")
colnames(z1) <- sample_ids

z1_df <- as.data.frame(cbind(z1,snp_id,chr,snp_pos))
#Load list of SNPs to exclude from the segmentation process:
#e.g. Noisy probes ( LRR SD > some cutoff), probes in common CNV regions, SNPs that failed to cluster, low intensity, probes too close together
exclude_probes <- read.csv("cnv_probe_exclusions.csv")
z1_df <- z1_df[!(z1_df$snp_id %in% exclude_probes$index_by_pos),]
#Remove Y and MT probes
z1_df <- z1_df[z1_df$chr <24,]
 
#Detect segments for each sample
for  (i in 1:sample_length) {
print(i)
CNA.object <- CNA(cbind(z1_df[i]),z1_df$chr,z1_df$snp_pos, data.type="logratio",sampleid=colnames(z1_df[i]),presorted = TRUE)
#Smooth single point outliers
smoothed.CNA.object <- smooth.CNA(CNA.object)
these_segments <- segment(smoothed.CNA.object,undo.splits="sdundo", undo.SD=2, verbose=0)
	 if(i == 1){
		output_comb <-  these_segments$output
	 }else{
		output_comb <- rbind(output_comb,these_segments$output)
	 }
}

#Create ID column without X at start
output_comb$SampleID <- as.numeric(gsub("X","",output_comb$ID))
write.csv(output_comb, "dnacopy_segments.csv")

#Use number of probes (num.mark) and mean z-score (seg.mean) to filter the list of segments down to possible CNVs 
cnvs <- output_comb[output_comb$num.mark >2 & output_comb$num.mark <200,]
cnvs$type <- ""
cnvs[cnvs$seg.mean < 0,"type"] <- "deletion"
cnvs[cnvs$seg.mean > 0,"type"] <- "duplication"
#Apply mean z-score cut-offs 
cnvs <- cnvs[(cnvs$seg.mean < -3.7 & cnvs$type == "deletion") | (cnvs$seg.mean > 2 & cnvs$type == "duplication"),]
#And upper bounds
cnvs <- cnvs[(cnvs$seg.mean > -14 & cnvs$type == "deletion") | (cnvs$seg.mean < 10 & cnvs$type == "duplication"),]

write.csv( cnvs, "called_cnvs_v1.csv")

#### For each sample count the total number of short segments  
# We will exclude noisy samples with an excessive number of segments
segs <- output_comb[output_comb$num.mark >2 & output_comb$num.mark <200,]
seg_counts <- as.data.frame(table(segs$SampleID))
#In a large sample set it is likely that segment counts distribution will have a heavy tail with many samples with excessive segments so transform data and use median for an exclusion cut-off
seg_counts$y <- sqrt(seg_counts$Freq) * 2 
cut_off <- median(seg_counts$y) + 3.5
#Transform back to get maximum number of segments
seg_cut_off <- (cut_off/2) * (cut_off/2)
write.csv( seg_counts, "segment_counts.csv")

##### STEP 6 #### 
##### Generate additional QC scores for each CNV using:
##### a. the shift in original LRR between the LRR for probes within the CNV and the mean LRR for this sample
##### b. the B Allele Frequency measurements within the CNVs
#####

lrr <- ncvar_get(nc,"LRR_orig",c(start_snp,start_sample), c(snp_length,sample_length))
baf <- ncvar_get(nc,"BAF",c(start_snp,start_sample), c(snp_length,sample_length))
rownames(baf) <- snp_id
gt <- ncvar_get(nc,"GT",c(start_snp,start_sample), c(snp_length,sample_length))
rownames(gt) <- snp_id

lrr_df <- as.data.frame(cbind(snp_id,chr,snp_pos,lrr))
baf_df <- as.data.frame(cbind(snp_id,chr,snp_pos,baf))
#Filter out excluded SNPs
lrr_df <- lrr_df[!(lrr_df$snp_id %in% exclude_probes$index_by_pos),]
baf_df <- baf_df[!(baf_df$snp_id %in% exclude_probes$index_by_pos),]


for  (i in 1:sample_length) {
	this_sample <- sample_ids[i]
	print(this_sample)
	this_gt <- gt[,i]
	this_baf <- baf[,i]
	include_snps <- names(this_gt[this_gt==1])
	#First get the mean and sd of the  LRR for each sample for each chromosome 
	#Select the heterozygous SNPs within the normal BAF range so that we are confident we are measuring two-copy probes 
	include_snps2 <- names(this_baf[this_baf > 0.45 & this_baf < 0.55])
	this_lrr <- lrr_df[lrr_df$snp_id %in% include_snps & lrr_df$snp_id %in% include_snps2,]
    #Filter to this sample
	this_lrr <- this_lrr[,c(2,i+3)]
	colnames(this_lrr) <- c("chr","lrr")
	this_mean <- aggregate( lrr ~ chr  ,this_lrr, mean )
	this_sd <-  aggregate( lrr ~ chr  ,this_lrr, sd )
	
	#Select LRR and BAF for all QC probes for this sample
	sample_lrr <- lrr_df[,c(1,2,3,i+3)]
	sample_baf <- baf_df[,c(1,2,3,i+3)]

	#Retrieve the CNVs for this sample
	these_cnvs <- cnvs[cnvs$SampleID == this_sample,]
	#Score each CNV for the shift between the CNV LRR and the sample's mean for the chromosome
	cnv_count <- dim(these_cnvs)[1]
	if(cnv_count > 0 ){
    for(j in 1:cnv_count){
	cnv_mean <- mean(sample_lrr[sample_lrr$chr == these_cnvs[j,"chrom"] & sample_lrr$snp_pos >= these_cnvs[j,"loc.start"] & sample_lrr$snp_pos <= these_cnvs[j,"loc.end"],4])
	chr_mean <- this_mean[this_mean$chr == these_cnvs[j,"chrom"],"lrr"]
	
		if(these_cnvs[j,"type"] == "deletion"){
			lrr_diff <-  chr_mean - cnv_mean 
		}else{
			lrr_diff <-  cnv_mean -  chr_mean
		}
	cnv_id <- rownames(these_cnvs[j,])
	cnvs[rownames(cnvs) == cnv_id,"lrr_shift_score"] <- lrr_diff
	
	cnv_baf <- sample_baf[sample_baf$chr == these_cnvs[j,"chrom"] & sample_baf$snp_pos >= these_cnvs[j,"loc.start"] & sample_baf$snp_pos <= these_cnvs[j,"loc.end"],4]
	probe_count <- length(cnv_baf)
			
		if(these_cnvs[j,"type"] == "deletion"){
		#For deletions BAF should be close to zero or one 
			bad_probes <- sum(cnv_baf > 0.1 & cnv_baf < 0.9)
			baf_problem_pc <-  bad_probes / probe_count 
			cnvs[rownames(cnvs) == cnv_id,"baf_negative_score"] <- baf_problem_pc
		}else{
			#For duplications look for probes with a BAF that looks like AAB or ABB
			dup_probes <- sum(cnv_baf > 0.2 & cnv_baf < 0.4) + sum(cnv_baf > 0.6 & cnv_baf < 0.8)
			baf_good_pc <-  dup_probes / probe_count 
			cnvs[rownames(cnvs) == cnv_id,"baf_positive_score"] <- baf_good_pc
		}
	
	}
	}
}

#Exclude the CNVS outside cut-offs
cnvs[cnvs$lrr_shift_score > 0.8 & cnvs$type == "deletion" ,"QC_exclusion"] <- "LRR_Shift"
cnvs[cnvs$lrr_shift_score < 0.2 & cnvs$type == "deletion" ,"QC_exclusion"] <- "LRR_Shift"
cnvs[cnvs$lrr_shift_score > 0.4 & cnvs$type == "duplication" ,"QC_exclusion"] <- "LRR_Shift"
cnvs[cnvs$lrr_shift_score < 0.1 & cnvs$type == "duplication" ,"QC_exclusion"] <- "LRR_Shift"

cnvs$SampleID <- as.numeric(gsub("X","",cnvs$ID))
#Write CNVs with potential exclusions
write.csv( cnvs, "called_cnvs_v2.csv")

# Note: using this very small set of samples for the PCA is likely to produce unreliable calls but it will detect some genuine CNVs e.g. 
#ID	chrom	loc.start	loc.end	num.mark	seg.mean
#170	3	60831012	60908013	9	-7.4408
#170	20	55408853	55421778	12	-5.606
#for sample NA11992 which we can see from Ensembl carries the variant esv3646204 at 20:55406186-55428243 and esv3596323 at 3:60826517-60918647
