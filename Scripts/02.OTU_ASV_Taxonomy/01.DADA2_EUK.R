## DADA2 for creating ASV tables from trimmed sequences

##############This script is executed in Astbury HPC#################
### HPC execution (in bash)
##Create a single directory for all trimmed files

mkdir all_run1_trimmed
mkdir all_run2_trimmed

#copy all run1 and run2 trimmed files into respective dir

cp GAB_EUK*/run1/*trim*.fq.gz all_run1_trimmed
cp GAB_EUK*/run2/*trim*.fq.gz all_run2_trimmed

###############################################################################
###############################################################################
############################ Packages #########################################



##need to be loaded at every session in HPC
library(devtools) #in conda MetabMaple env

#remotes::install_github("benjjneb/dada2", ref = "master") #Sample inference (ASV tables)
library(dada2)
library(ShortRead)

#install_github("tobiasgf/lulu") #clustering curation of ASV tables
library(lulu)

library(ggplot2)
library(dplyr)


###############################################################################
########################## PATH CONFIGURATION #################################
###############################################################################

BASE_DIR        <- "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001"
DATA_DIR        <- "/data/glennsdata/MAPLE/18S/EUK_data"

RAW_DATA_DIR      <- file.path(BASE_DIR, "01.RawData")
PATH_RUN1         <- file.path(RAW_DATA_DIR, "all_run1_trimmed")
PATH_RUN2         <- file.path(RAW_DATA_DIR, "all_run2_trimmed")
RWORKSPACE_DIR    <- file.path(DATA_DIR, "R_workspaces")
READ_TRACKING_DIR <- file.path(RWORKSPACE_DIR, "ReadTracking")
LULU_DIR          <- file.path(RAW_DATA_DIR, "LULU_fasta")
SWARM_DIR         <- file.path(RAW_DATA_DIR, "swarm")
NCBI_DIR          <- file.path(RAW_DATA_DIR, "ncbi_data")
TAXA_DB_DIR       <- file.path(RAW_DATA_DIR, "taxa_db")
SILVA_DB_PATH     <- file.path(TAXA_DB_DIR, "silva_132.18s.99_rep_set.dada2.fa.gz")

## Create directories if they don't exist
invisible(lapply(
  c(PATH_RUN1, PATH_RUN2, RWORKSPACE_DIR, READ_TRACKING_DIR, LULU_DIR, SWARM_DIR, NCBI_DIR, TAXA_DB_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

###############################################################################


### Tracking system
read_tracking = list()

read_tracking$metadata = list(
  date_created = Sys.time(),
  pipeline_version = "DADA2_EUK18S",
  runs= c("run1", "run2"),
  pooling_methods = c("nP", "psP", "fP")
)

readTracking1 = readRDS(file.path(READ_TRACKING_DIR, "1.input_readTrack.rds"))
readTracking2 = readRDS(file.path(READ_TRACKING_DIR, "2.filtered_readTrack.rds"))
readTracking3 = readRDS(file.path(READ_TRACKING_DIR, "3.denoised_readTrack.rds"))
readTracking4 = readRDS(file.path(READ_TRACKING_DIR, "4.merged_readTrack.rds"))
readTracking5 = readRDS(file.path(READ_TRACKING_DIR, "5.nochimera_readTrack.rds"))
readTracking6 = readRDS(file.path(READ_TRACKING_DIR, "6.merged_runs_readTrack.rds"))
readTracking7 = readRDS(file.path(READ_TRACKING_DIR, "7.lulu_readTrack.rds"))
readTracking8 = readRDS(file.path(READ_TRACKING_DIR, "8.swarm_readTrack_FINAL.rds"))




###############################################################################
########################### Directories ########################################

#Set directory where trimmed files are (in HPC)

setwd(RAW_DATA_DIR)
path_run1 = PATH_RUN1
path_run2 = PATH_RUN2


###############################################################################
############################### Filenames #####################################
# create separate lists for forward and reverse reads
#run1
fnFs.run1 = sort(list.files(path_run1, pattern = "trim1.fq.gz$", full.names = TRUE))
fnRs.run1 = sort(list.files(path_run1, pattern = "trim2.fq.gz$", full.names = TRUE))
#run2
fnFs.run2 = sort(list.files(path_run2, pattern = "trim1.fq.gz$", full.names = TRUE))
fnRs.run2 = sort(list.files(path_run2, pattern = "trim2.fq.gz$", full.names = TRUE))

#Function to remove trimming suffix
get_sample_name = function(fname) {
  sub("\\.trim[12]\\.fq\\.gz$", "", basename(fname))
}

#Extract sample names
#run1
sample.names.F_run1 = sapply(fnFs.run1, get_sample_name)
sample.names.R_run1 = sapply(fnRs.run1, get_sample_name)
#run2
sample.names.F_run2 = sapply(fnFs.run2, get_sample_name)
sample.names.R_run2 = sapply(fnRs.run2, get_sample_name)

############ Robust Quality control (not in dada2 tutorial) ##################

#Find samples with matching fwd and rv
common.samples_run1 = intersect(sample.names.F_run1, sample.names.R_run1)
common.samples_run2 = intersect(sample.names.F_run2, sample.names.R_run2)


#Create vectors for easy matching (subset by name rather than position)
#run1
fnFs.named_run1 = setNames(fnFs.run1, sample.names.F_run1)
fnRs.named_run1 = setNames(fnRs.run1, sample.names.R_run1)
#run2
fnFs.named_run2 = setNames(fnFs.run2, sample.names.F_run2)
fnRs.named_run2 = setNames(fnRs.run2, sample.names.R_run2)

#Keep only matched pairs
#run1
fnFs_r1 = fnFs.named_run1[common.samples_run1]
fnRs_r1 = fnRs.named_run1[common.samples_run1]
r1_sample.names = names(fnFs_r1) # Definitive names list for run1
#run2
fnFs_r2 = fnFs.named_run2[common.samples_run2]
fnRs_r2 = fnRs.named_run2[common.samples_run2]
r2_sample.names = names(fnFs_r2)

####Final safety check
#run1
stopifnot(length(fnFs_r1) == length(fnRs_r1))
stopifnot(all(names(fnFs_r1) == names(fnRs_r1)))
#run2
stopifnot(length(fnFs_r2) == length(fnRs_r2))
stopifnot(all(names(fnFs_r2) == names(fnRs_r2)))

#For peace of mind: confirm first few matches
head(data.frame(forward = basename(fnFs_r1), reverse = basename(fnRs_r1), sample = r1_sample.names))



######################### Read track check up ##################################


# Function to count reads in FASTQ files
count_fastq_reads <- function(file_path) {
  if (!file.exists(file_path)) return(0)
  # For .gz files, use gzfile
  if (grepl("\\.gz$", file_path)) {
    lines <- length(readLines(gzfile(file_path)))
  } else {
    lines <- length(readLines(file_path))
  }
  return(lines / 4)  # FASTQ has 4 lines per read
}

# Count reads in input files (trim1 files from cutadapt)
cat("Counting reads in input files...\n")

input_counts_r1 <- sapply(fnFs_r1, count_fastq_reads)
input_counts_r2 <- sapply(fnFs_r2, count_fastq_reads)

# Store in tracking object
read_tracking$step1_input <- list(
  run1 = data.frame(
    sample = r1_sample.names,
    reads = input_counts_r1,
    stringsAsFactors = FALSE
  ),
  run2 = data.frame(
    sample = r2_sample.names,
    reads = input_counts_r2,
    stringsAsFactors = FALSE
  )
)

# Calculate totals
read_tracking$step1_input$run1_total <- sum(input_counts_r1)
read_tracking$step1_input$run2_total <- sum(input_counts_r2)
read_tracking$step1_input$combined_total <- sum(input_counts_r1) + sum(input_counts_r2)

cat("Total input reads (run1):", read_tracking$step1_input$run1_total, "\n")
cat("Total input reads (run2):", read_tracking$step1_input$run2_total, "\n")
cat("Total input reads (combined):", read_tracking$step1_input$combined_total, "\n\n")

# Save intermediate tracking
saveRDS(read_tracking, file = file.path(READ_TRACKING_DIR, "1.input_readTrack.rds"))



###################### Quality Inspection #####################################
#Using only the forward reads we can visualize the quality profiles
# were gonna randomly sample across the dataset


library(ggplot2)

# plots are saves as pnj files
png("EUK_qualityProfiles_FWD_r1.png")
plotQualityProfile(fnFs_r1[1:4])
dev.off()
png("EUK_qualityProfiles_RV_r1.png")
plotQualityProfile(fnRs_r1[1:4])
dev.off()

png("EUK_qualityProfiles_FWD_r2.png")
plotQualityProfile(fnFs_r2[1:4])
dev.off()
png("EUK_qualityProfiles_RV_r2.png")
plotQualityProfile(fnRs_r2[1:4])
dev.off()

###
### quality seems good for that bunch of samples
###


###############################################################################
###################### Filtering and Trimming #################################

#Setting file paths for every single samples
filtFs_r1 = file.path(path_run1, "filtered", paste0(r1_sample.names, "_F_filt.fastq.gz"))
filtRs_r1 = file.path(path_run1, "filtered", paste0(r1_sample.names, "_R_filt.fastq.gz"))

filtFs_r2 = file.path(path_run2, "filtered", paste0(r2_sample.names, "_F_filt.fastq.gz"))
filtRs_r2 = file.path(path_run2, "filtered", paste0(r2_sample.names, "_R_filt.fastq.gz"))

#With short amplicons (PE150), no need for trimming, only filtering

out_run1 = filterAndTrim(
  fnFs_r1, #directory with the fastq files
  filtFs_r1, #output directory
  fnRs_r1,
  filtRs_r1,
  truncQ = 2, #truncate reads if quality score is less or equal to 2
  compress = TRUE, #output fastq files are gzipped
  maxN = 0, #reads with more than 0 ambiguous bases
  maxEE = c(2,2), #reads with more than 2 expected errors will be discarded
  rm.phix = TRUE, #discard reads that match agains phiX genome
  multithread = TRUE #parallel processing, ok in HPC, not on windows
)
#save output
#saveRDS(out_run1, file = file.path(RWORKSPACE_DIR, "out_run1.rds"))
#load output
out_run1 = readRDS(file.path(RWORKSPACE_DIR, "out_run1.rds"))



out_run2 = filterAndTrim(
  fnFs_r2, #directory with the fastq files
  filtFs_r2, #output directory
  fnRs_r2,
  filtRs_r2,
  truncQ = 2, #truncate reads if quality score is less or equal to 2
  compress = TRUE, #output fastq files are gzipped
  maxN = 0, #reads with more than 0 ambiguous bases
  maxEE = c(2,2), #reads with more than 2 expected errors will be discarded
  rm.phix = TRUE, #discard reads that match againts phiX genome
  multithread = TRUE #parallel processing, ok in HPC, not on windows
)
#save output
#saveRDS(out_run2, file = file.path(RWORKSPACE_DIR, "out_run2.rds"))
#load output
out_run2 = readRDS(file.path(RWORKSPACE_DIR, "out_run2.rds"))


######## Finding samples and removing that did not pass the filter #############


##DF with existing and non-existing files
#rows = sample, cols= TRUE/FALSE and path
df.fe_r1 = data.frame(theref = file.exists(filtFs_r1), #checks file path vector, returns TRUE if file exists
                      therer = file.exists(filtRs_r1), #same for reverse files
                      filef = filtFs_r1, #stores actual paths for forward
                      filer = filtRs_r1) #stores actual paths for reverse

df.fe_r2 = data.frame(theref = file.exists(filtFs_r2), therer = file.exists(filtRs_r2), filef = filtFs_r2, filer = filtRs_r2)

##Finding avec create vector that contains samples that have failed filtering
drops_r1 = as.numeric(rownames(subset(df.fe_r1, theref =="FALSE"))) #length = 452
drops_r2 = as.numeric(rownames(subset(df.fe_r2, theref =="FALSE"))) #length = 1160

##Cleaning up filtering matrix
#Now we only have rows for samples that successfully passes the filtering step
out_run1 = out_run1[file.exists(filtFs_r1),]
out_run2 = out_run2[file.exists(filtFs_r2),]

##Cleaning up all file vectors
#All vectors should then be the same length and only contain filtered samples
#run1
filtFs_r1 = filtFs_r1[-drops_r1]
filtRs_r1 = filtRs_r1[-drops_r1]
r1_sample.names = r1_sample.names[-drops_r1]
#run2
filtFs_r2 = filtFs_r2[-drops_r2]
filtRs_r2 = filtRs_r2[-drops_r2]
r2_sample.names = r2_sample.names[-drops_r2]

save.image(file = file.path(RWORKSPACE_DIR, "FiltTrim_dada2_EUK.RData"))

################# Read tracking check point ################################

# "Create tracking dataframe directly from out_run1"

filtered_counts_r1 <- data.frame(
  sample = rownames(out_run1),
  input = out_run1[, "reads.in"],
  filtered = out_run1[, "reads.out"],
  stringsAsFactors = FALSE
)

# Calculate losses
filtered_counts_r1$reads_lost <- filtered_counts_r1$input - filtered_counts_r1$filtered
filtered_counts_r1$percent_retained <- (filtered_counts_r1$filtered / filtered_counts_r1$input) * 100

# Repeat for run2
filtered_counts_r2 <- data.frame(
  sample = rownames(out_run2),
  input = out_run2[, "reads.in"],
  filtered = out_run2[, "reads.out"],
  stringsAsFactors = FALSE
)

filtered_counts_r2$reads_lost <- filtered_counts_r2$input - filtered_counts_r2$filtered
filtered_counts_r2$percent_retained <- (filtered_counts_r2$filtered / filtered_counts_r2$input) * 100


# Store in tracking object
read_tracking$step2_filtered <- list(
  run1 = filtered_counts_r1,
  run2 = filtered_counts_r2,
  run1_samples_failed = length(drops_r1),
  run2_samples_failed = length(drops_r2),
  run1_samples_passed = nrow(out_run1),
  run2_samples_passed = nrow(out_run2)
)

# Calculate totals
read_tracking$step2_filtered$run1_total_filtered <- sum(filtered_counts_r1$filtered)
read_tracking$step2_filtered$run2_total_filtered <- sum(filtered_counts_r2$filtered)
read_tracking$step2_filtered$combined_total_filtered <- 
  sum(filtered_counts_r1$filtered) + sum(filtered_counts_r2$filtered)

# Calculate overall retention
read_tracking$step2_filtered$run1_percent_retained <- 
  (read_tracking$step2_filtered$run1_total_filtered / read_tracking$step1_input$run1_total) * 100
read_tracking$step2_filtered$run2_percent_retained <- 
  (read_tracking$step2_filtered$run2_total_filtered / read_tracking$step1_input$run2_total) * 100
read_tracking$step2_filtered$combined_percent_retained <- 
  (read_tracking$step2_filtered$combined_total_filtered / read_tracking$step1_input$combined_total) * 100

# Print summary
cat("\n--- Filtering Summary ---\n")
cat("Run1:\n")
cat("  Samples that passed filtering:", read_tracking$step2_filtered$run1_samples_passed, "\n")
cat("  Samples that failed filtering:", read_tracking$step2_filtered$run1_samples_failed, "\n")
cat("  Total reads after filtering:", format(read_tracking$step2_filtered$run1_total_filtered, big.mark = ","), "\n")
cat("  Retention rate:", round(read_tracking$step2_filtered$run1_percent_retained, 2), "%\n\n")

cat("Run2:\n")
cat("  Samples that passed filtering:", read_tracking$step2_filtered$run2_samples_passed, "\n")
cat("  Samples that failed filtering:", read_tracking$step2_filtered$run2_samples_failed, "\n")
cat("  Total reads after filtering:", format(read_tracking$step2_filtered$run2_total_filtered, big.mark = ","), "\n")
cat("  Retention rate:", round(read_tracking$step2_filtered$run2_percent_retained, 2), "%\n\n")

cat("Combined:\n")
cat("  Total reads after filtering:", format(read_tracking$step2_filtered$combined_total_filtered, big.mark = ","), "\n")
cat("  Overall retention rate:", round(read_tracking$step2_filtered$combined_percent_retained, 2), "%\n\n")

# Save intermediate tracking
saveRDS(read_tracking, file = file.path(READ_TRACKING_DIR, "2.filtered_readTrack.rds"))



###############################################################################
###################### Learning Error Rates ##################################

#https://github.com/benjjneb/dada2/issues/1307
#install dads2 dev version with remotes::

## Finding the binned Q-scores specific to our data

fq_sample_r1 <- readFastq(fnFs_r1)
quality_scores_r1 <- as.vector(as(quality(fq_sample_r1), "matrix"))
table(quality_scores_r1)
#2-9-24-40

fq_sample_r2 <- readFastq(fnFs_r2)
quality_scores_r2 <- as.vector(as(quality(fq_sample_r2), "matrix"))
table(quality_scores_r2)
#2-9-24-10



binnedQs = c(2,9,24,40)
binnedQualErrfun = makeBinnedQualErrfun(binnedQs) # custom error estimation model

#Learn error rates with custom function

#run1
errF_r1 = learnErrors(filtFs_r1, errorEstimationFunction = binnedQualErrfun, nbases = 2e10, multi = TRUE)
saveRDS(errF_r1, file = file.path(RWORKSPACE_DIR, "errF_r1.rds"))
#nbases = 20 billion bases for learning
errR_r1 = learnErrors(filtRs_r1, errorEstimationFunction = binnedQualErrfun, nbases = 2e10, multi = TRUE)
saveRDS(errR_r1, file = file.path(RWORKSPACE_DIR, "errR_r1.rds"))

#run2
errF_r2 = learnErrors(filtFs_r2, errorEstimationFunction = binnedQualErrfun, nbases = 2e10, multi = TRUE)
saveRDS(errF_r2, file = file.path(RWORKSPACE_DIR, "errF_r2.rds"))
#nbases = 20 billion bases for learning
errR_r2 = learnErrors(filtRs_r2, errorEstimationFunction = binnedQualErrfun, nbases = 2e10, multi = TRUE)
saveRDS(errR_r2, file = file.path(RWORKSPACE_DIR, "errR_r2.rds"))



errF_r1 = readRDS(file.path(RWORKSPACE_DIR, "errF_r1.rds"))
errR_r1 = readRDS(file.path(RWORKSPACE_DIR, "errR_r1.rds"))
errF_r2 = readRDS(file.path(RWORKSPACE_DIR, "errF_r2.rds"))
errR_r2 = readRDS(file.path(RWORKSPACE_DIR, "errR_r2.rds"))

################### Quality control
dada2:::checkConvergence(errF_r1)
dada2:::checkConvergence(errR_r1)
dada2:::checkConvergence(errF_r2)
dada2:::checkConvergence(errR_r2)
# Good convergence on all 4 models !!

################## Visualizing error learning
pdf("Fwd_r1_binnedError.pdf")
plotErrors(errF_r1, nominalQ=TRUE)
dev.off()

pdf("Rv_r1_binnedError.pdf")
plotErrors(errR_r1, nominalQ=TRUE)
dev.off()

pdf("Fwd_r2_binnedError.pdf")
plotErrors(errF_r2, nominalQ=TRUE)
dev.off()

pdf("Rv_r2_binnedError.pdf")
plotErrors(errR_r2, nominalQ=TRUE)
dev.off()

################################################################################
######################### Dereplication ########################################

derepFwd_r1 = derepFastq(filtFs_r1, verbose = TRUE)
derepRev_r1 = derepFastq(filtRs_r1, verbose = TRUE)

derepFwd_r2 = derepFastq(filtFs_r2, verbose = TRUE)
derepRev_r2 = derepFastq(filtRs_r2, verbose = TRUE)

#Derep objects have no names, so here we assigned them sample names
names(derepFwd_r1) = r1_sample.names
names(derepRev_r1) = r1_sample.names

names(derepFwd_r2) = r2_sample.names
names(derepRev_r2) = r2_sample.names

saveRDS(derepFwd_r1, file = file.path(RWORKSPACE_DIR, "derepFwd_r1.rds"))
saveRDS(derepRev_r1, file = file.path(RWORKSPACE_DIR, "derepRev_r1.rds"))

saveRDS(derepFwd_r2, file = file.path(RWORKSPACE_DIR, "derepFwd_r2.rds"))
saveRDS(derepRev_r2, file = file.path(RWORKSPACE_DIR, "derepRev_r2.rds"))

derepFwd_r1 = readRDS(file.path(RWORKSPACE_DIR, "derepFwd_r1.rds"))
derepRev_r1 = readRDS(file.path(RWORKSPACE_DIR, "derepRev_r1.rds"))

derepFwd_r2 = readRDS(file.path(RWORKSPACE_DIR, "derepFwd_r2.rds"))
derepRev_r2 = readRDS(file.path(RWORKSPACE_DIR, "derepRev_r2.rds"))

################################################################################
####################### Sample Inference #######################################

## No pooling

r1.dadaFs_nP = dada(derepFwd_r1, err = errF_r1, multithread = TRUE, verbose = TRUE)
r1.dadaRs_nP = dada(derepRev_r1, err = errR_r1, multithread = TRUE, verbose = TRUE)

r2.dadaFs_nP = dada(derepFwd_r2, err = errF_r2, multithread = TRUE, verbose = TRUE)
r2.dadaRs_nP = dada(derepRev_r2, err = errR_r2, multithread = TRUE, verbose = TRUE)

saveRDS(r1.dadaFs_nP, file = file.path(RWORKSPACE_DIR, "r1.dadaFs_nP.rds"))
saveRDS(r1.dadaRs_nP, file = file.path(RWORKSPACE_DIR, "r1.dadaRs_nP.rds"))

saveRDS(r2.dadaFs_nP, file = file.path(RWORKSPACE_DIR, "r2.dadaFs_nP.rds"))
saveRDS(r2.dadaRs_nP, file = file.path(RWORKSPACE_DIR, "r2.dadaRs_nP.rds"))


## Pseudo-Pooling

r1.dadaFs_psP = dada(derepFwd_r1, err = errF_r1, multithread = TRUE, verbose = TRUE, pool = "pseudo")
r1.dadaRs_psP = dada(derepRev_r1, err = errR_r1, multithread = TRUE, verbose = TRUE, pool = "pseudo")

r2.dadaFs_psP = dada(derepFwd_r2, err = errF_r2, multithread = TRUE, verbose = TRUE, pool = "pseudo")
r2.dadaRs_psP = dada(derepRev_r2, err = errR_r2, multithread = TRUE, verbose = TRUE, pool = "pseudo")

saveRDS(r1.dadaFs_psP, file = file.path(RWORKSPACE_DIR, "r1.dadaFs_psP.rds"))
saveRDS(r1.dadaRs_psP, file = file.path(RWORKSPACE_DIR, "r1.dadaRs_psP.rds"))

saveRDS(r2.dadaFs_psP, file = file.path(RWORKSPACE_DIR, "r2.dadaFs_psP.rds"))
saveRDS(r2.dadaRs_psP, file = file.path(RWORKSPACE_DIR, "r2.dadaRs_psP.rds"))


## Full pooling

r1.dadaFs_fP = dada(derepFwd_r1, err = errF_r1, multithread = TRUE, verbose = TRUE, pool = TRUE)
r1.dadaRs_fP = dada(derepRev_r1, err = errR_r1, multithread = TRUE, verbose = TRUE, pool = TRUE)

r2.dadaFs_fP = dada(derepFwd_r2, err = errF_r2, multithread = TRUE, verbose = TRUE, pool = TRUE)
r2.dadaRs_fP = dada(derepRev_r2, err = errR_r2, multithread = TRUE, verbose = TRUE, pool = TRUE)

saveRDS(r1.dadaFs_fP, file = file.path(RWORKSPACE_DIR, "r1.dadaFs_fP.rds"))
saveRDS(r1.dadaRs_fP, file = file.path(RWORKSPACE_DIR, "r1.dadaRs_fP.rds"))

saveRDS(r2.dadaFs_fP, file = file.path(RWORKSPACE_DIR, "r2.dadaFs_fP.rds"))
saveRDS(r2.dadaRs_fP, file = file.path(RWORKSPACE_DIR, "r2.dadaRs_fP.rds"))

#############
#no pooling
r1.dadaFs_nP = readRDS(file.path(RWORKSPACE_DIR, "r1.dadaFs_nP.rds"))
r1.dadaRs_nP = readRDS(file.path(RWORKSPACE_DIR, "r1.dadaRs_nP.rds"))
r2.dadaFs_nP = readRDS(file.path(RWORKSPACE_DIR, "r2.dadaFs_nP.rds"))
r2.dadaRs_nP = readRDS(file.path(RWORKSPACE_DIR, "r2.dadaRs_nP.rds"))

#pseudo-pooling
r1.dadaFs_psP = readRDS(file.path(RWORKSPACE_DIR, "r1.dadaFs_psP.rds"))
r1.dadaRs_psP = readRDS(file.path(RWORKSPACE_DIR, "r1.dadaRs_psP.rds"))
r2.dadaFs_psP = readRDS(file.path(RWORKSPACE_DIR, "r2.dadaFs_psP.rds"))
r2.dadaRs_psP = readRDS(file.path(RWORKSPACE_DIR, "r2.dadaRs_psP.rds"))


#Full pooling
r1.dadaFs_fP = readRDS(file.path(RWORKSPACE_DIR, "r1.dadaFs_fP.rds"))
r1.dadaRs_fP = readRDS(file.path(RWORKSPACE_DIR, "r1.dadaRs_fP.rds"))
r2.dadaFs_fP = readRDS(file.path(RWORKSPACE_DIR, "r2.dadaFs_fP.rds"))
r2.dadaRs_fP = readRDS(file.path(RWORKSPACE_DIR, "r2.dadaRs_fP.rds"))

save(list = c("derepFwd_r1", "derepRev_r1", "derepFwd_r2", "derepRev_r2", "r1.dadaFs_nP", "r1.dadaRs_nP", "r2.dadaFs_nP", "r2.dadaRs_nP",
              "r1.dadaFs_psP", "r1.dadaRs_psP", "r2.dadaFs_psP", "r2.dadaRs_psP", "r1.dadaFs_fP", "r1.dadaRs_fP", "r2.dadaFs_fP", "r2.dadaRs_fP"),
     file = file.path(RWORKSPACE_DIR, "Inference_dadaPooling.RData"),
     compress = "gzip")


######################## Read Tracking check up ################################

# Helper function to extract read counts from dada objects
getN <- function(x) sum(getUniques(x))

# NO POOLING (nP)
cat("Processing no pooling (nP) method...\n")

denoised_nP_r1 <- data.frame(
  sample = names(r1.dadaFs_nP),
  filtered = filtered_counts_r1$filtered[match(names(r1.dadaFs_nP), 
                                               gsub("\\.trim1\\.fq\\.gz$", "", filtered_counts_r1$sample))],
  denoisedF = sapply(r1.dadaFs_nP, getN),
  denoisedR = sapply(r1.dadaRs_nP, getN),
  stringsAsFactors = FALSE
)

denoised_nP_r2 <- data.frame(
  sample = names(r2.dadaFs_nP),
  filtered = filtered_counts_r2$filtered[match(names(r2.dadaFs_nP), 
                                               gsub("\\.trim1\\.fq\\.gz$", "", filtered_counts_r2$sample))],
  denoisedF = sapply(r2.dadaFs_nP, getN),
  denoisedR = sapply(r2.dadaRs_nP, getN),
  stringsAsFactors = FALSE
)

# PSEUDO-POOLING (psP)
cat("Processing pseudo-pooling (psP) method...\n")

denoised_psP_r1 <- data.frame(
  sample = names(r1.dadaFs_psP),
  filtered = filtered_counts_r1$filtered[match(names(r1.dadaFs_psP), 
                                               gsub("\\.trim1\\.fq\\.gz$", "", filtered_counts_r1$sample))],
  denoisedF = sapply(r1.dadaFs_psP, getN),
  denoisedR = sapply(r1.dadaRs_psP, getN),
  stringsAsFactors = FALSE
)

denoised_psP_r2 <- data.frame(
  sample = names(r2.dadaFs_psP),
  filtered = filtered_counts_r2$filtered[match(names(r2.dadaFs_psP), 
                                               gsub("\\.trim1\\.fq\\.gz$", "", filtered_counts_r2$sample))],
  denoisedF = sapply(r2.dadaFs_psP, getN),
  denoisedR = sapply(r2.dadaRs_psP, getN),
  stringsAsFactors = FALSE
)

# FULL POOLING (fP)
cat("Processing full pooling (fP) method...\n")

denoised_fP_r1 <- data.frame(
  sample = names(r1.dadaFs_fP),
  filtered = filtered_counts_r1$filtered[match(names(r1.dadaFs_fP), 
                                               gsub("\\.trim1\\.fq\\.gz$", "", filtered_counts_r1$sample))],
  denoisedF = sapply(r1.dadaFs_fP, getN),
  denoisedR = sapply(r1.dadaRs_fP, getN),
  stringsAsFactors = FALSE
)

denoised_fP_r2 <- data.frame(
  sample = names(r2.dadaFs_fP),
  filtered = filtered_counts_r2$filtered[match(names(r2.dadaFs_fP), 
                                               gsub("\\.trim1\\.fq\\.gz$", "", filtered_counts_r2$sample))],
  denoisedF = sapply(r2.dadaFs_fP, getN),
  denoisedR = sapply(r2.dadaRs_fP, getN),
  stringsAsFactors = FALSE
)

# Store in tracking object
read_tracking$step3_denoised <- list(
  nP = list(
    run1 = denoised_nP_r1,
    run2 = denoised_nP_r2,
    run1_total_denoisedF = sum(denoised_nP_r1$denoisedF),
    run1_total_denoisedR = sum(denoised_nP_r1$denoisedR),
    run2_total_denoisedF = sum(denoised_nP_r2$denoisedF),
    run2_total_denoisedR = sum(denoised_nP_r2$denoisedR)
  ),
  psP = list(
    run1 = denoised_psP_r1,
    run2 = denoised_psP_r2,
    run1_total_denoisedF = sum(denoised_psP_r1$denoisedF),
    run1_total_denoisedR = sum(denoised_psP_r1$denoisedR),
    run2_total_denoisedF = sum(denoised_psP_r2$denoisedF),
    run2_total_denoisedR = sum(denoised_psP_r2$denoisedR)
  ),
  fP = list(
    run1 = denoised_fP_r1,
    run2 = denoised_fP_r2,
    run1_total_denoisedF = sum(denoised_fP_r1$denoisedF),
    run1_total_denoisedR = sum(denoised_fP_r1$denoisedR),
    run2_total_denoisedF = sum(denoised_fP_r2$denoisedF),
    run2_total_denoisedR = sum(denoised_fP_r2$denoisedR)
  )
)

# Print summary for each method
cat("\n--- Denoising Summary ---\n")
for (method in c("nP", "psP", "fP")) {
  cat("\nMethod:", method, "\n")
  cat("  Run1 - Forward reads denoised:", 
      format(read_tracking$step3_denoised[[method]]$run1_total_denoisedF, big.mark = ","), "\n")
  cat("  Run1 - Reverse reads denoised:", 
      format(read_tracking$step3_denoised[[method]]$run1_total_denoisedR, big.mark = ","), "\n")
  cat("  Run2 - Forward reads denoised:", 
      format(read_tracking$step3_denoised[[method]]$run2_total_denoisedF, big.mark = ","), "\n")
  cat("  Run2 - Reverse reads denoised:", 
      format(read_tracking$step3_denoised[[method]]$run2_total_denoisedR, big.mark = ","), "\n")
}

# Save intermediate tracking
saveRDS(read_tracking, file = file.path(READ_TRACKING_DIR, "3.denoised_readTrack.rds"))




################################################################################
####################### Merging paired-ends ####################################

#load(file.path(RWORKSPACE_DIR, "Inference_dadaPooling.RData"))


merge_r1_nP = mergePairs(r1.dadaFs_nP, derepFwd_r1, r1.dadaRs_nP, derepRev_r1, verbose = TRUE)
merge_r1_psP = mergePairs(r1.dadaFs_psP, derepFwd_r1, r1.dadaRs_psP, derepRev_r1, verbose = TRUE)
merge_r1_fP = mergePairs(r1.dadaFs_fP, derepFwd_r1, r1.dadaRs_fP, derepRev_r1, verbose = TRUE)

merge_r2_nP = mergePairs(r2.dadaFs_nP, derepFwd_r2, r2.dadaRs_nP, derepRev_r2, verbose = TRUE)
merge_r2_psP = mergePairs(r2.dadaFs_psP, derepFwd_r2, r2.dadaRs_psP, derepRev_r2, verbose = TRUE)
merge_r2_fP = mergePairs(r2.dadaFs_fP, derepFwd_r2, r2.dadaRs_fP, derepRev_r2, verbose = TRUE)



# Combine all merge results into a single list
merged_results <- list(
  merge_r1_nP = merge_r1_nP,
  merge_r1_psP = merge_r1_psP,
  merge_r1_fP = merge_r1_fP,
  merge_r2_nP = merge_r2_nP,
  merge_r2_psP = merge_r2_psP,
  merge_r2_fP = merge_r2_fP
)

# Save the list as an RDS file
saveRDS(merged_results, file = file.path(RWORKSPACE_DIR, "Merged_results.rds"))


#Load merged objects
merged_results <- readRDS(file.path(RWORKSPACE_DIR, "Merged_results.rds"))

# Unpack list elements into your workspace
list2env(merged_results, envir = .GlobalEnv)

######################### Read Tracking Check up ###############################

# Helper function to get merged read counts
getN <- function(x) sum(getUniques(x))

# NO POOLING (nP)
cat("Processing no pooling (nP) method...\n")

merged_nP_r1 <- data.frame(
  sample = names(merge_r1_nP),
  denoisedF = denoised_nP_r1$denoisedF[match(names(merge_r1_nP), denoised_nP_r1$sample)],
  denoisedR = denoised_nP_r1$denoisedR[match(names(merge_r1_nP), denoised_nP_r1$sample)],
  merged = sapply(merge_r1_nP, getN),
  stringsAsFactors = FALSE
)
merged_nP_r1$percent_merged <- (merged_nP_r1$merged / pmin(merged_nP_r1$denoisedF, merged_nP_r1$denoisedR)) * 100

merged_nP_r2 <- data.frame(
  sample = names(merge_r2_nP),
  denoisedF = denoised_nP_r2$denoisedF[match(names(merge_r2_nP), denoised_nP_r2$sample)],
  denoisedR = denoised_nP_r2$denoisedR[match(names(merge_r2_nP), denoised_nP_r2$sample)],
  merged = sapply(merge_r2_nP, getN),
  stringsAsFactors = FALSE
)
merged_nP_r2$percent_merged <- (merged_nP_r2$merged / pmin(merged_nP_r2$denoisedF, merged_nP_r2$denoisedR)) * 100

# PSEUDO-POOLING (psP)
cat("Processing pseudo-pooling (psP) method...\n")

merged_psP_r1 <- data.frame(
  sample = names(merge_r1_psP),
  denoisedF = denoised_psP_r1$denoisedF[match(names(merge_r1_psP), denoised_psP_r1$sample)],
  denoisedR = denoised_psP_r1$denoisedR[match(names(merge_r1_psP), denoised_psP_r1$sample)],
  merged = sapply(merge_r1_psP, getN),
  stringsAsFactors = FALSE
)
merged_psP_r1$percent_merged <- (merged_psP_r1$merged / pmin(merged_psP_r1$denoisedF, merged_psP_r1$denoisedR)) * 100

merged_psP_r2 <- data.frame(
  sample = names(merge_r2_psP),
  denoisedF = denoised_psP_r2$denoisedF[match(names(merge_r2_psP), denoised_psP_r2$sample)],
  denoisedR = denoised_psP_r2$denoisedR[match(names(merge_r2_psP), denoised_psP_r2$sample)],
  merged = sapply(merge_r2_psP, getN),
  stringsAsFactors = FALSE
)
merged_psP_r2$percent_merged <- (merged_psP_r2$merged / pmin(merged_psP_r2$denoisedF, merged_psP_r2$denoisedR)) * 100

# FULL POOLING (fP)
cat("Processing full pooling (fP) method...\n")

merged_fP_r1 <- data.frame(
  sample = names(merge_r1_fP),
  denoisedF = denoised_fP_r1$denoisedF[match(names(merge_r1_fP), denoised_fP_r1$sample)],
  denoisedR = denoised_fP_r1$denoisedR[match(names(merge_r1_fP), denoised_fP_r1$sample)],
  merged = sapply(merge_r1_fP, getN),
  stringsAsFactors = FALSE
)
merged_fP_r1$percent_merged <- (merged_fP_r1$merged / pmin(merged_fP_r1$denoisedF, merged_fP_r1$denoisedR)) * 100

merged_fP_r2 <- data.frame(
  sample = names(merge_r2_fP),
  denoisedF = denoised_fP_r2$denoisedF[match(names(merge_r2_fP), denoised_fP_r2$sample)],
  denoisedR = denoised_fP_r2$denoisedR[match(names(merge_r2_fP), denoised_fP_r2$sample)],
  merged = sapply(merge_r2_fP, getN),
  stringsAsFactors = FALSE
)
merged_fP_r2$percent_merged <- (merged_fP_r2$merged / pmin(merged_fP_r2$denoisedF, merged_fP_r2$denoisedR)) * 100

# Store in tracking object
read_tracking$step4_merged <- list(
  nP = list(
    run1 = merged_nP_r1,
    run2 = merged_nP_r2,
    run1_total_merged = sum(merged_nP_r1$merged),
    run2_total_merged = sum(merged_nP_r2$merged),
    run1_merge_rate = mean(merged_nP_r1$percent_merged, na.rm = TRUE),
    run2_merge_rate = mean(merged_nP_r2$percent_merged, na.rm = TRUE)
  ),
  psP = list(
    run1 = merged_psP_r1,
    run2 = merged_psP_r2,
    run1_total_merged = sum(merged_psP_r1$merged),
    run2_total_merged = sum(merged_psP_r2$merged),
    run1_merge_rate = mean(merged_psP_r1$percent_merged, na.rm = TRUE),
    run2_merge_rate = mean(merged_psP_r2$percent_merged, na.rm = TRUE)
  ),
  fP = list(
    run1 = merged_fP_r1,
    run2 = merged_fP_r2,
    run1_total_merged = sum(merged_fP_r1$merged),
    run2_total_merged = sum(merged_fP_r2$merged),
    run1_merge_rate = mean(merged_fP_r1$percent_merged, na.rm = TRUE),
    run2_merge_rate = mean(merged_fP_r2$percent_merged, na.rm = TRUE)
  )
)

# Print summary
cat("\n--- Merging Summary ---\n")
for (method in c("nP", "psP", "fP")) {
  cat("\nMethod:", method, "\n")
  cat("  Run1:\n")
  cat("    Total merged reads:", 
      format(read_tracking$step4_merged[[method]]$run1_total_merged, big.mark = ","), "\n")
  cat("    Average merge rate:", 
      round(read_tracking$step4_merged[[method]]$run1_merge_rate, 2), "%\n")
  cat("  Run2:\n")
  cat("    Total merged reads:", 
      format(read_tracking$step4_merged[[method]]$run2_total_merged, big.mark = ","), "\n")
  cat("    Average merge rate:", 
      round(read_tracking$step4_merged[[method]]$run2_merge_rate, 2), "%\n")
}

# Save intermediate tracking
saveRDS(read_tracking, file = file.path(READ_TRACKING_DIR, "4.merged_readTrack.rds"))


################################################################################
####################### Make ASV tables ########################################

r1.asvTable_nP = makeSequenceTable(merge_r1_nP)
r1.asvTable_psP = makeSequenceTable(merge_r1_psP)
r1.asvTable_fP = makeSequenceTable(merge_r1_fP)

r2.asvTable_nP = makeSequenceTable(merge_r2_nP)
r2.asvTable_psP = makeSequenceTable(merge_r2_psP)
r2.asvTable_fP = makeSequenceTable(merge_r2_fP)

saveRDS(r1.asvTable_nP,  file = file.path(RWORKSPACE_DIR, "r1.asvTable_nP.rds"))
saveRDS(r1.asvTable_psP, file = file.path(RWORKSPACE_DIR, "r1.asvTable_psP.rds"))
saveRDS(r1.asvTable_fP,  file = file.path(RWORKSPACE_DIR, "r1.asvTable_fP.rds"))

saveRDS(r2.asvTable_nP,  file = file.path(RWORKSPACE_DIR, "r2.asvTable_nP.rds"))
saveRDS(r2.asvTable_psP, file = file.path(RWORKSPACE_DIR, "r2.asvTable_psP.rds"))
saveRDS(r2.asvTable_fP,  file = file.path(RWORKSPACE_DIR, "r2.asvTable_fP.rds"))

r1.asvTable_nP  = readRDS(file.path(RWORKSPACE_DIR, "r1.asvTable_nP.rds"))
r1.asvTable_psP = readRDS(file.path(RWORKSPACE_DIR, "r1.asvTable_psP.rds"))
r1.asvTable_fP  = readRDS(file.path(RWORKSPACE_DIR, "r1.asvTable_fP.rds"))

r2.asvTable_nP  = readRDS(file.path(RWORKSPACE_DIR, "r2.asvTable_nP.rds"))
r2.asvTable_psP = readRDS(file.path(RWORKSPACE_DIR, "r2.asvTable_psP.rds"))
r2.asvTable_fP  = readRDS(file.path(RWORKSPACE_DIR, "r2.asvTable_fP.rds"))

######## Differences in Pooling methods

#dimension = how many ASVs in each

dim(r1.asvTable_nP) # 2651 26219
dim(r1.asvTable_psP) #  2651 27376
dim(r1.asvTable_fP) #  2651 11882


dim(r2.asvTable_nP) #  2715 39722
dim(r2.asvTable_psP) # 2715 41647
dim(r2.asvTable_fP) # 2715 18745

## Full pooling has less ASVs than no pooling or pseudo-pooling
#No pooling: Calls many sample-specific rare variants (possibly including noise)
#Full pooling: More stringent, only calls ASVs supported across the entire dataset

#sequence length distribution

table(nchar(getSequences(r1.asvTable_nP))) 
table(nchar(getSequences(r1.asvTable_psP))) 
table(nchar(getSequences(r1.asvTable_fP))) 

table(nchar(getSequences(r2.asvTable_nP))) 
table(nchar(getSequences(r2.asvTable_psP))) 
table(nchar(getSequences(r2.asvTable_fP))) 



################################################################################
##################### Remove chimeras ##########################################


r1.asvTable_nochim_nP = removeBimeraDenovo(r1.asvTable_nP, method = "consensus", multithread = TRUE, verbose = TRUE)
r1.asvTable_nochim_psP = removeBimeraDenovo(r1.asvTable_psP, method = "consensus", multithread = TRUE, verbose = TRUE)
r1.asvTable_nochim_fP = removeBimeraDenovo(r1.asvTable_fP, method = "consensus", multithread = TRUE, verbose = TRUE)

r2.asvTable_nochim_nP = removeBimeraDenovo(r2.asvTable_nP, method = "consensus", multithread = TRUE, verbose = TRUE)
r2.asvTable_nochim_psP = removeBimeraDenovo(r2.asvTable_psP, method = "consensus", multithread = TRUE, verbose = TRUE)
r2.asvTable_nochim_fP = removeBimeraDenovo(r2.asvTable_fP, method = "consensus", multithread = TRUE, verbose = TRUE)

saveRDS(r1.asvTable_nochim_nP,  file = file.path(RWORKSPACE_DIR, "r1.asvTable_nochim_nP.rds"))
saveRDS(r1.asvTable_nochim_psP, file = file.path(RWORKSPACE_DIR, "r1.asvTable_nochim_psP.rds"))
saveRDS(r1.asvTable_nochim_fP,  file = file.path(RWORKSPACE_DIR, "r1.asvTable_nochim_fP.rds"))

saveRDS(r2.asvTable_nochim_nP,  file = file.path(RWORKSPACE_DIR, "r2.asvTable_nochim_nP.rds"))
saveRDS(r2.asvTable_nochim_psP, file = file.path(RWORKSPACE_DIR, "r2.asvTable_nochim_psP.rds"))
saveRDS(r2.asvTable_nochim_fP,  file = file.path(RWORKSPACE_DIR, "r2.asvTable_nochim_fP.rds"))


r1.asvTable_nochim_nP  = readRDS(file.path(RWORKSPACE_DIR, "r1.asvTable_nochim_nP.rds"))
r1.asvTable_nochim_psP = readRDS(file.path(RWORKSPACE_DIR, "r1.asvTable_nochim_psP.rds"))
r1.asvTable_nochim_fP  = readRDS(file.path(RWORKSPACE_DIR, "r1.asvTable_nochim_fP.rds"))

r2.asvTable_nochim_nP  = readRDS(file.path(RWORKSPACE_DIR, "r2.asvTable_nochim_nP.rds"))
r2.asvTable_nochim_psP = readRDS(file.path(RWORKSPACE_DIR, "r2.asvTable_nochim_psP.rds"))
r2.asvTable_nochim_fP  = readRDS(file.path(RWORKSPACE_DIR, "r2.asvTable_nochim_fP.rds"))




######### Proportion of reads retained after chimera removal

sum(r1.asvTable_nochim_nP)/sum(r1.asvTable_nP) #  0.9979975
sum(r1.asvTable_nochim_psP)/sum(r1.asvTable_psP) # 0.9975914
sum(r1.asvTable_nochim_fP)/sum(r1.asvTable_fP) # 0.9978662

sum(r2.asvTable_nochim_nP)/sum(r2.asvTable_nP) #  0.9962997
sum(r2.asvTable_nochim_psP)/sum(r2.asvTable_psP) #  0.9952712
sum(r2.asvTable_nochim_fP)/sum(r2.asvTable_fP) # 0.9957366

#### ASV counts

dim(r1.asvTable_nochim_nP) 
dim(r1.asvTable_nochim_psP) 
dim(r1.asvTable_nochim_fP) 

dim(r2.asvTable_nochim_nP) 
dim(r2.asvTable_nochim_psP) 
dim(r2.asvTable_nochim_fP)

#################### Tracking read loss ########################################

# NO POOLING (nP)
cat("Processing no pooling (nP) method...\n")

# Get per-sample counts from ASV tables
chimera_nP_r1 <- data.frame(
  sample = rownames(r1.asvTable_nochim_nP),
  merged = merged_nP_r1$merged[match(rownames(r1.asvTable_nochim_nP), merged_nP_r1$sample)],
  nonchimeric = rowSums(r1.asvTable_nochim_nP),
  stringsAsFactors = FALSE
)
chimera_nP_r1$chimeric_removed <- chimera_nP_r1$merged - chimera_nP_r1$nonchimeric
chimera_nP_r1$percent_retained <- (chimera_nP_r1$nonchimeric / chimera_nP_r1$merged) * 100

chimera_nP_r2 <- data.frame(
  sample = rownames(r2.asvTable_nochim_nP),
  merged = merged_nP_r2$merged[match(rownames(r2.asvTable_nochim_nP), merged_nP_r2$sample)],
  nonchimeric = rowSums(r2.asvTable_nochim_nP),
  stringsAsFactors = FALSE
)
chimera_nP_r2$chimeric_removed <- chimera_nP_r2$merged - chimera_nP_r2$nonchimeric
chimera_nP_r2$percent_retained <- (chimera_nP_r2$nonchimeric / chimera_nP_r2$merged) * 100

# PSEUDO-POOLING (psP)
cat("Processing pseudo-pooling (psP) method...\n")

chimera_psP_r1 <- data.frame(
  sample = rownames(r1.asvTable_nochim_psP),
  merged = merged_psP_r1$merged[match(rownames(r1.asvTable_nochim_psP), merged_psP_r1$sample)],
  nonchimeric = rowSums(r1.asvTable_nochim_psP),
  stringsAsFactors = FALSE
)
chimera_psP_r1$chimeric_removed <- chimera_psP_r1$merged - chimera_psP_r1$nonchimeric
chimera_psP_r1$percent_retained <- (chimera_psP_r1$nonchimeric / chimera_psP_r1$merged) * 100

chimera_psP_r2 <- data.frame(
  sample = rownames(r2.asvTable_nochim_psP),
  merged = merged_psP_r2$merged[match(rownames(r2.asvTable_nochim_psP), merged_psP_r2$sample)],
  nonchimeric = rowSums(r2.asvTable_nochim_psP),
  stringsAsFactors = FALSE
)
chimera_psP_r2$chimeric_removed <- chimera_psP_r2$merged - chimera_psP_r2$nonchimeric
chimera_psP_r2$percent_retained <- (chimera_psP_r2$nonchimeric / chimera_psP_r2$merged) * 100

# FULL POOLING (fP)
cat("Processing full pooling (fP) method...\n")

chimera_fP_r1 <- data.frame(
  sample = rownames(r1.asvTable_nochim_fP),
  merged = merged_fP_r1$merged[match(rownames(r1.asvTable_nochim_fP), merged_fP_r1$sample)],
  nonchimeric = rowSums(r1.asvTable_nochim_fP),
  stringsAsFactors = FALSE
)
chimera_fP_r1$chimeric_removed <- chimera_fP_r1$merged - chimera_fP_r1$nonchimeric
chimera_fP_r1$percent_retained <- (chimera_fP_r1$nonchimeric / chimera_fP_r1$merged) * 100

chimera_fP_r2 <- data.frame(
  sample = rownames(r2.asvTable_nochim_fP),
  merged = merged_fP_r2$merged[match(rownames(r2.asvTable_nochim_fP), merged_fP_r2$sample)],
  nonchimeric = rowSums(r2.asvTable_nochim_fP),
  stringsAsFactors = FALSE
)
chimera_fP_r2$chimeric_removed <- chimera_fP_r2$merged - chimera_fP_r2$nonchimeric
chimera_fP_r2$percent_retained <- (chimera_fP_r2$nonchimeric / chimera_fP_r2$merged) * 100

# Store in tracking object
read_tracking$step5_nochimera <- list(
  nP = list(
    run1 = chimera_nP_r1,
    run2 = chimera_nP_r2,
    run1_total_nonchimeric = sum(chimera_nP_r1$nonchimeric),
    run2_total_nonchimeric = sum(chimera_nP_r2$nonchimeric),
    run1_total_chimeric_removed = sum(chimera_nP_r1$chimeric_removed),
    run2_total_chimeric_removed = sum(chimera_nP_r2$chimeric_removed),
    run1_retention_rate = (sum(chimera_nP_r1$nonchimeric) / sum(chimera_nP_r1$merged)) * 100,
    run2_retention_rate = (sum(chimera_nP_r2$nonchimeric) / sum(chimera_nP_r2$merged)) * 100,
    run1_n_asvs = ncol(r1.asvTable_nochim_nP),
    run2_n_asvs = ncol(r2.asvTable_nochim_nP)
  ),
  psP = list(
    run1 = chimera_psP_r1,
    run2 = chimera_psP_r2,
    run1_total_nonchimeric = sum(chimera_psP_r1$nonchimeric),
    run2_total_nonchimeric = sum(chimera_psP_r2$nonchimeric),
    run1_total_chimeric_removed = sum(chimera_psP_r1$chimeric_removed),
    run2_total_chimeric_removed = sum(chimera_psP_r2$chimeric_removed),
    run1_retention_rate = (sum(chimera_psP_r1$nonchimeric) / sum(chimera_psP_r1$merged)) * 100,
    run2_retention_rate = (sum(chimera_psP_r2$nonchimeric) / sum(chimera_psP_r2$merged)) * 100,
    run1_n_asvs = ncol(r1.asvTable_nochim_psP),
    run2_n_asvs = ncol(r2.asvTable_nochim_psP)
  ),
  fP = list(
    run1 = chimera_fP_r1,
    run2 = chimera_fP_r2,
    run1_total_nonchimeric = sum(chimera_fP_r1$nonchimeric),
    run2_total_nonchimeric = sum(chimera_fP_r2$nonchimeric),
    run1_total_chimeric_removed = sum(chimera_fP_r1$chimeric_removed),
    run2_total_chimeric_removed = sum(chimera_fP_r2$chimeric_removed),
    run1_retention_rate = (sum(chimera_fP_r1$nonchimeric) / sum(chimera_fP_r1$merged)) * 100,
    run2_retention_rate = (sum(chimera_fP_r2$nonchimeric) / sum(chimera_fP_r2$merged)) * 100,
    run1_n_asvs = ncol(r1.asvTable_nochim_fP),
    run2_n_asvs = ncol(r2.asvTable_nochim_fP)
  )
)

# Print summary
cat("\n--- Chimera Removal Summary ---\n")
for (method in c("nP", "psP", "fP")) {
  cat("\nMethod:", method, "\n")
  cat("  Run1:\n")
  cat("    Non-chimeric reads:", 
      format(read_tracking$step5_nochimera[[method]]$run1_total_nonchimeric, big.mark = ","), "\n")
  cat("    Chimeric reads removed:", 
      format(read_tracking$step5_nochimera[[method]]$run1_total_chimeric_removed, big.mark = ","), "\n")
  cat("    Retention rate:", 
      round(read_tracking$step5_nochimera[[method]]$run1_retention_rate, 2), "%\n")
  cat("    Number of ASVs:", 
      format(read_tracking$step5_nochimera[[method]]$run1_n_asvs, big.mark = ","), "\n")
  cat("  Run2:\n")
  cat("    Non-chimeric reads:", 
      format(read_tracking$step5_nochimera[[method]]$run2_total_nonchimeric, big.mark = ","), "\n")
  cat("    Chimeric reads removed:", 
      format(read_tracking$step5_nochimera[[method]]$run2_total_chimeric_removed, big.mark = ","), "\n")
  cat("    Retention rate:", 
      round(read_tracking$step5_nochimera[[method]]$run2_retention_rate, 2), "%\n")
  cat("    Number of ASVs:", 
      format(read_tracking$step5_nochimera[[method]]$run2_n_asvs, big.mark = ","), "\n")
}

# Save intermediate tracking
saveRDS(read_tracking, file = file.path(READ_TRACKING_DIR, "5.nochimera_readTrack.rds"))

################################################################################


rm("derepFwd_r1", "derepRev_r1", "derepFwd_r2", "derepRev_r2", "r1.dadaFs_nP", "r1.dadaRs_nP", "r2.dadaFs_psP", "r2.dadaRs_psP", "r1.dadaFs_fP", "r1.dadaRs_fP", 
   "r2.dadaFs_nP", "r2.dadaRs_nP", "r1.dadaFs_psP", "r1.dadaRs_psP", "r2.dadaFs_fP", "r2.dadaRs_fP")
gc()


################################################################################
################# Merging sequencing Runs ######################################

# Function to extract biological sample ID (removes sequencing run info)
extract_bio_sample <- function(full_name) {
  # Extract everything before _MKD (sample name) and the replicate info
  # EB_bc12_MKDL250008975-1A_22WW55LT4_L6__.rep4 -> EB_bc12_rep4
  
  # First, extract sample name before _MKD
  sample_part <- sub("_MKD.*", "", full_name)
  
  # Extract replicate number (rep1-8)
  rep_part <- regmatches(full_name, regexpr("rep[0-9]+", full_name))
  
  # Combine them
  bio_id <- paste(sample_part, rep_part, sep = "_")
  
  return(bio_id)
}


# NO POOLING (nP)
cat("Merging nP method...\n")
# Rename rows to biological IDs
rownames(r1.asvTable_nochim_nP) <- sapply(rownames(r1.asvTable_nochim_nP), extract_bio_sample)
rownames(r2.asvTable_nochim_nP) <- sapply(rownames(r2.asvTable_nochim_nP), extract_bio_sample)

# Merge using DADA2 function (repeats="sum" will sum counts for duplicate sample names)
nP.seqMerged <- mergeSequenceTables(r1.asvTable_nochim_nP, r2.asvTable_nochim_nP, repeats = "sum", orderBy = NULL)
cat("  nP merged:", nrow(nP.seqMerged), "samples,", ncol(nP.seqMerged), "sequences\n")

# nP merged: 2769 samples, 42987 sequences

# PSEUDO-POOLING (psP)
cat("Merging psP method...\n")
rownames(r1.asvTable_nochim_psP) <- sapply(rownames(r1.asvTable_nochim_psP), extract_bio_sample)
rownames(r2.asvTable_nochim_psP) <- sapply(rownames(r2.asvTable_nochim_psP), extract_bio_sample)

psP.seqMerged <- mergeSequenceTables(r1.asvTable_nochim_psP, r2.asvTable_nochim_psP, repeats = "sum", orderBy = NULL)
cat("  psP merged:", nrow(psP.seqMerged), "samples,", ncol(psP.seqMerged), "sequences\n")

# psP merged: 2769 samples, 44709 sequences

# FULL POOLING (fP)
cat("Merging fP method...\n")
rownames(r1.asvTable_nochim_fP) <- sapply(rownames(r1.asvTable_nochim_fP), extract_bio_sample)
rownames(r2.asvTable_nochim_fP) <- sapply(rownames(r2.asvTable_nochim_fP), extract_bio_sample)

fP.seqMerged <- mergeSequenceTables(r1.asvTable_nochim_fP, r2.asvTable_nochim_fP, repeats = "sum", orderBy = NULL)
cat("  fP merged:", nrow(fP.seqMerged), "samples,", ncol(fP.seqMerged), "sequences\n")

#  fP merged: 2769 samples, 19676 sequences

saveRDS(nP.seqMerged,  file = file.path(RWORKSPACE_DIR, "nP_seqMerged.rds"))
saveRDS(psP.seqMerged, file = file.path(RWORKSPACE_DIR, "psP_seqMerged.rds"))
saveRDS(fP.seqMerged,  file = file.path(RWORKSPACE_DIR, "fP_seqMerged.rds"))

nP.seqMerged  = readRDS(file.path(RWORKSPACE_DIR, "nP_seqMerged.rds"))
psP.seqMerged = readRDS(file.path(RWORKSPACE_DIR, "psP_seqMerged.rds"))
fP.seqMerged  = readRDS(file.path(RWORKSPACE_DIR, "fP_seqMerged.rds"))


save.image("MAPLE_18SEUK_dada2.RData")

############################# Read tracking check up ###########################


# NO POOLING (nP)
merged_nP_stats <- list(
  # Before merging (separate runs)
  r1_samples = nrow(r1.asvTable_nochim_nP),
  r2_samples = nrow(r2.asvTable_nochim_nP),
  r1_reads = sum(r1.asvTable_nochim_nP),
  r2_reads = sum(r2.asvTable_nochim_nP),
  r1_sequences = ncol(r1.asvTable_nochim_nP),
  r2_sequences = ncol(r2.asvTable_nochim_nP),
  
  # After merging
  merged_samples = nrow(nP.seqMerged),
  merged_reads = sum(nP.seqMerged),
  merged_sequences = ncol(nP.seqMerged),
  
  # Verify no reads lost
  total_reads_before = sum(r1.asvTable_nochim_nP) + sum(r2.asvTable_nochim_nP),
  reads_retained_pct = (sum(nP.seqMerged) / (sum(r1.asvTable_nochim_nP) + sum(r2.asvTable_nochim_nP))) * 100
)

# PSEUDO-POOLING (psP)
merged_psP_stats <- list(
  # Before merging
  r1_samples = nrow(r1.asvTable_nochim_psP),
  r2_samples = nrow(r2.asvTable_nochim_psP),
  r1_reads = sum(r1.asvTable_nochim_psP),
  r2_reads = sum(r2.asvTable_nochim_psP),
  r1_sequences = ncol(r1.asvTable_nochim_psP),
  r2_sequences = ncol(r2.asvTable_nochim_psP),
  
  # After merging
  merged_samples = nrow(psP.seqMerged),
  merged_reads = sum(psP.seqMerged),
  merged_sequences = ncol(psP.seqMerged),
  
  # Verify
  total_reads_before = sum(r1.asvTable_nochim_psP) + sum(r2.asvTable_nochim_psP),
  reads_retained_pct = (sum(psP.seqMerged) / (sum(r1.asvTable_nochim_psP) + sum(r2.asvTable_nochim_psP))) * 100
)

# FULL POOLING (fP)
merged_fP_stats <- list(
  # Before merging
  r1_samples = nrow(r1.asvTable_nochim_fP),
  r2_samples = nrow(r2.asvTable_nochim_fP),
  r1_reads = sum(r1.asvTable_nochim_fP),
  r2_reads = sum(r2.asvTable_nochim_fP),
  r1_sequences = ncol(r1.asvTable_nochim_fP),
  r2_sequences = ncol(r2.asvTable_nochim_fP),
  
  # After merging
  merged_samples = nrow(fP.seqMerged),
  merged_reads = sum(fP.seqMerged),
  merged_sequences = ncol(fP.seqMerged),
  
  # Verify
  total_reads_before = sum(r1.asvTable_nochim_fP) + sum(r2.asvTable_nochim_fP),
  reads_retained_pct = (sum(fP.seqMerged) / (sum(r1.asvTable_nochim_fP) + sum(r2.asvTable_nochim_fP))) * 100
)

# Store in tracking object
read_tracking$step7_merged_runs <- list(
  nP = merged_nP_stats,
  psP = merged_psP_stats,
  fP = merged_fP_stats
)

# Print summary
cat("\n--- Merging Runs Summary ---\n")
for (method in c("nP", "psP", "fP")) {
  cat("\nMethod:", method, "\n")
  cat("  Before merging:\n")
  cat("    Run1:", 
      read_tracking$step7_merged_runs[[method]]$r1_samples, "samples,",
      format(read_tracking$step7_merged_runs[[method]]$r1_reads, big.mark = ","), "reads,",
      format(read_tracking$step7_merged_runs[[method]]$r1_sequences, big.mark = ","), "sequences\n")
  cat("    Run2:", 
      read_tracking$step7_merged_runs[[method]]$r2_samples, "samples,",
      format(read_tracking$step7_merged_runs[[method]]$r2_reads, big.mark = ","), "reads,",
      format(read_tracking$step7_merged_runs[[method]]$r2_sequences, big.mark = ","), "sequences\n")
  cat("    Total reads:", 
      format(read_tracking$step7_merged_runs[[method]]$total_reads_before, big.mark = ","), "\n")
  
  cat("  After merging:\n")
  cat("    Merged:", 
      read_tracking$step7_merged_runs[[method]]$merged_samples, "unique biological samples,",
      format(read_tracking$step7_merged_runs[[method]]$merged_reads, big.mark = ","), "reads,",
      format(read_tracking$step7_merged_runs[[method]]$merged_sequences, big.mark = ","), "sequences\n")
  cat("    Reads retained:", 
      sprintf("%.2f%%", read_tracking$step7_merged_runs[[method]]$reads_retained_pct), "\n")
}

# Save intermediate tracking
saveRDS(read_tracking, file = file.path(READ_TRACKING_DIR, "6.merged_runs_readTrack.rds"))



#################### Prepping for LULU #########################################

########## Prepping files

## Export sequences as FASTA files

uniquesToFasta(nP.seqMerged,  file.path(LULU_DIR, "nP_dada2.fasta"),
               ids = paste0("OTU",seq(length(getSequences(nP.seqMerged)))))
uniquesToFasta(psP.seqMerged, file.path(LULU_DIR, "psP_dada2.fasta"),
               ids = paste0("OTU",seq(length(getSequences(psP.seqMerged)))))
uniquesToFasta(fP.seqMerged,  file.path(LULU_DIR, "fP_dada2.fasta"),
               ids = paste0("OTU",seq(length(getSequences(fP.seqMerged)))))


## dada asv -to-> Lulu OTU tables (row = OTUs, col = samples)

nP.lulu_ready = nP.seqMerged
colnames(nP.lulu_ready) = paste0("OTU", seq(length(getSequences(nP.seqMerged))))
nP.lulu_ready = t(nP.lulu_ready)

psP.lulu_ready = psP.seqMerged
colnames(psP.lulu_ready) = paste0("OTU", seq(length(getSequences(psP.seqMerged))))
psP.lulu_ready = t(psP.lulu_ready)

fP.lulu_ready = fP.seqMerged
colnames(fP.lulu_ready) = paste0("OTU", seq(length(getSequences(fP.seqMerged))))
fP.lulu_ready = t(fP.lulu_ready)


saveRDS(nP.lulu_ready,  file = file.path(RWORKSPACE_DIR, "nP_lulu_ready.rds"))
saveRDS(psP.lulu_ready, file = file.path(RWORKSPACE_DIR, "psP_lulu_ready.rds"))
saveRDS(fP.lulu_ready,  file = file.path(RWORKSPACE_DIR, "fP_lulu_ready.rds"))

#################################### Astbury HPC Terminal
##BASH
'''
cd data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/01.RawData/LULU_fasta/
  
  # Making BLAST database with OTUs
  
makeblastdb -in fP_dada2.fasta -parse_seqids -dbtype nucl
makeblastdb -in nP_dada2.fasta -parse_seqids -dbtype nucl
makeblastdb -in psP_dada2.fasta -parse_seqids -dbtype nucl

# Blast OTUs against new database (independent datasets)

blastn -db fP_dada2.fasta -outfmt '6 qseqid sseqid pident' -out fP_dada2_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query fP_dada2.fasta -num_threads 16
blastn -db nP_dada2.fasta -outfmt '6 qseqid sseqid pident' -out nP_dada2_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query nP_dada2.fasta -num_threads 16
blastn -db psP_dada2.fasta  -outfmt '6 qseqid sseqid pident' -out psP_dada2_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query psP_dada2.fasta -num_threads 16

'''
################################################################################
######################## LULU algorithm in R ###################################
################################################################################


nP.lulu_ready  = readRDS(file.path(RWORKSPACE_DIR, "nP_lulu_ready.rds"))
psP.lulu_ready = readRDS(file.path(RWORKSPACE_DIR, "psP_lulu_ready.rds"))
fP.lulu_ready  = readRDS(file.path(RWORKSPACE_DIR, "fP_lulu_ready.rds"))

load("MAPLE_18SEUK_dada2.RData")

setwd(LULU_DIR)

library(devtools)
#devtools::install_github("benjjneb/dada2", ref = "master")
library(dada2)
#install_github("tobiasgf/lulu")
library(lulu)
#install.packages("dplyr")
library(dplyr)

#reading txt files from dir
nP_dada2_match_list  = read.table(file.path(LULU_DIR, "nP_dada2_match_list.txt"))
psP_dada2_match_list = read.table(file.path(LULU_DIR, "psP_dada2_match_list.txt"))
fP_dada2_match_list  = read.table(file.path(LULU_DIR, "fP_dada2_match_list.txt"))

#Applying lulu algorithm
nP.nochim.curated_result  = lulu(as.data.frame(nP.lulu_ready),  nP_dada2_match_list)
psP.nochim.curated_result = lulu(as.data.frame(psP.lulu_ready), psP_dada2_match_list)
fP.nochim.curated_result  = lulu(as.data.frame(fP.lulu_ready),  fP_dada2_match_list)


saveRDS(nP.nochim.curated_result,  file = file.path(LULU_DIR, "nP_LULU_curated.rds"))
saveRDS(psP.nochim.curated_result, file = file.path(LULU_DIR, "psP_LULU_curated.rds"))
saveRDS(fP.nochim.curated_result,  file = file.path(LULU_DIR, "fP_LULU_curated.rds"))


nP.nochim.curated_result  = readRDS(file.path(LULU_DIR, "nP_LULU_curated.rds"))
psP.nochim.curated_result = readRDS(file.path(LULU_DIR, "psP_LULU_curated.rds"))
fP.nochim.curated_result  = readRDS(file.path(LULU_DIR, "fP_LULU_curated.rds"))


#Checking LULU's effect
print(paste0("Not Pooled: ", "OTUs after Lulu: ", nP.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(nP.nochim.curated_result$original_table)))
#"Not Pooled: OTUs after Lulu: 27460 --- OTUs before Lulu: 42987
print(paste0("Pooled: ", "OTUs after Lulu: ", fP.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(fP.nochim.curated_result$original_table)))
#Pooled: OTUs after Lulu: 16337 --- OTUs before Lulu: 19676
print(paste0("PseudoPooled: ", "OTUs after Lulu: ", psP.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(psP.nochim.curated_result$original_table)))
#PseudoPooled: OTUs after Lulu: 28719 --- OTUs before Lulu: 44709



######################## read loss tracking check up ###########################

# NO POOLING (nP)
cat("Processing no pooling (nP) method...\n")

# LULU works on the combined runs, so we track at the combined level
# Calculate total reads before and after LULU
lulu_nP_stats <- list(
  # Before LULU (from merged runs)
  asvs_before = ncol(nP.seqMerged),
  reads_before = sum(nP.seqMerged),
  
  
  # After LULU
  otus_after = nP.nochim.curated_result$curated_count,
  reads_after = sum(nP.nochim.curated_result$curated_table),
  asvs_discarded = nrow(nP.nochim.curated_result$original_table) - nP.nochim.curated_result$curated_count,
  reads_retained_pct = (sum(nP.nochim.curated_result$curated_table) / sum(nP.nochim.curated_result$original_table)) * 100
)

# PSEUDO-POOLING (psP)
cat("Processing pseudo-pooling (psP) method...\n")

lulu_psP_stats <- list(
  # Before LULU
  asvs_before = ncol(psP.seqMerged),
  reads_before = sum(psP.seqMerged),
  
  # After LULU
  otus_after = psP.nochim.curated_result$curated_count,
  reads_after = sum(psP.nochim.curated_result$curated_table),
  asvs_discarded = nrow(psP.nochim.curated_result$original_table) - psP.nochim.curated_result$curated_count,
  reads_retained_pct = (sum(psP.nochim.curated_result$curated_table) / sum(psP.nochim.curated_result$original_table)) * 100
)

# FULL POOLING (fP)
cat("Processing full pooling (fP) method...\n")

lulu_fP_stats <- list(
  # Before LULU
  asvs_before = ncol(fP.seqMerged),
  reads_before = sum(fP.seqMerged),
  
  # After LULU
  otus_after = fP.nochim.curated_result$curated_count,
  reads_after = sum(fP.nochim.curated_result$curated_table),
  asvs_discarded = nrow(fP.nochim.curated_result$original_table) - fP.nochim.curated_result$curated_count,
  reads_retained_pct = (sum(fP.nochim.curated_result$curated_table) / sum(fP.nochim.curated_result$original_table)) * 100
)

# Store in tracking object
read_tracking$step6_lulu <- list(
  nP = lulu_nP_stats,
  psP = lulu_psP_stats,
  fP = lulu_fP_stats
)

# Print summary
cat("\n--- LULU Curation Summary ---\n")
for (method in c("nP", "psP", "fP")) {
  cat("\nMethod:", method, "\n")
  cat("  Run1:\n")
  cat("    ASVs before LULU:", 
      format(read_tracking$step6_lulu[[method]]$asvs_before, big.mark = ","), "\n")
  cat("    OTUs after LULU:", 
      format(read_tracking$step6_lulu[[method]]$otus_after, big.mark = ","), "\n")
  cat("    ASVs discarded:", 
      format(read_tracking$step6_lulu[[method]]$asvs_discarded, big.mark = ","), 
      sprintf(" (%.1f%%)", (read_tracking$step6_lulu[[method]]$asvs_discarded / 
                              read_tracking$step6_lulu[[method]]$asvs_before) * 100), "\n")
  cat("    Reads retained:", 
      format(read_tracking$step6_lulu[[method]]$reads_after, big.mark = ","),
      sprintf(" (%.2f%%)", read_tracking$step6_lulu[[method]]$reads_retained_pct), "\n")
  
}

# Save intermediate tracking
saveRDS(read_tracking, file = file.path(READ_TRACKING_DIR, "7.lulu_readTrack.rds"))

################################################################################

##### Examine structure of curated tables
rownames(nP.nochim.curated_result$curated_table) ## OTUs names after lulu are not sequencial !! They keep they're original name

#### Extract OTU names to numeric indices -> this returns vectors with only numeric indices of OTUs
nP.keptOTU  = as.numeric(gsub("OTU", "", rownames(nP.nochim.curated_result$curated_table),  ignore.case = FALSE, perl = TRUE))
psP.keptOTU = as.numeric(gsub("OTU", "", rownames(psP.nochim.curated_result$curated_table), ignore.case = FALSE, perl = TRUE))
fP.keptOTU  = as.numeric(gsub("OTU", "", rownames(fP.nochim.curated_result$curated_table),  ignore.case = FALSE, perl = TRUE))


## Return to DADA2 OTU tables format: OTU = columns
nP.lulu = t(nP.nochim.curated_result$curated_table)
colnames(nP.lulu) = colnames(nP.seqMerged[, nP.keptOTU])

psP.lulu = t(psP.nochim.curated_result$curated_table)
colnames(psP.lulu) = colnames(psP.seqMerged[, psP.keptOTU])

fP.lulu = t(fP.nochim.curated_result$curated_table)
colnames(fP.lulu) = colnames(fP.seqMerged[, fP.keptOTU])


saveRDS(nP.lulu,  file = file.path(RWORKSPACE_DIR, "nP_lulu.rds"))
saveRDS(psP.lulu, file = file.path(RWORKSPACE_DIR, "psP_lulu.rds"))
saveRDS(fP.lulu,  file = file.path(RWORKSPACE_DIR, "fP_lulu.rds"))


save.image(file.path(RAW_DATA_DIR, "MAPLE_18SEUK_LULUdone.RData"))


################################################################################
## Export to FASTA

# Function to write sequences to FASTA
write_fasta <- function(sequences, filename) {
  # Get unique sequences (column names are the sequences)
  unique_seqs <- colnames(sequences)
  
  # Create sequence IDs
  seq_ids <- paste0("ASV_", 1:length(unique_seqs))
  
  # Write FASTA file
  fasta_lines <- c()
  for(i in 1:length(unique_seqs)) {
    fasta_lines <- c(fasta_lines, paste0(">", seq_ids[i]), unique_seqs[i])
  }
  
  writeLines(fasta_lines, filename)
  cat("Written", length(unique_seqs), "sequences to", filename, "\n")
}

# Convert your tables to FASTA
write_fasta(nP.lulu,  file.path(LULU_DIR, "nP_notSwarmed.fasta"))
write_fasta(psP.lulu, file.path(LULU_DIR, "psP_notSwarmed.fasta"))
write_fasta(fP.lulu,  file.path(LULU_DIR, "fP_notSwarmed.fasta"))

###############################################################################
###################### SWARM clustering ########################################

setwd(RAW_DATA_DIR)

# IMPROVED SWARM WORKFLOW FOR METABARCODING DATA
# ================================================

# 1. IMPROVED FASTA FORMATTER
# ---------------------------
swarm_fasta_formatter <- function(seqtab, filename) {
  abundances <- colSums(seqtab)
  sequences <- colnames(seqtab)
  
  # Create fasta with abundance info in swarm format
  fasta_lines <- c()
  for(i in 1:length(sequences)) {
    # Use swarm format: >ASV_ID_abundance
    header <- paste0(">ASV_", i, "_", abundances[i])
    fasta_lines <- c(fasta_lines, header, sequences[i])
  }
  
  writeLines(fasta_lines, filename)
  cat("Written", length(sequences), "sequences to", filename, "\n")
}

# Apply function to files
#dir.create("swarm", showWarnings = FALSE)
swarm_fasta_formatter(nP.lulu,  file.path(SWARM_DIR, "nP_swarm.fasta"))
swarm_fasta_formatter(psP.lulu, file.path(SWARM_DIR, "psP_swarm.fasta"))
swarm_fasta_formatter(fP.lulu,  file.path(SWARM_DIR, "fP_swarm.fasta"))

# 2. FUNCTION TO CREATE SEQUENCE-TO-INDEX MAPPING
# -----------------------------------------------
create_sequence_mapping <- function(seqtab) {
  sequences <- colnames(seqtab)
  mapping <- data.frame(
    asv_id = paste0("ASV_", 1:length(sequences)),
    sequence = sequences,
    abundance = colSums(seqtab),
    stringsAsFactors = FALSE
  )
  return(mapping)
}

# Create mappings for later use
nP_mapping  <- create_sequence_mapping(nP.lulu)
psP_mapping <- create_sequence_mapping(psP.lulu)
fP_mapping  <- create_sequence_mapping(fP.lulu)

# 4. IMPROVED FUNCTION TO CONVERT SWARM CLUSTERS TO OTU TABLE
# -----------------------------------------------------------
swarm_to_otu <- function(cluster_file, original_seqtab, rep_seq_file = NULL, otu_prefix = "OTU") {
  
  # Check if files exist
  if(!file.exists(cluster_file)) {
    stop("Cluster file not found: ", cluster_file)
  }
  
  # Read cluster file
  clusters <- readLines(cluster_file)
  clusters <- clusters[clusters != ""]  # Remove empty lines
  
  if(length(clusters) == 0) {
    stop("No clusters found in file: ", cluster_file)
  }
  
  # Initialize OTU table
  otu_table <- matrix(0, nrow = nrow(original_seqtab), ncol = length(clusters))
  rownames(otu_table) <- rownames(original_seqtab)
  colnames(otu_table) <- paste0(otu_prefix, "_", sprintf("%04d", 1:length(clusters)))
  
  # Store representative sequences for taxonomy assignment
  rep_sequences <- character(length(clusters))
  names(rep_sequences) <- colnames(otu_table)
  
  # Process each cluster
  for(i in 1:length(clusters)) {
    # Split cluster into individual ASV names
    asv_names <- trimws(strsplit(clusters[i], "[ \t]+")[[1]])
    asv_names <- asv_names[asv_names != ""]
    
    if(length(asv_names) == 0) next
    
    # Extract ASV indices and find sequences
    asv_indices <- c()
    cluster_sequences <- c()
    
    for(asv_name in asv_names) {
      # Extract ASV number from "ASV_X_abundance" format
      parts <- strsplit(asv_name, "_")[[1]]
      if(length(parts) >= 2) {
        asv_idx <- as.numeric(parts[2])
        if(!is.na(asv_idx) && asv_idx <= ncol(original_seqtab)) {
          asv_indices <- c(asv_indices, asv_idx)
          cluster_sequences <- c(cluster_sequences, colnames(original_seqtab)[asv_idx])
        }
      }
    }
    
    if(length(asv_indices) > 0) {
      # Sum abundances across all ASVs in this cluster
      otu_table[, i] <- rowSums(original_seqtab[, asv_indices, drop = FALSE])
      
      # Use the most abundant sequence as representative (first in swarm output)
      rep_sequences[i] <- cluster_sequences[1]
    }
  }
  
  # Remove empty OTUs
  non_zero_otus <- colSums(otu_table) > 0
  otu_table <- otu_table[, non_zero_otus, drop = FALSE]
  rep_sequences <- rep_sequences[non_zero_otus]
  
  # Return list with OTU table and representative sequences
  result <- list(
    otu_table = otu_table,
    representative_sequences = rep_sequences,
    n_otus = ncol(otu_table),
    n_original_asvs = ncol(original_seqtab)
  )
  
  return(result)
}

# 4. FUNCTION TO WRITE REPRESENTATIVE SEQUENCES FOR TAXONOMY
# ----------------------------------------------------------
write_rep_sequences <- function(rep_sequences, filename) {
  fasta_lines <- c()
  for(i in 1:length(rep_sequences)) {
    header <- paste0(">", names(rep_sequences)[i])
    fasta_lines <- c(fasta_lines, header, rep_sequences[i])
  }
  writeLines(fasta_lines, filename)
  cat("Written", length(rep_sequences), "representative sequences to", filename, "\n")
}

# 5. PROCESS YOUR DATASETS (run after swarm clustering)
# bash script: 02.SwarmClustering_EUK.sh in repo


# AFTER CLUSTERING IN BASH
# -----------------------------------------------------
# Convert your datasets and extract representative sequences
nP_result  <- swarm_to_otu(file.path(SWARM_DIR, "nP_post_swarm_clusters.txt"),  nP.lulu,  otu_prefix = "nP_OTU")
psP_result <- swarm_to_otu(file.path(SWARM_DIR, "psP_post_swarm_clusters.txt"), psP.lulu, otu_prefix = "psP_OTU")
fP_result  <- swarm_to_otu(file.path(SWARM_DIR, "fP_post_swarm_clusters.txt"),  fP.lulu,  otu_prefix = "fP_OTU")

# Extract OTU tables
nP_otu_table  <- nP_result$otu_table
psP_otu_table <- psP_result$otu_table
fP_otu_table  <- fP_result$otu_table

# Write representative sequences for taxonomic assignment
write_rep_sequences(nP_result$representative_sequences,  file.path(SWARM_DIR, "nP_rep_sequences.fasta"))
write_rep_sequences(psP_result$representative_sequences, file.path(SWARM_DIR, "psP_rep_sequences.fasta"))
write_rep_sequences(fP_result$representative_sequences,  file.path(SWARM_DIR, "fP_rep_sequences.fasta"))

# 6. SUMMARY STATISTICS
# ---------------------
cat("SWARM CLUSTERING RESULTS:\n")
cat("=========================\n")
cat("nP dataset: ", ncol(nP.lulu), "ASVs ->", nP_result$n_otus, "OTUs\n")
cat("psP dataset:", ncol(psP.lulu), "ASVs ->", psP_result$n_otus, "OTUs\n")
cat("fP dataset: ", ncol(fP.lulu), "ASVs ->", fP_result$n_otus, "OTUs\n")

# nP dataset:  22470 ASVs -> 14662 OTUs
# psP dataset: 23157 ASVs -> 15163 OTUs
# fP dataset:  16337 ASVs -> 14930 OTUs


# 7. SAVE RESULTS
# -------------------------
# Save OTU tables
write.csv(nP_otu_table,  file.path(SWARM_DIR, "nP_otu_table.csv"))
write.csv(psP_otu_table, file.path(SWARM_DIR, "psP_otu_table.csv"))
write.csv(fP_otu_table,  file.path(SWARM_DIR, "fP_otu_table.csv"))

save(nP_result, psP_result, fP_result, 
     nP_otu_table, psP_otu_table, fP_otu_table,
     file = file.path(RWORKSPACE_DIR, "swarm_results.RData"))

save.image(file.path(RAW_DATA_DIR, "MAPLE_18SEUK_SWARMdone.RData"))

############################## Read loss tracking ##############################

# NO POOLING (nP)
cat("Processing no pooling (nP) method...\n")

swarm_nP_stats <- list(
  # Before SWARM (ASV level after LULU)
  asvs_before = ncol(nP.lulu),
  reads_before = sum(nP.lulu),
  samples_before = nrow(nP.lulu),
  
  # After SWARM (OTU level)
  otus_after = nP_result$n_otus,
  reads_after = sum(nP_result$otu_table),
  samples_after = nrow(nP_result$otu_table),
  
  # Calculate reduction
  asv_to_otu_reduction = ncol(nP.lulu) - nP_result$n_otus,
  asv_to_otu_reduction_pct = ((ncol(nP.lulu) - nP_result$n_otus) / ncol(nP.lulu)) * 100,
  reads_retained_pct = (sum(nP_result$otu_table) / sum(nP.lulu)) * 100
)

# PSEUDO-POOLING (psP)
cat("Processing pseudo-pooling (psP) method...\n")

swarm_psP_stats <- list(
  # Before SWARM
  asvs_before = ncol(psP.lulu),
  reads_before = sum(psP.lulu),
  samples_before = nrow(psP.lulu),
  
  # After SWARM
  otus_after = psP_result$n_otus,
  reads_after = sum(psP_result$otu_table),
  samples_after = nrow(psP_result$otu_table),
  
  # Calculate reduction
  asv_to_otu_reduction = ncol(psP.lulu) - psP_result$n_otus,
  asv_to_otu_reduction_pct = ((ncol(psP.lulu) - psP_result$n_otus) / ncol(psP.lulu)) * 100,
  reads_retained_pct = (sum(psP_result$otu_table) / sum(psP.lulu)) * 100
)

# FULL POOLING (fP)
cat("Processing full pooling (fP) method...\n")

swarm_fP_stats <- list(
  # Before SWARM
  asvs_before = ncol(fP.lulu),
  reads_before = sum(fP.lulu),
  samples_before = nrow(fP.lulu),
  
  # After SWARM
  otus_after = fP_result$n_otus,
  reads_after = sum(fP_result$otu_table),
  samples_after = nrow(fP_result$otu_table),
  
  # Calculate reduction
  asv_to_otu_reduction = ncol(fP.lulu) - fP_result$n_otus,
  asv_to_otu_reduction_pct = ((ncol(fP.lulu) - fP_result$n_otus) / ncol(fP.lulu)) * 100,
  reads_retained_pct = (sum(fP_result$otu_table) / sum(fP.lulu)) * 100
)

# Store in tracking object
read_tracking$step8_swarm <- list(
  nP = swarm_nP_stats,
  psP = swarm_psP_stats,
  fP = swarm_fP_stats
)

# Print summary
cat("\n--- SWARM Clustering Summary ---\n")
for (method in c("nP", "psP", "fP")) {
  cat("\nMethod:", method, "\n")
  cat("  Before SWARM (ASV level):\n")
  cat("    Samples:", read_tracking$step8_swarm[[method]]$samples_before, "\n")
  cat("    ASVs:", format(read_tracking$step8_swarm[[method]]$asvs_before, big.mark = ","), "\n")
  cat("    Reads:", format(read_tracking$step8_swarm[[method]]$reads_before, big.mark = ","), "\n")
  
  cat("  After SWARM (OTU level):\n")
  cat("    Samples:", read_tracking$step8_swarm[[method]]$samples_after, "\n")
  cat("    OTUs:", format(read_tracking$step8_swarm[[method]]$otus_after, big.mark = ","), "\n")
  cat("    Reads:", format(read_tracking$step8_swarm[[method]]$reads_after, big.mark = ","), "\n")
  
  cat("  Clustering effect:\n")
  cat("    ASVs clustered into OTUs:", 
      format(read_tracking$step8_swarm[[method]]$asv_to_otu_reduction, big.mark = ","),
      sprintf(" (%.1f%% reduction)", read_tracking$step8_swarm[[method]]$asv_to_otu_reduction_pct), "\n")
  cat("    Reads retained:", 
      sprintf("%.2f%%", read_tracking$step8_swarm[[method]]$reads_retained_pct), "\n")
}

# Save final tracking object
saveRDS(read_tracking, file = file.path(READ_TRACKING_DIR, "8.swarm_readTrack_FINAL.rds"))



################################################################################
##################### Assign Taxonomy ##########################################



###################### DADA2 taxonomy assignment
#using Silva Eukaryotic 18S (v132)

#WITHOUT SWARM

taxa132_nP  = assignTaxonomy(nP.lulu,  SILVA_DB_PATH, multithread = 16, verbose = TRUE)
saveRDS(taxa132_nP,  file = file.path(TAXA_DB_DIR, "taxa132_nP.rds"))

taxa132_psP = assignTaxonomy(psP.lulu, SILVA_DB_PATH, multithread = 16, verbose = TRUE)
saveRDS(taxa132_psP, file = file.path(TAXA_DB_DIR, "taxa132_psP.rds"))

taxa132_fP  = assignTaxonomy(fP.lulu,  SILVA_DB_PATH, multithread = 16, verbose = TRUE)
saveRDS(taxa132_fP,  file = file.path(TAXA_DB_DIR, "taxa132_fP.rds"))


taxa132_nP  = readRDS(file.path(TAXA_DB_DIR, "taxa132_nP.rds"))
taxa132_psP = readRDS(file.path(TAXA_DB_DIR, "taxa132_psP.rds"))
taxa132_fP  = readRDS(file.path(TAXA_DB_DIR, "taxa132_fP.rds"))

#WITH SWARM

load(file.path(RWORKSPACE_DIR, "swarm_results.RData"))

taxa132_sw_nP  = assignTaxonomy(nP_result$representative_sequences,  SILVA_DB_PATH, multithread = 16, verbose = TRUE)
saveRDS(taxa132_sw_nP,  file = file.path(TAXA_DB_DIR, "taxa132_swarmed_nP.rds"))

taxa132_sw_psP = assignTaxonomy(psP_result$representative_sequences, SILVA_DB_PATH, multithread = 16, verbose = TRUE)
saveRDS(taxa132_sw_psP, file = file.path(TAXA_DB_DIR, "taxa132_swarmed_psP.rds"))

taxa132_sw_fP  = assignTaxonomy(fP_result$representative_sequences,  SILVA_DB_PATH, multithread = 16, verbose = TRUE)
saveRDS(taxa132_sw_fP,  file = file.path(TAXA_DB_DIR, "taxa132_swarmed_fP.rds"))


taxa132_sw_nP  = readRDS(file.path(TAXA_DB_DIR, "taxa132_swarmed_nP.rds"))
taxa132_sw_psP = readRDS(file.path(TAXA_DB_DIR, "taxa132_swarmed_psP.rds"))
taxa132_sw_fP  = readRDS(file.path(TAXA_DB_DIR, "taxa132_swarmed_fP.rds"))


########################## NCBI BLAST taxonomy assignment

#Export representative sequences as FASTA files
# For each dataset, write FASTA files
library(Biostrings)

# Export representative sequences
writeXStringSet(DNAStringSet(nP_result$representative_sequences),  file.path(NCBI_DIR, "nP_repseq.fasta"))
writeXStringSet(DNAStringSet(psP_result$representative_sequences), file.path(NCBI_DIR, "psP_repseq.fasta"))
writeXStringSet(DNAStringSet(fP_result$representative_sequences),  file.path(NCBI_DIR, "fP_repseq.fasta"))



