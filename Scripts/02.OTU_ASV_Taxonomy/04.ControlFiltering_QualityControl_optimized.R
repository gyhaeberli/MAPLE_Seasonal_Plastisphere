################################################################################
############## OTU TABLE FILTERING PIPELINE - v2
############## optimized by ClaudeAI (Anthropic) Sonnet 4.6
################################################################################
#
# PIPELINE OVERVIEW:
#   1.  Load raw OTU tables (LULU + SWARM x nP / psP / fP)
#   2.  Demultiplex EUK1 library sample names
#   3.  Build metadata (index.info)
#   4.  Filter positive controls          → remove OTUs present in 2+ reps of any positive
#   5.  Subtract NTC contamination        → per-OTU max NTC read count subtracted
#   6.  Correct extraction blanks         → per cycling group, per sample type
#   7.  Keep environmental samples only   → remove controls from OTU tables
#   8.  QUALITY FILTER                    → remove PCR replicates with < min_OTUs
#                                            OR reads < min_pct_sibling_median
#   9.  Rep filters (50%, 75%)            → applied to quality-filtered clean tables
#  10.  Collapse PCR replicates           → sum across reps per bio replicate
#  11.  Remove zero-abundance OTUs        → from all versions
#  12.  Save all output tables
#  13.  Filtering diagnostics + tracking
#
# Changes from v1:
#   - All 6 methods processed through shared function loops (no repeated blocks)
#   - Quality filter (step 8) added BEFORE rep filters so rep fractions are
#     calculated on clean replicates only
#   - Collapsed tables derived from quality-filtered clean tables (same source
#     as rep-filtered tables), ensuring consistency across all versions
#   - Zero-OTU removal consolidated into one loop
#   - Diagnostic summary uses a single generalized function
#
################################################################################

################################################################################
# PATH CONFIGURATION
# All project paths are defined here. Edit this section to relocate the project.
################################################################################

# ---- Project root ----
proj_root  <- "~/GitHub/MAPLE_Seasonal_Plastisphere"

# ---- Working directory ----
path_scripts_dir        <- file.path(proj_root, "Scripts/02.OTU_ASV_Taxonomy")

# ---- Input: raw OTU tables ----
path_nP_lulu_rds        <- file.path(proj_root, "Raw_data/OTU_tables/nP_lulu.rds")
path_psP_lulu_rds       <- file.path(proj_root, "Raw_data/OTU_tables/psP_lulu.rds")
path_fP_lulu_rds        <- file.path(proj_root, "Raw_data/OTU_tables/fP_lulu.rds")
path_swarm_rdata        <- file.path(proj_root, "Raw_data/OTU_tables/swarm_results.RData")

# ---- Input: metadata Excel ----
path_euk1_metadata      <- "EUK1_metadata.xlsx"   # relative to working directory

# ---- Output: demultiplexed OTU tables ----
dir_demult              <- file.path(proj_root, "Processed_data/OTU_tables")
path_demult_nP_lulu     <- file.path(dir_demult, "nP.lulu_original_demult.rds")
path_demult_psP_lulu    <- file.path(dir_demult, "psP.lulu_original_demult.rds")
path_demult_fP_lulu     <- file.path(dir_demult, "fP.lulu_original_demult.rds")
path_demult_nP_sw       <- file.path(dir_demult, "nP.sw_original_demult.rds")
path_demult_psP_sw      <- file.path(dir_demult, "psP.sw_original_demult.rds")
path_demult_fP_sw       <- file.path(dir_demult, "fP.sw_original_demult.rds")

# ---- Output: metadata RDS ----
dir_metadata            <- file.path(proj_root, "Raw_data/Metadata")
path_meta_nP_lulu       <- file.path(dir_metadata, "nP.lulu_metadata.rds")
path_meta_psP_lulu      <- file.path(dir_metadata, "psP.lulu_metadata.rds")
path_meta_fP_lulu       <- file.path(dir_metadata, "fP.lulu_metadata.rds")
path_meta_nP_sw         <- file.path(dir_metadata, "nP.swarm_metadata.rds")
path_meta_psP_sw        <- file.path(dir_metadata, "psP.swarm_metadata.rds")
path_meta_fP_sw         <- file.path(dir_metadata, "fP.swarm_metadata.rds")

# ---- Output: filtered OTU tables ----
dir_filtered            <- file.path(proj_root, "Processed_data/OTU_tables/filtered")
path_qf_nP_lulu         <- file.path(dir_filtered, "LULU_nP_qualityfilt.rds")
path_qf_psP_lulu        <- file.path(dir_filtered, "LULU_psP_qualityfilt.rds")
path_qf_fP_lulu         <- file.path(dir_filtered, "LULU_fP_qualityfilt.rds")
path_qf_nP_sw           <- file.path(dir_filtered, "SWARM_nP_qualityfilt.rds")
path_qf_psP_sw          <- file.path(dir_filtered, "SWARM_psP_qualityfilt.rds")
path_qf_fP_sw           <- file.path(dir_filtered, "SWARM_fP_qualityfilt.rds")

# ---- Output: sanity check / diagnostics ----
dir_sanity              <- file.path(proj_root, "Results/SanityChecks")
path_pos_ctrl_log       <- file.path(dir_sanity, "positive_control_filter_log.csv")
path_qf_removed_log     <- file.path(dir_sanity, "quality_filter_removed_pcr_reps.csv")
path_tag_jump_csv       <- file.path(dir_sanity, "tag_jumping_rates.csv")
path_rep_sim_summary    <- file.path(dir_sanity, "replicate_similarity_summary.csv")

# ---- Output: decontam library size plots ----
dir_decontam_plots      <- file.path(dir_sanity, "decontam_plots")

# ---- Output: read tracking ----
dir_read_tracking       <- file.path(dir_sanity, "Read_tracking")
path_filter_summary_csv <- file.path(dir_read_tracking, "filter_summary_ALL_METHODS.csv")
path_filter_tracking    <- file.path(dir_read_tracking, "filter_tracking_object.rds")


# ---- Output: workspace image ----
path_rdata              <- "OTU_filtering_v2.RData"   # relative to working directory

# ---- Create all output directories ----
dirs_to_create <- c(
  dir_demult,
  dir_metadata,
  dir_filtered,
  dir_sanity,
  dir_decontam_plots,
  dir_read_tracking,
  dir_rep_sim
)
for (d in dirs_to_create) dir.create(d, showWarnings = FALSE, recursive = TRUE)

################################################################################

setwd(path_scripts_dir)

library(dplyr)
library(stringr)
library(readxl)
library(ggplot2)
library(knitr)
library(kableExtra)
library(tidyr)

out_dir <- dir_filtered

################################################################################
# 1.  LOAD RAW OTU TABLES
################################################################################

nP.lulu  <- readRDS(path_nP_lulu_rds)
psP.lulu <- readRDS(path_psP_lulu_rds)
fP.lulu  <- readRDS(path_fP_lulu_rds)

load(path_swarm_rdata)
nP.sw  <- nP_result$otu_table;  colnames(nP.sw)  <- nP_result$representative_sequences[colnames(nP.sw)]
psP.sw <- psP_result$otu_table; colnames(psP.sw) <- psP_result$representative_sequences[colnames(psP.sw)]
fP.sw  <- fP_result$otu_table;  colnames(fP.sw)  <- fP_result$representative_sequences[colnames(fP.sw)]

################################################################################
# 2.  DEMULTIPLEX EUK1 LIBRARY
################################################################################

EUK1_sample_map        <- read_excel(path_euk1_metadata, sheet = "Sheet1")
colnames(EUK1_sample_map) <- c("wells", "sample_names", "rep_sets")

demultiplex_names <- function(matrix_obj, sample_map) {
  rnames   <- rownames(matrix_obj)
  gab_rows <- grep("GAB_EUK", rnames)
  
  for (i in gab_rows) {
    old_name  <- rnames[i]
    well      <- sub(".*_([A-H][0-9]{1,2})_.*", "\\1", old_name)
    rep_num   <- as.numeric(sub(".*_rep([0-9]+)$", "\\1", old_name))
    rep_set   <- ifelse(rep_num %in% 1:4, "set1", "set2")
    match_idx <- which(sample_map$wells == well & sample_map$rep_sets == rep_set)
    
    if (length(match_idx) > 0) {
      rnames[i] <- sub("GAB_EUK1_[A-H][0-9]{1,2}",
                       sample_map$sample_names[match_idx[1]], old_name)
    } else {
      warning(paste("No match found for well", well, "rep", rep_num, rep_set, "in row", i))
    }
  }
  rownames(matrix_obj) <- rnames
  return(matrix_obj)
}

# Ensure unique row names before demultiplexing
make_unique_rownames <- function(matrix_obj) {
  rnames   <- rownames(matrix_obj)
  dup_names <- unique(rnames[duplicated(rnames)])
  for (dup in dup_names) {
    positions <- which(rnames == dup)
    for (i in seq_along(positions)) rnames[positions[i]] <- paste0(dup, "_", i)
  }
  rownames(matrix_obj) <- rnames
  return(matrix_obj)
}

raw_list <- list(
  nP.lulu  = nP.lulu,  psP.lulu = psP.lulu, fP.lulu = fP.lulu,
  nP.sw    = nP.sw,    psP.sw   = psP.sw,   fP.sw   = fP.sw
)

demult_list <- lapply(raw_list, function(mat) {
  mat <- demultiplex_names(mat, EUK1_sample_map)
  make_unique_rownames(mat)
})

# Check that no GAB_EUK rows remain after demultiplexing
cat("=== DEMULTIPLEX CHECK ===\n")
for (nm in names(demult_list)) {
  remaining <- grep("GAB_EUK", rownames(demult_list[[nm]]))
  cat(nm, "- GAB_EUK rows remaining:", length(remaining), "\n")
  if (length(remaining) > 0) print(rownames(demult_list[[nm]])[remaining])
}

demult_paths <- list(
  nP.lulu  = path_demult_nP_lulu,
  psP.lulu = path_demult_psP_lulu,
  fP.lulu  = path_demult_fP_lulu,
  nP.sw    = path_demult_nP_sw,
  psP.sw   = path_demult_psP_sw,
  fP.sw    = path_demult_fP_sw
)
for (nm in names(demult_list)) saveRDS(demult_list[[nm]], demult_paths[[nm]])
cat("Demultiplexed OTU tables saved.\n")


################################################################################
# 3.  BUILD METADATA
################################################################################

plate_map     <- read_excel(path_euk1_metadata, sheet = "Sheet4")
cyc_group_map <- read_excel(path_euk1_metadata, sheet = "Sheet3")
date_metadata <- readxl::read_excel(path_euk1_metadata, sheet = "Sheet5")

library_map <- c(P1 = "GAB_EUK1", P2 = "GAB_EUK1",
                 P3 = "GAB_EUK2", P4 = "GAB_EUK2",
                 P5 = "GAB_EUK3", P6 = "GAB_EUK3",
                 P7 = "GAB_EUK4")

safe_date_convert <- function(x) {
  if (inherits(x, "Date"))                    return(x)
  if (inherits(x, "POSIXct") | inherits(x, "POSIXt")) return(as.Date(x))
  numeric_dates <- suppressWarnings(as.numeric(x))
  if (!all(is.na(numeric_dates)) &&
      all(numeric_dates[!is.na(numeric_dates)] > 1 &
          numeric_dates[!is.na(numeric_dates)] < 100000))
    return(as.Date(numeric_dates, origin = "1899-12-30"))
  return(as.Date(rep(NA_character_, length(x))))
}

date_metadata <- date_metadata %>%
  mutate(
    deployment_date = safe_date_convert(deployment_date),
    retrieval_date  = safe_date_convert(retrieval_date),
    site = case_when(
      site == "SEL" ~ "SELVA",
      TRUE          ~ site
    )
  )

index.info <- function(x, plate_map = NULL, library_map = NULL,
                       cyc_group_map = NULL, date_metadata = NULL) {
  
  y            <- data.frame(matrix(NA, nrow = nrow(x), ncol = 0))
  sample_names <- rownames(x)
  
  y$clean_sample_names <- sub("_rep[0-9]+$", "", sample_names)
  
  y$site <- case_when(
    grepl("^TBS_", sample_names) ~ "TBS",
    grepl("^SEL_", sample_names) ~ "SELVA",
    grepl("^EB_|^MB_", sample_names) ~ "LAB",
    TRUE ~ NA_character_)
  
  y$sample_type <- case_when(
    grepl("B(0[1-9]|1[0-2])($|[^0-9])", sample_names)          ~ "Blank",
    grepl("NTC", sample_names, ignore.case = TRUE)               ~ "NTC",
    grepl("_[0-9]+_Fil", sample_names)                          ~ "FilterBlank",
    grepl("blk|blank|Xblk", sample_names, ignore.case = TRUE)   ~ "ExtractBlank",
    grepl("FISH", sample_names, ignore.case = TRUE)              ~ "PositiveControl",
    grepl("^EB_|^MB_", sample_names)                            ~ "DegradationSample",
    grepl("^SEL_|^TBS_", sample_names)                          ~ "EnvironmentalSample",
    TRUE ~ "Unknown")
  
  y$substrate <- case_when(
    grepl("_G[0-9]",    sample_names) ~ "Glass",
    grepl("_wPE[0-9]",  sample_names) ~ "Weathered_PE",
    grepl("_PE[0-9]",   sample_names) ~ "PE",
    grepl("_wPET[0-9]", sample_names) ~ "Weathered_PET",
    grepl("_PET[0-9]",  sample_names) ~ "PET",
    grepl("_F[0-9]",    sample_names) ~ "Filter",
    grepl("_bc*",       sample_names) ~ "sterilized_water",
    grepl("_iso*",      sample_names) ~ "13C-PE",
    TRUE ~ NA_character_)
  
  y$bio_replicate <- as.numeric(stringr::str_extract(
    sample_names, "(?<=(_G|_PE|_wPE|_PET|_wPET|_F))[0-9]+"))
  special_idx <- grepl("^(EB|MB)_(bc|iso)", sample_names, ignore.case = TRUE)
  y$bio_replicate[special_idx] <- as.numeric(
    stringr::str_extract(sample_names[special_idx], "(?<=_(bc|iso))[0-9]"))
  y$bio_replicate <- dplyr::case_when(
    grepl("^(EB|MB)_", sample_names, ignore.case = TRUE) &
      y$bio_replicate %in% 1:6 ~ y$bio_replicate,
    y$bio_replicate %in% 1:3   ~ y$bio_replicate,
    TRUE ~ NA_real_)
  
  y$PCR_rep <- stringr::str_extract(sample_names, "rep[0-9]+")
  y$adapt_primer <- case_when(
    y$PCR_rep == "rep1" ~ "A1", y$PCR_rep == "rep2" ~ "B2",
    y$PCR_rep == "rep3" ~ "C3", y$PCR_rep == "rep4" ~ "D4",
    y$PCR_rep == "rep5" ~ "E5", y$PCR_rep == "rep6" ~ "F6",
    y$PCR_rep == "rep7" ~ "G7", y$PCR_rep == "rep8" ~ "H8",
    TRUE ~ NA_character_)
  
  y$sampling_event <- dplyr::case_when(
    grepl("^TBS_[0-9]+_|^SEL_[0-9]+_", sample_names) ~
      as.numeric(stringr::str_extract(sample_names, "(?<=TBS_|SEL_)[0-9]+")),
    grepl("_WIN_|_WIN$",  sample_names) ~ 3,
    grepl("_SPR_|_SPR$",  sample_names) ~ 6,
    grepl("_SUM_|_SUM$",  sample_names) ~ 8,
    grepl("_FAL_|_FAL$",  sample_names) ~ 10,
    grepl("_WIN2_|_WIN2$",sample_names) ~ 12,
    TRUE ~ NA_real_)
  
  y$season <- case_when(
    grepl("_SPR_|_SPR$",  sample_names) ~ "Spring",
    grepl("_SUM_|_SUM$",  sample_names) ~ "Summer",
    grepl("_FAL_|_FAL$",  sample_names) ~ "Fall",
    grepl("_WIN2_|_WIN2$",sample_names) ~ "Winter2",
    grepl("_WIN_|_WIN$",  sample_names) ~ "Winter",
    y$sampling_event == 3 & grepl("^TBS_3_|^SEL_3_", sample_names) ~ "Winter",
    TRUE ~ NA_character_)
  
  y$experiment <- case_when(
    y$sample_type %in% c("Blank","NTC","FilterBlank","ExtractBlank","PositiveControl") ~ NA_character_,
    grepl("^EB_|^MB_", sample_names)                                   ~ "Degradation",
    grepl("_F[0-9]", sample_names) & y$sample_type == "EnvironmentalSample" ~ "Water",
    !is.na(y$season)                                                   ~ "ST",
    grepl("^SEL_|^TBS_", sample_names) & is.na(y$season) &
      y$sample_type != "FilterBlank"                                   ~ "LT",
    TRUE ~ NA_character_)
  
  y$st_role <- dplyr::case_when(
    y$substrate != "Filter"                          ~ NA_character_,
    y$sampling_event %in% c(1, 5, 7, 9, 11)         ~ "Deployment",
    y$sampling_event %in% c(3, 6, 8, 10, 12)        ~ "Retrieval",
    TRUE                                             ~ NA_character_
  )
  
  y$season <- dplyr::case_when(
    !is.na(y$season)                        ~ y$season,
    y$sampling_event %in% c(1, 3)           ~ "Winter",
    y$sampling_event %in% c(5, 6)           ~ "Spring",
    y$sampling_event %in% c(7, 8)           ~ "Summer",
    y$sampling_event %in% c(9, 10)          ~ "Fall",
    y$sampling_event %in% c(11, 12)         ~ "Winter2",
    TRUE                                     ~ NA_character_
  )
  
  if (!is.null(plate_map)) {
    y$plate <- sapply(sample_names, function(sample) {
      core_name <- sub("_rep[0-9]+(_[0-9]+)?$", "", sample)
      for (plate_col in colnames(plate_map)) {
        if (core_name %in% plate_map[[plate_col]])
          return(gsub("plate", "P", plate_col))
      }
      return(NA_character_)
    })
    y$well_position <- sapply(sample_names, function(sample) {
      core_name <- sub("_rep[0-9]+(_[0-9]+)?$", "", sample)
      for (plate_col in colnames(plate_map)) {
        if (plate_col == "well_position") next
        if (core_name %in% plate_map[[plate_col]]) {
          row_idx <- which(plate_map[[plate_col]] == core_name)
          return(plate_map$well_position[row_idx])
        }
      }
      return(NA_character_)
    })
  } else {
    y$plate <- NA_character_
    y$well_position <- NA_character_
  }
  
  y$library <- if (!is.null(library_map)) library_map[y$plate] else
    stringr::str_extract(sample_names, "GAB_EUK[0-9]+")
  
  if (!is.null(cyc_group_map)) {
    y$cycling_group <- sapply(sample_names, function(sample) {
      core_name <- sub("_rep[0-9]+(_[0-9]+)?$", "", sample)
      for (group_col in colnames(cyc_group_map)) {
        if (core_name %in% cyc_group_map[[group_col]])
          return(sub("group_", "", group_col))
      }
      return(NA_character_)
    })
  } else {
    y$cycling_group <- NA_character_
  }
  
  # ---- DATE MATCHING ----
  if (!is.null(date_metadata)) {
    
    retrieval_event_lookup <- c(
      "1" = 3,  "3" = 3,
      "5" = 6,  "6" = 6,
      "7" = 8,  "8" = 8,
      "9" = 10, "10" = 10,
      "11" = 12,"12" = 12
    )
    
    deployment_event_lookup <- c(
      "1" = 1,  "3" = 1,
      "5" = 5,  "6" = 5,
      "7" = 7,  "8" = 7,
      "9" = 9,  "10" = 9,
      "11" = 11,"12" = 11
    )
    
    base_sample <- gsub("[_-](G|PE|wPE|PET|wPET|F)[0-9]+$", "", y$clean_sample_names)
    y$deployment_date <- rep(as.Date(NA_character_), length(base_sample))
    y$retrieval_date  <- rep(as.Date(NA_character_), length(base_sample))
    
    for (i in seq_along(base_sample)) {
      
      # Determine which experiment label to use for matching
      # LT samples match LT rows, ST biofilm matches ST rows, water matches Water rows
      exp_label <- dplyr::case_when(
        y$experiment[i] == "LT"    ~ "LT",
        y$experiment[i] == "ST"    ~ "ST",
        y$experiment[i] == "Water" ~ "Water",
        TRUE                       ~ NA_character_
      )
      
      if (is.na(exp_label) || is.na(y$site[i]) || is.na(y$sampling_event[i])) next
      
      event_key <- as.character(y$sampling_event[i])
      
      # --- Get deployment date ---
      dep_event <- deployment_event_lookup[event_key]
      if (!is.na(dep_event)) {
        dep_row <- which(date_metadata$sampling_event == dep_event &
                           date_metadata$site         == y$site[i] &
                           date_metadata$experiment   == exp_label)
        if (length(dep_row) > 0) {
          dep_dates <- date_metadata$deployment_date[dep_row]
          y$deployment_date[i] <- dep_dates[!is.na(dep_dates)][1]
        }
      }
      
      # --- Get retrieval date ---
      ret_event <- retrieval_event_lookup[event_key]
      if (!is.na(ret_event)) {
        ret_row <- which(date_metadata$sampling_event == ret_event &
                           date_metadata$site         == y$site[i] &
                           date_metadata$experiment   == exp_label)
        if (length(ret_row) > 0) {
          ret_dates <- date_metadata$retrieval_date[ret_row]
          y$retrieval_date[i] <- ret_dates[!is.na(ret_dates)][1]
        }
      }
    }
  }
  
  
  y$time_in_water_days <- as.numeric(difftime(y$retrieval_date, y$deployment_date, units = "days"))
  y$full   <- sample_names
  y$totseq <- rowSums(x)
  y$OTUs   <- rowSums(x > 0)
  rownames(y) <- sample_names
  return(y)
}

index_list <- lapply(demult_list, function(mat) {
  index.info(mat, plate_map = plate_map, library_map = library_map,
             cyc_group_map = cyc_group_map, date_metadata = date_metadata)
})

# Save metadata
meta_paths <- list(
  nP.lulu  = path_meta_nP_lulu,
  psP.lulu = path_meta_psP_lulu,
  fP.lulu  = path_meta_fP_lulu,
  nP.sw    = path_meta_nP_sw,
  psP.sw   = path_meta_psP_sw,
  fP.sw    = path_meta_fP_sw
)
for (nm in names(index_list)) saveRDS(index_list[[nm]], meta_paths[[nm]])






################################################################################
# SHARED GROUPING HELPER
# Groups are defined as cycling_group x plate combinations.
# This is the correct unit for all control corrections since cycling groups
# 12, 15, and 26 span multiple plates, each with their own controls.
################################################################################

get_groups <- function(index_data, rows) {
  cg    <- index_data[rows, "cycling_group"]
  plate <- index_data[rows, "plate"]
  ifelse(!is.na(cg) & !is.na(plate), paste(cg, plate, sep = "_"), NA_character_)
}


################################################################################
# 4.  FILTER POSITIVE CONTROLS
# Per cycling_group x plate group:
#   Step 1 — reproducibility filter: flag OTUs present in 2+ PCR replicates
#            of any positive control sample within the group
#   Step 2 — directionality filter: of those, only remove OTUs where mean
#            reads in the positive control are >= min_pos_ctrl_ratio times
#            higher than mean reads in environmental samples of the same group.
#            This protects abundant environmental taxa that leaked into the
#            positive control via tag jumping — if an OTU is far more abundant
#            in the environment than in the positive control, the contamination
#            went that direction, not the other way.
# Positive control rows are removed entirely after correction.
################################################################################

filter_positive <- function(mat, index_data, min_pos_ctrl_ratio = 10,
                            method_name = "unnamed") {
  rows         <- rownames(mat)
  sample_types <- index_data[rows, "sample_type"]
  groups       <- get_groups(index_data, rows)
  is_positive  <- sample_types == "PositiveControl"
  
  if (sum(is_positive) == 0) {
    warning("No positive control samples found — skipping filter")
    return(list(mat = mat[!is_positive, , drop = FALSE], log = NULL))
  }
  
  unique_groups <- unique(groups[is_positive & !is.na(groups)])
  mat_corrected <- mat
  is_env        <- sample_types %in% c("EnvironmentalSample", "DegradationSample")
  
  removal_log <- data.frame(
    method           = character(),
    cycling_group    = character(),
    otu              = character(),
    max_rep_count    = integer(),
    pos_mean_reads   = numeric(),
    env_mean_reads   = numeric(),
    ratio            = numeric(),
    n_env_samples    = integer(),
    reads_zeroed     = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (grp in unique_groups) {
    in_group <- !is.na(groups) & groups == grp
    is_pos_g <- is_positive & in_group
    is_env_g <- is_env & in_group
    
    if (sum(is_pos_g) == 0) next
    
    pos_rows <- rows[is_pos_g]
    pos_base <- sub("_rep[0-9]+$", "", pos_rows)
    pos_mat  <- mat[is_pos_g, , drop = FALSE]
    
    # Step 1: reproducibility filter
    base_counts   <- rowsum((pos_mat > 0) * 1L, pos_base)
    max_rep_count <- apply(base_counts, 2, max)
    repro_flagged <- names(max_rep_count)[max_rep_count >= 2]
    
    if (length(repro_flagged) == 0) next
    
    # Step 2: directionality filter
    if (sum(is_env_g) > 0) {
      pos_mean_reads <- colMeans(pos_mat[, repro_flagged, drop = FALSE])
      env_mean_reads <- colMeans(mat[is_env_g, repro_flagged, drop = FALSE])
      ratio          <- ifelse(env_mean_reads == 0, Inf,
                               pos_mean_reads / env_mean_reads)
      truly_flagged  <- repro_flagged[ratio >= min_pos_ctrl_ratio]
      
      cat(sprintf(
        "  Group %s: %d repro-flagged | %d pass directionality (ratio >= %dx) | zeroed in %d env samples\n",
        grp, length(repro_flagged), length(truly_flagged),
        min_pos_ctrl_ratio, sum(is_env_g)
      ))
      
      if (length(truly_flagged) > 0) {
        reads_zeroed <- colSums(mat[is_env_g, truly_flagged, drop = FALSE])
        mat_corrected[is_env_g, truly_flagged] <- 0
        
        removal_log <- rbind(removal_log, data.frame(
          method           = method_name,
          cycling_group    = grp,
          otu              = truly_flagged,
          max_rep_count    = max_rep_count[truly_flagged],
          pos_mean_reads   = round(pos_mean_reads[truly_flagged], 1),
          env_mean_reads   = round(env_mean_reads[truly_flagged], 1),
          ratio            = round(ratio[truly_flagged], 2),
          n_env_samples    = sum(is_env_g),
          reads_zeroed     = reads_zeroed,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  list(mat = mat_corrected[!is_positive, , drop = FALSE],
       log = removal_log)
}

# Run across all methods and collect logs
pos_results <- lapply(names(demult_list), function(nm) {
  filter_positive(demult_list[[nm]], index_list[[nm]],
                  min_pos_ctrl_ratio = 10,
                  method_name = nm)
}) %>% setNames(names(demult_list))

# Extract filtered matrices and combined log
pos_list       <- lapply(pos_results, `[[`, "mat")
pos_filter_log <- do.call(rbind, lapply(pos_results, `[[`, "log"))

# Save log
write.csv(pos_filter_log, path_pos_ctrl_log, row.names = FALSE)

cat("\n=== POSITIVE CONTROL FILTER CHECK ===\n")
for (nm in names(demult_list)) {
  reads_before <- sum(demult_list[[nm]])
  reads_after  <- sum(pos_list[[nm]])
  cat(sprintf(
    "  %-12s | samples: %d -> %d | reads: %s -> %s (lost: %.1f%%)\n",
    nm,
    nrow(demult_list[[nm]]), nrow(pos_list[[nm]]),
    format(reads_before, big.mark = ","),
    format(reads_after,  big.mark = ","),
    (reads_before - reads_after) / reads_before * 100
  ))
}

################################################################################
# 5.  NTC SUBTRACTION
# Per cycling_group x plate: subtract per-OTU maximum NTC read count from
# all samples in the same group only.
################################################################################

ntc_subtraction <- function(mat, index_data) {
  rows         <- rownames(mat)
  sample_types <- index_data[rows, "sample_type"]
  groups       <- get_groups(index_data, rows)
  is_ntc       <- sample_types == "NTC"
  
  if (sum(is_ntc) == 0) {
    warning("No NTC samples found!")
    return(mat)
  }
  
  unique_groups <- unique(groups[!is.na(groups)])
  mat_corrected <- mat
  
  for (grp in unique_groups) {
    in_group  <- !is.na(groups) & groups == grp
    is_ntc_g  <- is_ntc & in_group
    is_samp_g <- !is_ntc & in_group
    
    if (sum(is_ntc_g) == 0) next
    if (sum(is_samp_g) == 0) next
    
    ntc_max <- apply(mat[is_ntc_g, , drop = FALSE], 2,
                     function(y) max(y, na.rm = TRUE))
    ntc_max[is.infinite(ntc_max)] <- 0
    
    mat_corrected[is_samp_g, ] <- pmax(
      sweep(mat[is_samp_g, , drop = FALSE], 2, ntc_max), 0
    )
  }
  
  mat_corrected
}

ntc_list <- lapply(names(pos_list), function(nm) {
  ntc_subtraction(pos_list[[nm]], index_list[[nm]])
}) %>% setNames(names(pos_list))

cat("\n=== NTC SUBTRACTION CHECK ===\n")
for (nm in names(pos_list)) {
  reads_before <- sum(pos_list[[nm]])
  reads_after  <- sum(ntc_list[[nm]])
  cat(sprintf(
    "  %-12s | reads: %s -> %s (removed: %s) | negative values: %d\n",
    nm,
    format(reads_before, big.mark = ","),
    format(reads_after,  big.mark = ","),
    format(reads_before - reads_after, big.mark = ","),
    sum(ntc_list[[nm]] < 0)
  ))
}


################################################################################
# 6.  EXTRACTION BLANK CORRECTION
# Per cycling_group x plate: FilterBlanks correct water/filter samples;
# ExtractBlanks correct biofilm samples; Xblk correct degradation samples.
################################################################################

xblank_correction <- function(mat, index_data) {
  rows         <- rownames(mat)
  groups       <- get_groups(index_data, rows)
  unique_groups <- unique(groups[!is.na(groups)])
  
  corrected_list <- lapply(unique_groups, function(grp) {
    group_rows   <- rows[!is.na(groups) & groups == grp]
    mat_group    <- mat[group_rows, , drop = FALSE]
    sample_types <- index_data[group_rows, "sample_type"]
    substrate    <- index_data[group_rows, "substrate"]
    
    is_filter_blank  <- sample_types == "FilterBlank"
    is_extract_blank <- sample_types == "ExtractBlank"
    is_env           <- !is.na(sample_types) & sample_types == "EnvironmentalSample"
    is_filtered      <- !is.na(substrate) & substrate == "Filter"
    is_degradation   <- sample_types == "DegradationSample"
    
    mat_c <- mat_group
    
    # Filter blanks correct water/filter environmental samples
    to_correct <- is_filtered & is_env
    if (sum(is_filter_blank) > 0 && sum(to_correct) > 0) {
      fb_max <- apply(mat_group[is_filter_blank, , drop = FALSE], 2, max, na.rm = TRUE)
      fb_max[is.infinite(fb_max)] <- 0
      mat_c[to_correct, ] <- pmax(
        sweep(mat_group[to_correct, , drop = FALSE], 2, fb_max), 0
      )
    }
    
    # Extract blanks correct non-filter environmental samples
    to_correct <- is_env & !is_filtered
    if (sum(is_extract_blank) > 0 && sum(to_correct) > 0) {
      eb_max <- apply(mat_group[is_extract_blank, , drop = FALSE], 2, max, na.rm = TRUE)
      eb_max[is.infinite(eb_max)] <- 0
      mat_c[to_correct, ] <- pmax(
        sweep(mat_group[to_correct, , drop = FALSE], 2, eb_max), 0
      )
    }
    
    # Xblk correct degradation samples
    if (sum(is_degradation) > 0) {
      is_xblk <- grepl("Xblk", group_rows, ignore.case = TRUE)
      if (sum(is_xblk) > 0) {
        xb_max <- apply(mat_group[is_xblk, , drop = FALSE], 2, max, na.rm = TRUE)
        xb_max[is.infinite(xb_max)] <- 0
        mat_c[is_degradation, ] <- pmax(
          sweep(mat_group[is_degradation, , drop = FALSE], 2, xb_max), 0
        )
      }
    }
    
    mat_c
  })
  
  mat_final <- do.call(rbind, corrected_list)
  mat_final[rownames(mat), ]
}

xtr_list <- lapply(names(ntc_list), function(nm) {
  xblank_correction(ntc_list[[nm]], index_list[[nm]])
}) %>% setNames(names(ntc_list))

cat("\n=== EXTRACTION BLANK CORRECTION CHECK ===\n")
for (nm in names(ntc_list)) {
  reads_before <- sum(ntc_list[[nm]])
  reads_after  <- sum(xtr_list[[nm]])
  cat(sprintf(
    "  %-12s | reads: %s -> %s (removed: %s) | negative values: %d\n",
    nm,
    format(reads_before, big.mark = ","),
    format(reads_after,  big.mark = ","),
    format(reads_before - reads_after, big.mark = ","),
    sum(xtr_list[[nm]] < 0)
  ))
}


################################################################################
# 7.  KEEP ENVIRONMENTAL SAMPLES ONLY
################################################################################

clean_list <- lapply(names(xtr_list), function(nm) {
  mat        <- xtr_list[[nm]]
  keep_types <- c("EnvironmentalSample", "DegradationSample")
  keep_rows  <- index_list[[nm]][rownames(mat), "sample_type"] %in% keep_types
  mat[keep_rows, , drop = FALSE]
}) %>% setNames(names(xtr_list))


cat("\n=== ENV SAMPLE FILTER CHECK ===\n")
for (nm in names(xtr_list)) {
  remaining_types <- index_list[[nm]][rownames(clean_list[[nm]]), "sample_type"]
  cat(nm, "- remaining sample types:\n")
  print(table(remaining_types, useNA = "ifany"))
}

################################################################################
# 8.  QUALITY FILTER — SIBLING READ DEPTH + OTU COUNT
#
# Within each sibling group (same bio replicate = same base sample name),
# a PCR replicate is removed if EITHER:
#   (a) OTU count < min_OTUs  (hard floor — absolute failure)
#   (b) read count < min_pct_median * median read count of its siblings
#       (relative floor — outlier within otherwise healthy group), unless it has more than 2000 reads !
#   (c) Absolute minimum floor -- remove pcr replicate if total read is <50
#   (d) after the 3 first checks, sinlgetons pcr replicates with less than 2000 reads are also removed
#
# Applied BEFORE rep filters so that rep50/rep75 fractions are calculated
# on clean replicates only.
#
# Removed replicates are logged per method for inspection.
# Bio replicates that lose ALL their PCR replicates are flagged separately
# — these will appear as missing cells in the balance check downstream.
################################################################################

quality_filter_pcr <- function(mat, index_data,
                               min_OTUs            = 3,
                               min_pct_median      = 0.20,
                               min_reads_absolute  = 2000,
                               min_reads_hard      = 50,
                               min_reads_singleton = 2000,
                               method_name         = "unnamed") {
  
  base_names     <- sub("_rep[0-9]+$", "", rownames(mat))
  unique_samples <- unique(base_names)
  
  removed_log <- data.frame(
    sample      = character(), base_name   = character(),
    reads       = integer(),   OTUs        = integer(),
    sibling_med = numeric(),   read_thresh = numeric(),
    reason      = character(), stringsAsFactors = FALSE
  )
  
  keep_rows <- rep(TRUE, nrow(mat))
  names(keep_rows) <- rownames(mat)
  
  # --- First pass: Rules 1, 2, 3 ---
  for (bs in unique_samples) {
    sib_idx     <- which(base_names == bs)
    sib_reads   <- rowSums(mat[sib_idx, , drop = FALSE])
    sib_otus    <- rowSums(mat[sib_idx, , drop = FALSE] > 0)
    sib_median  <- if (length(sib_idx) >= 2) median(sib_reads) else NA
    read_thresh <- if (!is.na(sib_median)) sib_median * min_pct_median else NA
    is_singleton <- length(sib_idx) == 1
    
    for (i in sib_idx) {
      nm      <- rownames(mat)[i]
      r       <- sib_reads[nm]
      o       <- sib_otus[nm]
      reasons <- character(0)
      
      # Rule 1: OTU hard floor — always applies
      if (o < min_OTUs)
        reasons <- c(reasons, "low_OTU")
      
      # Rule 2: relative read floor — skipped for singletons
      if (!is_singleton && !is.na(sib_median) && r < read_thresh && r < min_reads_absolute)
        reasons <- c(reasons, "low_reads_relative")
      
      # Rule 3: absolute read hard floor — always applies
      if (r < min_reads_hard)
        reasons <- c(reasons, "low_reads_hard")
      
      if (length(reasons) > 0) {
        keep_rows[nm] <- FALSE
        removed_log <- rbind(removed_log, data.frame(
          sample      = nm,        base_name   = bs,
          reads       = r,         OTUs        = o,
          sibling_med = if (!is.na(sib_median)) sib_median else NA_real_,
          read_thresh = if (!is.na(read_thresh)) round(read_thresh, 1) else NA_real_,
          reason      = paste(reasons, collapse = "+"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # --- Second pass: Rule 4 - remove singletons below read threshold ---
  # Catches both original singletons and singletons created by Rules 1-3
  for (bs in unique_samples) {
    sib_idx   <- which(base_names == bs)
    surviving <- sib_idx[keep_rows[rownames(mat)[sib_idx]]]
    
    if (length(surviving) == 1) {
      nm <- rownames(mat)[surviving]
      r  <- rowSums(mat[nm, , drop = FALSE])
      o  <- rowSums(mat[nm, , drop = FALSE] > 0)
      
      if (r < min_reads_singleton) {
        keep_rows[nm] <- FALSE
        removed_log <- rbind(removed_log, data.frame(
          sample      = nm,
          base_name   = bs,
          reads       = r,
          OTUs        = o,
          sibling_med = NA_real_,
          read_thresh = NA_real_,
          reason      = "singleton_low_reads",
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  mat_filtered <- mat[keep_rows, , drop = FALSE]
  n_removed    <- sum(!keep_rows)
  n_singleton  <- sum(removed_log$reason == "singleton_low_reads")
  
  cat(sprintf(
    "\n[%s] Quality filter: %d PCR replicates removed\n  Rule 1 (always):    OTUs < %d\n  Rule 2 (relative):  reads < %.0f%% sibling median AND reads < %d\n  Rule 3 (always):    reads < %d\n  Rule 4 (singleton): single surviving rep with reads < %d\n",
    method_name, n_removed, min_OTUs, min_pct_median * 100, min_reads_absolute,
    min_reads_hard, min_reads_singleton
  ))
  
  if (n_removed > 0) {
    cat("  Removed samples:\n")
    print(removed_log[, c("sample","reads","OTUs","sibling_med","read_thresh","reason")],
          row.names = FALSE)
    
    remaining_bases <- sub("_rep[0-9]+$", "", rownames(mat_filtered))
    lost_bioreps    <- setdiff(unique(base_names), unique(remaining_bases))
    if (length(lost_bioreps) > 0) {
      cat(sprintf("  WARNING: %d bio replicate(s) lost ALL PCR replicates:\n", length(lost_bioreps)))
      for (br in lost_bioreps) cat("    -", br, "\n")
    }
  }
  
  return(list(filtered = mat_filtered, log = removed_log))
}




# Apply to all 6 methods
qf_results <- lapply(names(clean_list), function(nm) {
  quality_filter_pcr(clean_list[[nm]], index_list[[nm]],
                     min_OTUs = 3, min_pct_median = 0.20,
                     min_reads_singleton = 2000,
                     method_name = nm)
}) %>% setNames(names(clean_list))

# Extract filtered tables and logs
qf_list <- lapply(qf_results, `[[`, "filtered")
qf_logs <- lapply(qf_results, `[[`, "log")

# Save removal logs
qf_log_all <- do.call(rbind, lapply(names(qf_logs), function(nm) {
  if (nrow(qf_logs[[nm]]) > 0) { qf_logs[[nm]]$method <- nm; qf_logs[[nm]] } else NULL
}))
if (!is.null(qf_log_all) && nrow(qf_log_all) > 0) {
  write.csv(qf_log_all, path_qf_removed_log, row.names = FALSE)
  cat("\nQuality filter log saved.\n")
}



cat("=== QUALITY FILTER DIAGNOSTIC SUMMARY ===\n\n")

qf_summary <- do.call(rbind, lapply(names(qf_list), function(nm) {
  mat        <- qf_list[[nm]]
  base_names <- sub("_rep[0-9]+$", "", rownames(mat))
  
  # Per-sample summary
  reads    <- rowSums(mat)
  otus     <- rowSums(mat > 0)
  n_pcr    <- table(base_names)
  
  data.frame(
    method         = nm,
    n_pcr_reps     = nrow(mat),
    n_bio_reps     = length(unique(base_names)),
    min_reads      = min(reads),
    median_reads   = median(reads),
    max_reads      = max(reads),
    min_otus       = min(otus),
    median_otus    = median(otus),
    max_otus       = max(otus),
    n_singletons   = sum(n_pcr == 1),
    n_pairs        = sum(n_pcr == 2),
    n_triplets     = sum(n_pcr == 3),
    n_quadruplets  = sum(n_pcr == 4)
  )
}))

print(kable(qf_summary, format = "simple"))

# Also flag any remaining low-depth or low-OTU samples after filtering
cat("\n--- Remaining samples with reads < 1500 ---\n")
low_after <- do.call(rbind, lapply(names(qf_list), function(nm) {
  mat   <- qf_list[[nm]]
  reads <- rowSums(mat)
  otus  <- rowSums(mat > 0)
  low   <- reads < 1500
  if (sum(low) == 0) return(NULL)
  data.frame(
    method = nm,
    sample = rownames(mat)[low],
    reads  = reads[low],
    otus   = otus[low]
  )
}))

if (is.null(low_after) || nrow(low_after) == 0) {
  cat("✓ No remaining samples with reads < 1500\n")
} else {
  print(low_after)
}




################################################################################
# 8b.  REMOVE ZERO-READ SAMPLES AND ZERO-ABUNDANCE OTUs
# Applied immediately after quality filtering, before rep filters.
# - Zero-read samples: PCR replicates that lost all reads (e.g. from blank
#   correction) but weren't caught by the quality filter (singletons, or
#   samples where all OTUs were subtracted to 0)
# - Zero-abundance OTUs: OTUs with no reads remaining across all samples
#   after quality filter removed bad replicates
################################################################################

qf_list <- lapply(names(qf_list), function(nm) {
  mat <- qf_list[[nm]]
  
  # Remove zero-read samples
  zero_samples <- rownames(mat)[rowSums(mat) == 0]
  if (length(zero_samples) > 0) {
    cat(nm, "- removing", length(zero_samples), "zero-read samples after quality filter\n")
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  }
  
  # Remove zero-abundance OTUs
  zero_otus <- sum(colSums(mat) == 0)
  if (zero_otus > 0) {
    cat(nm, "- removing", zero_otus, "zero-abundance OTUs after quality filter\n")
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
  }
  
  mat
}) %>% setNames(names(qf_list))

# Save the quality-filtered (pre-rep-filter) tables — useful for diagnostics
qf_save_paths <- list(
  nP.lulu  = path_qf_nP_lulu,
  psP.lulu = path_qf_psP_lulu,
  fP.lulu  = path_qf_fP_lulu,
  nP.sw    = path_qf_nP_sw,
  psP.sw   = path_qf_psP_sw,
  fP.sw    = path_qf_fP_sw
)
for (nm in names(qf_list)) saveRDS(qf_list[[nm]], qf_save_paths[[nm]])
cat("Quality-filtered tables saved.\n")


################################################################################
# 9.  REPLICATE FILTERS (rep50, rep75)
# Applied to quality-filtered clean tables.
# OTUs zeroed out if present in fewer than min_fraction of PCR replicates.
################################################################################

replicate_filter <- function(mat, min_fraction) {
  base_names     <- sub("_rep[0-9]+$", "", rownames(mat))
  unique_samples <- unique(base_names)
  mat_filtered   <- mat
  
  for (sample in unique_samples) {
    sample_rows    <- which(base_names == sample)
    n_pcr_reps     <- length(sample_rows)
    min_reps       <- ceiling(n_pcr_reps * min_fraction)
    presence_counts <- colSums(mat[sample_rows, , drop = FALSE] > 0)
    mat_filtered[sample_rows, presence_counts < min_reps] <- 0
  }
  return(mat_filtered)
}

rep50_list <- lapply(qf_list, replicate_filter, min_fraction = 0.50)
rep75_list <- lapply(qf_list, replicate_filter, min_fraction = 0.75)

################################################################################
# 10.  COLLAPSE PCR REPLICATES
# Sum reads across PCR replicates per bio replicate (base sample name).
################################################################################

collapse_reps <- function(mat) {
  base_names     <- sub("_rep[0-9]+$", "", rownames(mat))
  unique_samples <- unique(base_names)
  collapsed <- do.call(rbind, lapply(unique_samples, function(sample) {
    colSums(mat[base_names == sample, , drop = FALSE])
  }))
  rownames(collapsed) <- unique_samples
  return(collapsed)
}

collap_norep_list <- lapply(qf_list,   collapse_reps)
collap_rep50_list <- lapply(rep50_list, collapse_reps)
collap_rep75_list <- lapply(rep75_list, collapse_reps)


################################################################################
# 11b.  REMOVE ZERO-ABUNDANCE OTUs AND ZERO-SUM SAMPLES
# Applied after each transformation that can introduce new zeros:
#   - Rep filters zero out OTU values → new zero-sum OTUs possible in rep lists
#   - Collapsing sums across (potentially rep-filtered) replicates → new
#     zero-sum samples and OTUs possible in collapsed lists
# qf_list is already clean from Step 8b and does not need sample removal here.
################################################################################

remove_zero_otus <- function(mat) mat[, colSums(mat) > 0, drop = FALSE]
remove_zero_both <- function(mat) {
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]  # zero-sum samples
  mat <- mat[, colSums(mat) > 0, drop = FALSE]  # zero-sum OTUs
  mat
}

# qf_list: zero-OTU removal only (samples already clean from Step 8b)
qf_list <- lapply(qf_list, remove_zero_otus)

# rep-filtered lists: rep filter may have created new zero-sum OTUs
rep50_list <- lapply(rep50_list, remove_zero_both)
rep75_list <- lapply(rep75_list, remove_zero_both)

# collapsed lists: collapsing can produce both zero-sum samples and OTUs
collap_norep_list <- lapply(collap_norep_list, remove_zero_both)
collap_rep50_list <- lapply(collap_rep50_list, remove_zero_both)
collap_rep75_list <- lapply(collap_rep75_list, remove_zero_both)

cat("\n=== STEP 11 ZERO REMOVAL SUMMARY ===\n")
all_lists <- list(
  qf           = qf_list,
  rep50        = rep50_list,
  rep75        = rep75_list,
  collap_norep = collap_norep_list,
  collap_rep50 = collap_rep50_list,
  collap_rep75 = collap_rep75_list
)

for (list_name in names(all_lists)) {
  for (nm in names(all_lists[[list_name]])) {
    mat <- all_lists[[list_name]][[nm]]
    cat(sprintf(
      "  %-20s | %-12s | %d samples | %d OTUs | %s total reads | median %s reads/sample | min %s reads/sample\n",
      list_name, nm,
      nrow(mat),
      ncol(mat),
      format(sum(mat),                   big.mark = ","),
      format(round(median(rowSums(mat))), big.mark = ","),
      format(min(rowSums(mat)),           big.mark = ",")
    ))
  }
}




################################################################################
# 11C.  REMOVE genuinely low-read samples
# remove biological replicates with 
# insufficient total reads after collapsing PCR replicates.
# This catches samples where ALL PCR replicates were consistently low,
# which the PCR-level relative filter cannot detect.
# Note: low OTU count with HIGH reads is kept — that is real biology.
################################################################################

post_collapse_quality_filter <- function(mat, min_reads = 1000, 
                                         min_otus = 3,
                                         method_name = "unnamed") {
  reads      <- rowSums(mat)
  otus       <- rowSums(mat > 0)
  fail_reads <- reads < min_reads
  fail_otus  <- otus  < min_otus
  fail_any   <- fail_reads | fail_otus
  
  removed_df <- data.frame(
    sample     = rownames(mat)[fail_any],
    reads      = reads[fail_any],
    n_otus     = otus[fail_any],
    fail_reads = fail_reads[fail_any],
    fail_otus  = fail_otus[fail_any]
  )
  
  cat(sprintf(
    "\n[%s] Post-collapse filter: %d / %d biological replicates removed\n  reads < %d: %d samples\n  OTUs < %d:  %d samples\n",
    method_name, sum(fail_any), nrow(mat),
    min_reads, sum(fail_reads),
    min_otus,  sum(fail_otus)
  ))
  if (nrow(removed_df) > 0) print(removed_df)
  
  mat[!fail_any, , drop = FALSE]
}

# Apply to all three collapsed versions
collap_norep_list <- lapply(names(collap_norep_list), function(nm) {
  post_collapse_quality_filter(collap_norep_list[[nm]],
                               min_reads = 1000, min_otus = 3,
                               method_name = paste(nm, "norep"))
}) %>% setNames(names(collap_norep_list))

collap_rep50_list <- lapply(names(collap_rep50_list), function(nm) {
  post_collapse_quality_filter(collap_rep50_list[[nm]],
                               min_reads = 1000, min_otus = 3,
                               method_name = paste(nm, "rep50"))
}) %>% setNames(names(collap_rep50_list))

collap_rep75_list <- lapply(names(collap_rep75_list), function(nm) {
  post_collapse_quality_filter(collap_rep75_list[[nm]],
                               min_reads = 1000, min_otus = 3,
                               method_name = paste(nm, "rep75"))
}) %>% setNames(names(collap_rep75_list))


cat("=== POST-COLLAPSE QUALITY FILTER DIAGNOSTIC ===\n\n")

all_collapsed <- list(
  norep = collap_norep_list,
  rep50 = collap_rep50_list,
  rep75 = collap_rep75_list
)

collapse_summary <- do.call(rbind, lapply(names(all_collapsed), function(version) {
  do.call(rbind, lapply(names(all_collapsed[[version]]), function(nm) {
    mat   <- all_collapsed[[version]][[nm]]
    reads <- rowSums(mat)
    otus  <- rowSums(mat > 0)
    data.frame(
      version      = version,
      method       = nm,
      n_samples    = nrow(mat),
      min_reads    = min(reads),
      median_reads = median(reads),
      max_reads    = max(reads),
      min_otus     = min(otus),
      median_otus  = median(otus),
      max_otus     = max(otus)
    )
  }))
}))

print(kable(collapse_summary, format = "simple"))

################################################################################
# 12.  SAVE OUTPUT TABLES
################################################################################

# Helper: save one OTU table with a given naming template
save_otu <- function(mat, lulu_sw, nP_psP_fP, suffix) {
  fname <- sprintf("%s/%s_%s_%s.rds", out_dir, lulu_sw, nP_psP_fP, suffix)
  saveRDS(mat, fname)
  cat("Saved ->", basename(fname), "\n")
}

# Method key: list names → (clustering, pool) labels
method_key <- list(
  nP.lulu  = c("LULU",  "nP"),
  psP.lulu = c("LULU",  "psP"),
  fP.lulu  = c("LULU",  "fP"),
  nP.sw    = c("SWARM", "nP"),
  psP.sw   = c("SWARM", "psP"),
  fP.sw    = c("SWARM", "fP")
)

suffixes <- list(
  norepfilt          = qf_list,
  repfilt50          = rep50_list,
  repfilt75          = rep75_list,
  collap_norepfilt   = collap_norep_list,
  CollapsedRepFilt50 = collap_rep50_list,
  CollapsedRepFilt75 = collap_rep75_list
)

cat("\n### SAVING OTU TABLES ###\n\n")
for (suffix in names(suffixes)) {
  for (nm in names(method_key)) {
    lulu_sw  <- method_key[[nm]][1]
    pool_tag <- method_key[[nm]][2]
    save_otu(suffixes[[suffix]][[nm]], lulu_sw, pool_tag, suffix)
  }
}




################################################################################
# 13.  DIAGNOSTICS — DECONTAM LIBRARY SIZE PLOTS
################################################################################

cat("\n### GENERATING DECONTAM LIBRARY SIZE PLOTS ###\n\n")

generate_decontam_plots <- function(index_data, method_name) {
  cycling_groups <- unique(index_data$cycling_group)
  cycling_groups <- cycling_groups[!is.na(cycling_groups)]
  
  for (group in cycling_groups) {
    df <- index_data[index_data$cycling_group == group, ]
    df <- df[order(df$totseq), ]
    df$index <- seq(nrow(df))
    df$alpha  <- ifelse(df$sample_type == "EnvironmentalSample" |
                          df$sample_type == "DegradationSample" |
                          df$sample_type == "PositiveControl", 0.8, 1)
    
    p <- ggplot(df, aes(x = index, y = log(totseq + 1), color = sample_type)) +
      geom_point(aes(alpha = alpha)) +
      facet_grid(cols = vars(PCR_rep)) +
      ggtitle(paste("Log library size -", method_name, "- cycling group", group))
    
    pdf(file.path(dir_decontam_plots, paste0(method_name, "_cycling_group_", group, ".pdf")),
        width = 7, height = 5)
    print(p)
    dev.off()
  }
}

for (nm in names(index_list)) generate_decontam_plots(index_list[[nm]], nm)




################################################################################
# 13b.  TAG JUMPING ASSESSMENT
################################################################################

cat("\n### TAG JUMPING ASSESSMENT ###\n\n")

tag_jump_rates <- do.call(rbind, lapply(names(index_list), function(nm) {
  idx <- index_list[[nm]]
  
  by_plate <- lapply(split(idx, idx$plate), function(x) {
    total    <- sum(x$totseq, na.rm = TRUE)
    blanks   <- sum(x[x$sample_type == "Blank", "totseq"], na.rm = TRUE)
    rate     <- if (total > 0) blanks / total else NA
    data.frame(method = nm, group_type = "plate", group = unique(x$plate)[1],
               total_reads = total, blank_reads = blanks, jump_rate = rate)
  })
  
  by_cyc <- lapply(split(idx, idx$cycling_group), function(x) {
    total  <- sum(x$totseq, na.rm = TRUE)
    blanks <- sum(x[x$sample_type == "Blank", "totseq"], na.rm = TRUE)
    rate   <- if (total > 0) blanks / total else NA
    data.frame(method = nm, group_type = "cycling_group", group = unique(x$cycling_group)[1],
               total_reads = total, blank_reads = blanks, jump_rate = rate)
  })
  
  do.call(rbind, c(by_plate, by_cyc))
}))

tag_jump_rates$per_100k <- round(tag_jump_rates$jump_rate * 1e5, 2)
print(tag_jump_rates[order(tag_jump_rates$jump_rate, decreasing = TRUE), ])

write.csv(tag_jump_rates, path_tag_jump_csv, row.names = FALSE)

################################################################################
# 13c.  FILTERING SUMMARY TABLE
################################################################################

cat("\n### FILTERING SUMMARY ###\n\n")

get_stats <- function(mat) {
  list(samples  = nrow(mat),
       reads    = sum(mat),
       otus     = sum(colSums(mat) > 0))
}

filter_summary <- do.call(rbind, lapply(names(method_key), function(nm) {
  lulu_sw  <- method_key[[nm]][1]
  pool_tag <- method_key[[nm]][2]
  method_label <- paste0(pool_tag, "_", lulu_sw)
  
  steps <- list(
    "1_Demultiplexed"        = demult_list[[nm]],
    "2_Positive_filtered"    = pos_list[[nm]],
    "3_NTC_subtracted"       = ntc_list[[nm]],
    "4_Extraction_corrected" = xtr_list[[nm]],
    "5_EnvSamples_only"      = clean_list[[nm]],
    "6_QualityFiltered"      = qf_list[[nm]],
    "7_RepFilter_50pct"      = rep50_list[[nm]],
    "8_RepFilter_75pct"      = rep75_list[[nm]],
    "9_Collapsed_NoFilter"   = collap_norep_list[[nm]],
    "10_Collapsed_50pct"     = collap_rep50_list[[nm]],
    "11_Collapsed_75pct"     = collap_rep75_list[[nm]]
  )
  
  stats_df <- do.call(rbind, lapply(names(steps), function(step) {
    s <- get_stats(steps[[step]])
    data.frame(method = method_label, step = step,
               samples = s$samples, reads = s$reads, otus = s$otus,
               stringsAsFactors = FALSE)
  }))
  
  stats_df$reads_retained_pct <- round(stats_df$reads / stats_df$reads[1] * 100, 2)
  stats_df$otus_retained_pct  <- round(stats_df$otus  / stats_df$otus[1]  * 100, 2)
  stats_df
}))

print(filter_summary)
write.csv(filter_summary, path_filter_summary_csv, row.names = FALSE)

saveRDS(list(
  demult         = demult_list,
  pos_filtered   = pos_list,
  ntc_filtered   = ntc_list,
  xtr_filtered   = xtr_list,
  env_only       = clean_list,
  quality_filt   = qf_list,
  rep50          = rep50_list,
  rep75          = rep75_list,
  collap_norep   = collap_norep_list,
  collap_rep50   = collap_rep50_list,
  collap_rep75   = collap_rep75_list,
  index          = index_list,
  qf_logs        = qf_logs,
  filter_summary = filter_summary
), path_filter_tracking)

save.image(path_rdata)

cat("\n### PIPELINE COMPLETE ###\n")
cat("Quality filter thresholds: min OTUs =", 3,
    "| min read depth = 20% of sibling median\n")
cat("Output tables saved to:", out_dir, "\n")





################################################################################
# VISUALISATION — READS AND OTU LOSS THROUGHOUT THE PIPELINE
################################################################################

plot_dir <- dir_read_tracking

# ---- Prepare data ----

# Clean step labels for axis (strip the number prefix)
filter_summary_plot <- filter_summary %>%
  mutate(
    step_label = sub("^[0-9]+_", "", step),
    step_label = factor(step_label, levels = unique(step_label)),
    # Separate clustering and pooling for faceting
    clustering = ifelse(grepl("LULU", method), "LULU", "SWARM"),
    pooling    = sub("_(LULU|SWARM)$", "", method)
  )

# Pivot to long format for reads and OTUs together
filter_long <- filter_summary_plot %>%
  select(method, clustering, pooling, step_label,
         reads_retained_pct, otus_retained_pct) %>%
  pivot_longer(
    cols      = c(reads_retained_pct, otus_retained_pct),
    names_to  = "metric",
    values_to = "pct_retained"
  ) %>%
  mutate(metric = case_when(
    metric == "reads_retained_pct" ~ "Reads",
    metric == "otus_retained_pct"  ~ "OTUs"
  ))

# ---- Plot 1: Retention % through pipeline, faceted by pooling x clustering ----

p1 <- ggplot(filter_long,
             aes(x = step_label, y = pct_retained,
                 color = metric, group = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(rows = vars(pooling), cols = vars(clustering)) +
  scale_color_manual(values = c("Reads" = "#2166ac", "OTUs" = "#d6604d")) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Read and OTU retention through the filtering pipeline",
    subtitle = "Each step shown as % of starting demultiplexed counts",
    x        = NULL,
    y        = "% retained",
    color    = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_dir, "pipeline_retention_by_method.pdf"),
       p1, width = 12, height = 10)
print(p1)


# ---- Plot 2 revised: absolute reads and OTUs, line plot, no rep filter steps ----

steps_no_repfilt <- c("Demultiplexed", "Positive_filtered", "NTC_subtracted",
                      "Extraction_corrected", "EnvSamples_only", "QualityFiltered",
                      "Collapsed_NoFilter", "Collapsed_50pct", "Collapsed_75pct")

filter_line_abs <- filter_summary_plot %>%
  filter(step_label %in% steps_no_repfilt) %>%
  mutate(step_label = factor(step_label, levels = steps_no_repfilt)) %>%
  pivot_longer(
    cols      = c(reads, otus),
    names_to  = "metric",
    values_to = "count"
  ) %>%
  mutate(metric = case_when(
    metric == "reads" ~ "Reads",
    metric == "otus"  ~ "OTUs"
  ))

p2<- ggplot(filter_line_abs,
            aes(x = step_label, y = count,
                color = metric, group = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(rows = vars(metric, pooling), cols = vars(clustering),
             scales = "free_y")+
  scale_color_manual(values = c("Reads" = "#2166ac", "OTUs" = "#d6604d")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Read and OTU counts through the filtering pipeline",
    subtitle = "Absolute counts — control filtering and quality filter only",
    x        = NULL,
    y        = "Count",
    color    = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_dir, "pipeline_counts_lineplot_no_repfilt.pdf"),
       p2, width = 12, height = 10)
print(p2)


# ---- Plot 3: LULU vs SWARM comparison at key endpoints ----
# Useful for seeing if clustering method changes how much is lost

endpoints <- c("EnvSamples_only", "QualityFiltered",
               "RepFilter_50pct", "RepFilter_75pct",
               "Collapsed_NoFilter", "Collapsed_50pct", "Collapsed_75pct")

filter_endpoints <- filter_summary_plot %>%
  filter(step_label %in% endpoints) %>%
  mutate(step_label = factor(step_label, levels = endpoints))

p3 <- ggplot(filter_endpoints,
             aes(x = step_label, y = otus_retained_pct,
                 color = clustering, shape = pooling, group = interaction(clustering, pooling))) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("LULU" = "#1b7837", "SWARM" = "#762a83")) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "OTU retention at key endpoints — LULU vs SWARM",
    subtitle = "Comparing clustering method and pooling strategy",
    x        = NULL,
    y        = "OTUs retained (%)",
    color    = "Clustering",
    shape    = "Pooling"
  ) +
  theme_bw(base_size = 11) +
  
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_dir, "pipeline_lulu_vs_swarm_endpoints.pdf"),
       p3, width = 9, height = 5)
print(p3)




