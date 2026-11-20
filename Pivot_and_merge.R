# ============================================================
# CMI-Flu OMICS Project - Group 2
# Pivot_and_merge.R
#
# Purpose:
# Prepare shared datasets for downstream analyses by:
# 1. Selecting Day 7 plasmablast frequency (outcome)
# 2. Selecting Day 0 PBMC RNA-seq (predictors)
# 3. Pivoting RNA-seq from long to wide format
# 4. Merging expression data with participant information
# 5. Saving separate raw-count, TPM, and metadata datasets
#
# This script does NOT normalize, filter genes, or fit models.
#
# Outputs:
#   data/processed/analysis_counts.rds --> Raw RNA-seq counts + 
#                                 Day 7 plasmablast frequency + demographics 
#                                 (Andrew) for TMM normalization; Tony for limma-voom
#   data/processed/analysis_tpm.rds --> Pre-normalized TPM gene expression + Day 7
#                                     plasmablast frequency + demographics 
#                                       (Genevieve for elastic net)
#   data/processed/analysis_metadata.rds --> Participant info, vaccination history, and Day 7
#                                       plasmablast frequency, without gene expression
#                                       (everyone for Table 1 and exploratory analysis)
#
# Run from the RStudio project root.
# ============================================================


# 1. Load packages and data -------------------------------

library(tidyverse)
library(edgeR)
library(sva)

participants <- readRDS("data/participants_cleaned.rds")
rnaseq <- readRDS("data/rnaseq.rds")
vaccine_hist <- readRDS("data/vaccine_hist.rds")

# Check dimensions of original datasets
dim(participants)
# 29, 19
dim(rnaseq)
# 2289028, 9
dim(vaccine_hist)
# 204, 3


# 2. PREPARE DAY 7 PLASMABLAST OUTCOME -------------------------

# The participants_cleaned dataset contains the Day 7
# antibody-secreting cell (ASC) frequency for 29 participants.
#
# Retain demographic variables and previous vaccination history.

outcome <- participants %>%
  filter(
    timepoint == 7,
    name == "Antibody-secreting_cells_(ASC)"
  ) %>%
  select(
    participant_id,
    subject,
    biological_sex,
    age,
    plasmablast_frequency = value,
    n_prior_vax
  ) %>%
  mutate(
    transf_plasmablast_frequency = log(plasmablast_frequency)
  )


# Check outcome observations
dim(outcome)
# 29, 6
summary(outcome$plasmablast_frequency)

# Confirm each participant has only one outcome measurement
stopifnot(!anyDuplicated(outcome$participant_id))

# Ensure outcome is numeric and nonmissing
stopifnot(
  is.numeric(outcome$plasmablast_frequency),
  !anyNA(outcome$plasmablast_frequency)
)

# 3. CHECK VACCINATION HISTORY --------------------------------

# participants_cleaned already contains n_prior_vax,
# representing the number of previously recorded vaccinations.
#
# I previously confirmed that n_prior_vax agrees with the
# number of "Yes" vaccination records for all 27 participants
# with vaccination-history records.
#
# Two participants (ID_443 and ID_683) have n_prior_vax = 0
# but no vaccination-history records. We preserve this
# distinction rather than assuming confirmed non-vaccination.

outcome <- outcome %>%
  mutate(
    vax_history_available = subject %in% vaccine_hist$subject
  )

# Check vaccination history availability
table(outcome$vax_history_available)
# False = 2, True = 27

# The detailed vaccine-history table is not merged because
# it contains multiple yearly records per participant.
# Directly joining it would create duplicate participant rows.

# 4. SELECT BASELINE RNA-SEQ ----------------------------------

# Keep only Day 0 PBMC RNA-seq observations.
#
# counts = raw gene counts
# tpm    = transcripts per million (already normalized
#          for gene length and sequencing depth)
#
# We retain BOTH to allow different preprocessing approaches.

rnaseq_d0 <- rnaseq %>%
  filter(
    timepoint == 0,
    material == "PBMCs",
    study_accession == "2024_UGA"
  ) %>%
  select(
    participant_id,
    ensembl_gene_id,
    counts,
    tpm
  )

# Check number of unique participants and genes
n_distinct(rnaseq_d0$participant_id)
# 29
n_distinct(rnaseq_d0$ensembl_gene_id)
# 78932

# Check whether expression values are missing
sum(is.na(rnaseq_d0$counts))
# 0 
sum(is.na(rnaseq_d0$tpm))
# 0 

# 5. CHECK FOR DUPLICATE GENE MEASUREMENTS ---------------------

# Each participant should have one Day 0 measurement per gene.
#
# Do not automatically average or remove duplicates:
# they could represent repeated samples or technical issues.

duplicates <- rnaseq_d0 %>%
  count(participant_id, ensembl_gene_id) %>%
  filter(n > 1)

print(duplicates)
# 0 rows

if (nrow(duplicates) > 0) {
  stop("Duplicate participant-gene combinations detected.")
}

# Check gene identifiers
stopifnot(
  !anyNA(rnaseq_d0$participant_id),
  !anyNA(rnaseq_d0$ensembl_gene_id)
)

# 6. CHECK PARTICIPANT OVERLAP --------------------------------

# Determine which participants have both:
# - Day 0 RNA-seq
# - Day 7 plasmablast frequency

common_ids <- intersect(
  outcome$participant_id,
  rnaseq_d0$participant_id
)

cat("Participants with Day 7 outcome:",
    n_distinct(outcome$participant_id), "\n")
# 29

cat("Participants with Day 0 RNA-seq:",
    n_distinct(rnaseq_d0$participant_id), "\n")
# 29

cat("Participants with both measurements:",
    length(common_ids), "\n")
# 29

# Check which outcome participants lack RNA-seq
missing_rnaseq <- setdiff(
  outcome$participant_id,
  rnaseq_d0$participant_id
)

print(missing_rnaseq)
# none 

if (length(common_ids) == 0) {
  stop("No matching participants found.")
}

# The expected analysis population is 29 participants.
# Report a warning rather than silently assuming all match.

if (length(common_ids) != 29) {
  warning("Matched population differs from expected n = 29.")
}

# 7. PIVOT RNA-SEQ INTO WIDE FORMAT AND RENORMALIZE ----------------------------

# Currently RNA-seq is in long format:
# participant_id | ensembl_gene_id | counts | tpm
#
# For modeling, we need:
# participant_id | Gene1 | Gene2 | Gene3 | ... 
#
# We create two separate matrices to avoid mixing
# raw counts and TPM values in the same model.

# A. Raw-count matrix

rnaseq_counts_wide <- rnaseq_d0 %>%
  select(participant_id, ensembl_gene_id, counts) %>%
  pivot_wider(
    names_from = ensembl_gene_id,
    values_from = counts
  )

# B. TPM matrix

rnaseq_tpm_wide <- rnaseq_d0 %>%
  select(participant_id, ensembl_gene_id, tpm) %>%
  pivot_wider(
    names_from = ensembl_gene_id,
    values_from = tpm
  )

# C. Add n_prior_vax and vax_history_available to the RNA-seq count matrix

rnaseq_counts_tmm_wide <- rnaseq_counts_wide %>%
  left_join(
    outcome %>%
      select(participant_id, n_prior_vax),
    by = "participant_id"
  )


# Check dimensions
dim(rnaseq_counts_wide)
# 29, 78933
dim(rnaseq_tpm_wide)
# 29, 78933

# Each row should represent one participant.
# Each gene should be a separate column.

gene_cols <- c(2:ncol(rnaseq_counts_tmm_wide)) # columns with gene counts

# Normalizing raw-count matrix with TMM method

dge_tmm <- DGEList(t(rnaseq_counts_tmm_wide[,gene_cols]))

keep <- filterByExpr(dge_tmm)

dge_tmm <- dge_tmm[keep, , keep.lib.sizes=FALSE] %>%
  calcNormFactors(method = "TMM")

counts_log2_cpm_tmm <- cpm(dge_tmm, log = TRUE) %>%
  t() %>%
  as.data.frame()

counts_log2_cpm_tmm$participant_id = rnaseq_counts_wide$participant_id


# 8. MERGE WITH PARTICIPANT INFORMATION -----------------------

# inner_join keeps only participants present in BOTH datasets.
#
# This ensures that all included participants have
# baseline RNA-seq and Day 7 outcome data.


# A. Raw-count analysis dataset

analysis_counts <- outcome %>%
  inner_join(
    rnaseq_counts_wide,
    by = "participant_id"
  )


# B. TPM analysis dataset

analysis_tpm <- outcome %>%
  inner_join(
    rnaseq_tpm_wide,
    by = "participant_id"
  )

# C. TMM analysis dataset

analysis_tmm <- outcome %>%
  inner_join(
    counts_log2_cpm_tmm,
    by = "participant_id"
  )

# D. Participant metadata only

analysis_metadata <- outcome %>%
  filter(participant_id %in% common_ids)



# 9. VERIFY MATCHING AND DATA STRUCTURE -----------------------

# Ensure one row per participant
stopifnot(
  !anyDuplicated(analysis_counts$participant_id),
  !anyDuplicated(analysis_tpm$participant_id),
  !anyDuplicated(analysis_tmm$participant_id)
)

# Check that both datasets contain the same participants
stopifnot(
  setequal(
    analysis_counts$participant_id,
    analysis_tpm$participant_id
  )
)

stopifnot(
  setequal(
    analysis_counts$participant_id,
    analysis_tmm$participant_id
  )
)


# Check that both expression matrices contain the same genes
gene_cols <- 
setdiff(
    names(rnaseq_counts_wide),
    "participant_id"
  )

stopifnot(
  setequal(
    gene_cols,
    setdiff(names(counts_log2_cpm_tmm), "participant_id")
  )
)

# Align participant ordering across all three outputs
analysis_counts <- analysis_counts %>%
  arrange(participant_id)

analysis_tpm <- analysis_tpm %>%
  arrange(participant_id)

analysis_tmm <- analysis_tmm %>%
  arrange(participant_id)

analysis_metadata <- analysis_metadata %>%
  arrange(participant_id)

# Confirm identical participant ordering
stopifnot(
  identical(
    analysis_counts$participant_id,
    analysis_tpm$participant_id
  ),
  identical(
    analysis_counts$participant_id,
    analysis_metadata$participant_id
  ),
  identical(
    analysis_counts$participant_id,
    analysis_tmm$participant_id
  )
)

# Confirm identical outcome values
stopifnot(
  identical(
    analysis_counts$plasmablast_frequency,
    analysis_tpm$plasmablast_frequency
  ),
  identical(
    analysis_counts$plasmablast_frequency,
    analysis_tmm$plasmablast_frequency
  ),
)

# Report missing gene-expression values
cat("Missing raw-count values:",
    sum(is.na(analysis_counts[, gene_cols])), "\n")
# 0

cat("Missing TPM values:",
    sum(is.na(analysis_tpm[, gene_cols])), "\n")
# 0

cat("Missing TPM values:",
    sum(is.na(analysis_tmm[, gene_cols])), "\n")
# 0

# Missing expression entries are NOT automatically set to zero.
# Their meaning should be investigated before modeling.


# 10. FINAL DATASET SUMMARY -----------------------------------

cat("\n--- Final analysis datasets ---\n")

cat("Participants:",
    nrow(analysis_metadata), "\n")
# 29

cat("Genes:",
    length(gene_cols), "\n")
# 78932

cat("Raw-count dataset dimensions:",
    dim(analysis_counts), "\n")
# 29, 78939

cat("TPM dataset dimensions:",
    dim(analysis_tpm), "\n")
# 29, 78939

cat("TPM dataset dimensions:",
    dim(analysis_tmm), "\n")
# 29, 78939

# Outcome summary
summary(analysis_metadata$plasmablast_frequency)

# Check demographic variables
analysis_metadata %>%
  select(
    participant_id,
    age,
    biological_sex,
    n_prior_vax,
    vax_history_available,
    plasmablast_frequency
  ) %>%
  head()


# 11. SAVE SHARED DATASETS ------------------------------------

# Create output directory if needed
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  analysis_counts,
  "data/processed/analysis_counts.rds"
)

saveRDS(
  analysis_tpm,
  "data/processed/analysis_tpm.rds"
)

saveRDS(
  analysis_tpm,
  "data/processed/analysis_tmm.rds"
)

saveRDS(
  dge_tmm,
  "data/processed/dge_tmm.rds"
)

saveRDS(
  analysis_metadata,
  "data/processed/analysis_metadata.rds"
)

cat("\nData preparation completed successfully.\n")



