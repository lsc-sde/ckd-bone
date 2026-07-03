# ============================================================================
# FEASIBILITY 2: Frailty Fractures in Patients with CKD (eGFR < threshold)
# ============================================================================
# Round 2 grant application feasibility analysis.
#
# Study definition:
#   Base cohort : Patients aged ≥50 with frailty fracture in study period
#   Key question: Of those, who has eGFR within ±4 weeks of fracture?
#                 Of those, who has eGFR < threshold (default 44)?
#   Then:         Characterise creatinine history, demographics, missingness
#
# Cohort flow (nested):
#   N1: Fracture cohort (age ≥50, study period)
#   N2: Of N1 — have eGFR within ±egfr_window_days of fracture
#   N3: Of N2 — eGFR < egfr_threshold (ANALYSIS COHORT)
#   N4: Of N3 — have ≥1 creatinine in 5yr pre-fracture
#   N5: Of N3 — eGFR or CrCl < 30
#
# Outputs:
#   - Attrition waterfall (N1 → N2 → N3 → N4 → N5)
#   - Creatinine characterisation (median count, median value pre-fracture)
#   - Demographic missingness (age, sex, ethnicity, IMD, BMI)
#   - Comorbidity prevalence
#   - Frailty score characterisation (Rockwood/CFS)
#   - Exported codelists CSV for partner sites
#   - Consolidated feasibility report (HTML)
# ============================================================================

# ----------------------------------------------------------------------------
# 0. LIBRARIES & PARAMETERS
# ----------------------------------------------------------------------------
library(CDMConnector)
library(CodelistGenerator)
library(CohortConstructor)
library(PatientProfiles)
library(CohortCharacteristics)
library(omopgenerics)
library(PhenotypeR)
library(OmopSketch)
library(visOmopResults)

library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(glue)
library(here)
library(flextable)
library(gt)
library(purrr)

# --- Study parameters --------------------------------------------------------
study_start            <- as.Date("2022-01-01")
study_end              <- as.Date("2023-12-31")
min_age                <- 50
egfr_threshold         <- 44        # eGFR cut-off for analysis cohort
egfr_window_days       <- 28        # ±4 weeks around fracture date
lookback_creatinine_days <- 1825    # 5 years
suppression_threshold  <- 5         # Cell suppression: counts < this are masked (i.e. 1-4)

# --- Toggle: generate codelists or use cached? ------------------------------
# Set FALSE after first run to skip expensive codelist generation.
# NOTE: Set TRUE to regenerate after adding thyrotoxicosis keywords and calcium codelist.
generate_codelists <- FALSE

# --- Known concept IDs for this OMOP instance --------------------------------
egfr_known     <- 40771922L
creat_known    <- 3020564L
frailty_known  <- 40483383L   # LTHT frailty score concept

# --- Cell suppression helper --------------------------------------------------
# Masks counts below suppression_threshold to prevent disclosure of small cells.
# Returns "<5" (or appropriate label) for suppressed values, otherwise the count as string.
suppress_count <- function(n, threshold = suppression_threshold) {
  if_else(n > 0 & n < threshold, paste0("<", threshold), as.character(n))
}
suppress_pct <- function(n, pct, threshold = suppression_threshold) {
  if_else(n > 0 & n < threshold, "-", as.character(pct))
}

# --- Client-side type helper --------------------------------------------------
patch_int64 <- function(df) {
  df |> mutate(across(where(bit64::is.integer64), as.integer))
}

# Output directory
dir.create(here("Results", "Feasibility2"), recursive = TRUE, showWarnings = FALSE)
out <- here("Results", "Feasibility2")

# ----------------------------------------------------------------------------
# 1. DATABASE CONNECTION
# ----------------------------------------------------------------------------
# Connection details are read from config.yml (not committed to git).
# Copy config.yml.example to config.yml and fill in your site-specific values.
cfg <- config::get()
db_name <- cfg$database$db_name

con <- DBI::dbConnect(
  odbc::odbc(),
  Driver   = cfg$database$driver,
  Server   = cfg$database$server,
  Database = cfg$database$database,
  trusted_connection = "yes"
)

cdm <- cdmFromCon(
  con         = con,
  cdmSchema   = cfg$database$cdm_schema,
  writeSchema = cfg$database$write_schema,
  writePrefix = cfg$database$write_prefix
)

message(glue("Connected to CDM: {cdmName(cdm)} | Version: {cdmVersion(cdm)}"))

# --- Server-side type patches for LTHT  -------------------------------------------------
if ("observation_period" %in% names(cdm)) {
  cdm$observation_period <- cdm$observation_period |>
    mutate(
      observation_period_start_date = as.Date(observation_period_start_date),
      observation_period_end_date = as.Date(observation_period_end_date)
    )
}
if ("measurement" %in% names(cdm)) {
  cdm$measurement <- cdm$measurement |>
    mutate(
      measurement_date = as.Date(measurement_date),
      measurement_time = as.character(measurement_time),
      measurement_event_id = sql("TRY_CAST(measurement_event_id AS BIGINT)")
    )
}
if ("observation" %in% names(cdm)) {
  cdm$observation <- cdm$observation |>
    mutate(
      observation_event_id = sql("TRY_CAST(observation_event_id AS BIGINT)")
    )
}
message("Server-side type patches applied.")

# --- Database overview -------------------------------------------------------
snapshot <- summariseOmopSnapshot(cdm)
tbl_snapshot <- tableOmopSnapshot(snapshot)
print(tbl_snapshot)
tbl_snapshot |> gtsave(file.path(out, "00_db_snapshot.html"))
message("Database snapshot saved.")

# --- Pre-flight validation ---------------------------------------------------
message("\n--- Pre-flight checks ---")
preflight_ok <- TRUE

for (tbl_name in c("person", "observation_period", "measurement", "condition_occurrence", "concept")) {
  n_check <- tryCatch(
    cdm[[tbl_name]] |> head(1) |> collect() |> nrow(),
    error = function(e) { message(glue("  FAIL: cannot access cdm${tbl_name}: {e$message}")); -1L }
  )
  if (n_check < 0) preflight_ok <- FALSE
}
if (preflight_ok) message("  OK: Core CDM tables accessible")

for (cid in c(egfr_known, creat_known)) {
  found <- cdm$concept |> filter(concept_id == !!cid) |> collect() |> nrow()
  if (found == 0) {
    message(glue("  FAIL: concept_id {cid} not found. Update egfr_known/creat_known."))
    preflight_ok <- FALSE
  }
}
if (preflight_ok) message(glue("  OK: Hardcoded concept IDs ({egfr_known}, {creat_known}) exist"))

obs_n <- cdm$observation_period |> summarise(n = n()) |> collect() |> pull(n)
if (obs_n == 0) { message("  FAIL: observation_period empty"); preflight_ok <- FALSE
} else { message(glue("  OK: observation_period has {obs_n} rows")) }

death_accessible <- tryCatch({ cdm$death |> head(1) |> collect(); TRUE }, error = function(e) FALSE)
if (death_accessible) { message("  OK: death table accessible")
} else { message("  WARN: death table not accessible") }

# patch_int64 check
tryCatch({
  test_row <- cdm$person |> head(1) |> select(person_id) |> collect()
  if (any(sapply(test_row, bit64::is.integer64))) {
    patch_int64(test_row); message("  OK: integer64 detected, patch_int64() works")
  } else { message("  OK: No integer64 (patch_int64 still safe)") }
}, error = function(e) message(glue("  WARN: patch_int64 check: {e$message}")))

# Write schema check
write_ok <- tryCatch({
  test_df <- tibble(subject_id = 1L, cohort_start_date = as.Date("2020-01-01"))
  cdm <- insertTable(cdm, name = "preflight_test", table = test_df)
  DBI::dbRemoveTable(con, DBI::Id(schema = cfg$database$write_schema, table = paste0(cfg$database$write_prefix, "preflight_test")))
  TRUE
}, error = function(e) { message(glue("  FAIL: writeSchema: {e$message}")); FALSE })
if (write_ok) { message("  OK: writeSchema writable") } else { preflight_ok <- FALSE }

if (!preflight_ok) stop("Pre-flight checks FAILED.")
message("  All pre-flight checks passed.\n")

# ============================================================================
# 2. CODELIST GENERATION (or load from cache)
# ============================================================================

if (generate_codelists) {
  message("\n========== GENERATING CODELISTS ==========")

  # --- 2a. Fracture concepts ---
  fracture_ancestor_ids <- c(4138412, 4129394, 4300192, 1245444, 4133610)
  thoracic_fx_codes <- getCandidateCodes(
    cdm = cdm, keywords = c("thoracic vertebra fracture", "thoracic spine fracture"),
    domains = "Condition", includeDescendants = TRUE
  )
  fracture_descendants <- getDescendants(cdm, conceptId = fracture_ancestor_ids)
  fracture_concept_ids <- unique(c(fracture_descendants$concept_id, thoracic_fx_codes$concept_id))
  fracture_codelist <- newCodelist(list("frailty_fracture" = fracture_concept_ids))
  message(glue("Fracture codelist: {length(fracture_concept_ids)} concepts"))

  # --- 2b. eGFR ---
  egfr_codes_search <- tryCatch(
    getCandidateCodes(cdm = cdm, keywords = c("estimated glomerular filtration rate", "egfr"),
                      domains = "Measurement", includeDescendants = TRUE) |> pull(concept_id),
    error = function(e) integer(0)
  )
  egfr_codes <- unique(c(egfr_known, egfr_codes_search))
  egfr_codelist <- newCodelist(list("egfr" = egfr_codes))
  message(glue("eGFR codelist: {length(egfr_codes)} concepts"))

  # --- 2c. Creatinine (exclude creatine kinase) ---
  creatinine_codes_search <- tryCatch({
    candidates <- getCandidateCodes(cdm = cdm, keywords = c("creatinine"),
                                    domains = "Measurement", includeDescendants = TRUE)
    candidates |> filter(!str_detect(str_to_lower(concept_name), "kinase|creatine k")) |> pull(concept_id)
  }, error = function(e) integer(0))
  creatinine_codes <- unique(c(creat_known, creatinine_codes_search))
  creatinine_codelist <- newCodelist(list("creatinine" = creatinine_codes))
  message(glue("Creatinine codelist: {length(creatinine_codes)} concepts"))

  # --- 2d. CKD diagnosis ---
  ckd_diag_ids <- getDescendants(cdm, conceptId = 46271022) |> pull(concept_id)
  ckd_diag_codelist <- newCodelist(list("ckd_diagnosis" = ckd_diag_ids))
  message(glue("CKD diagnosis codelist: {length(ckd_diag_ids)} concepts"))

  # --- 2e. Comorbidities ---
  comorbidity_keywords <- list(
    "Hypertension"               = c("essential hypertension", "hypertensive disease"),
    "Diabetes"                   = c("diabetes mellitus"),
    "Ischemic_Heart_Disease"     = c("ischemic heart disease", "coronary artery disease"),
    "Heart_Failure"              = c("heart failure", "congestive heart failure"),
    "Peripheral_Vascular_Disease"= c("peripheral vascular disease"),
    "Cerebrovascular_Disease"    = c("cerebrovascular disease", "stroke"),
    "Cancer"                     = c("malignant neoplasm"),
    "Osteoporosis"               = c("osteoporosis"),
    "Arthritis"                  = c("rheumatoid arthritis", "inflammatory arthritis"),
    "Hyperthyroidism"            = c("hyperthyroidism", "thyrotoxicosis"),
    "Hypoparathyroidism"         = c("hypoparathyroidism"),
    "Hyperparathyroidism"        = c("hyperparathyroidism"),
    "Undernutrition"             = c("malnutrition", "undernutrition"),
    "Dementia"                   = c("dementia"),
    "Hypocalcaemia"              = c("hypocalcemia", "hypocalcaemia"),
    "Hypercalcaemia"             = c("hypercalcemia", "hypercalcaemia"),
    "Calciphylaxis"              = c("calciphylaxis"),
    "Systemic_Inflammatory"      = c("systemic lupus", "vasculitis", "sarcoidosis"),
    "Chronic_Infections"         = c("hiv", "hepatitis c", "chronic hepatitis b"),
    "Substance_Abuse"            = c("substance abuse", "drug dependence", "alcohol dependence")
  )

  comorbidity_codes <- purrr::imap(comorbidity_keywords, function(kw, nm) {
    ids <- tryCatch(
      getCandidateCodes(cdm = cdm, keywords = kw, domains = "Condition",
                        includeDescendants = TRUE) |> pull(concept_id) |> unique(),
      error = function(e) { message(glue("  Warning: '{nm}' failed: {e$message}")); integer(0) }
    )
    message(glue("  {nm}: {length(ids)} concepts"))
    ids
  }) |> newCodelist()

  # --- 2f. Frailty Score (known concept + generated) ---
  # Uses LTHT-specific concept 40483383 plus keyword search for frailty/rockwood
  frailty_codes_search <- tryCatch(
    getCandidateCodes(cdm = cdm, keywords = c("frailty", "rockwood frailty score"),
                      domains = "Measurement", includeDescendants = TRUE) |> pull(concept_id),
    error = function(e) integer(0)
  )
  frailty_codes <- unique(c(frailty_known, frailty_codes_search))
  frailty_codelist <- newCodelist(list("frailty_score" = frailty_codes))
  message(glue("Frailty score codelist: {length(frailty_codes)} concepts (includes known {frailty_known})"))

  # --- 2g. IMD ---
  imd_source_concept_id <- 35812882

  # --- 2h. BMI ---
  bmi_codes <- tryCatch(
    getCandidateCodes(cdm = cdm, keywords = c("body mass index", "bmi"),
                      domains = "Measurement", includeDescendants = TRUE) |> pull(concept_id) |> unique(),
    error = function(e) integer(0)
  )
  height_codes <- tryCatch(
    getCandidateCodes(cdm = cdm, keywords = c("body height"),
                      domains = "Measurement", includeDescendants = TRUE) |> pull(concept_id) |> unique(),
    error = function(e) integer(0)
  )
  weight_codes <- tryCatch(
    getCandidateCodes(cdm = cdm, keywords = c("body weight"),
                      domains = "Measurement", includeDescendants = TRUE) |> pull(concept_id) |> unique(),
    error = function(e) integer(0)
  )
  message(glue("BMI: {length(bmi_codes)} concepts | Height: {length(height_codes)} | Weight: {length(weight_codes)}"))

  # --- 2i. Calcium (for lab-based hypercalcaemia definition) ---
  # Hypercalcaemia may not be coded as a condition in all datasets.
  # Lab-based definition: serum/plasma calcium > 2.6 mmol/L (upper limit of normal).
  # Excludes ionised calcium and urine calcium to avoid false positives.
  calcium_codes_search <- tryCatch({
    candidates <- getCandidateCodes(cdm = cdm, keywords = c("calcium"),
                                    domains = "Measurement", includeDescendants = TRUE)
    candidates |>
      filter(!str_detect(str_to_lower(concept_name), "ionized|ionised|urine|24.hour|24h|ratio|clearance|score")) |>
      pull(concept_id) |> unique()
  }, error = function(e) { message(glue("  Warning: Calcium codelist failed: {e$message}")); integer(0) })
  calcium_codelist <- newCodelist(list("calcium" = calcium_codes_search))
  message(glue("Calcium (serum/plasma) codelist: {length(calcium_codes_search)} concepts"))

  # --- Export codelists ---
  codelist_export <- bind_rows(
    tibble(codelist = "frailty_fracture", concept_id = fracture_concept_ids),
    tibble(codelist = "egfr",            concept_id = egfr_codes),
    tibble(codelist = "creatinine",      concept_id = creatinine_codes),
    tibble(codelist = "ckd_diagnosis",   concept_id = ckd_diag_ids),
    tibble(codelist = "bmi",             concept_id = bmi_codes),
    tibble(codelist = "height",          concept_id = height_codes),
    tibble(codelist = "weight",          concept_id = weight_codes),
    tibble(codelist = "calcium",         concept_id = calcium_codes_search),
    tibble(codelist = "frailty_score",   concept_id = frailty_codes),
    purrr::imap_dfr(as.list(comorbidity_codes), ~ tibble(codelist = .y, concept_id = .x))
  )
  write.csv(codelist_export, file.path(out, "all_codelists_for_partners.csv"), row.names = FALSE)
  message(glue("Codelists exported: {nrow(codelist_export)} rows, {n_distinct(codelist_export$codelist)} lists"))

  message("========== CODELIST GENERATION COMPLETE ==========\n")

} else {
  # --- Load cached codelists ---
  message("Loading cached codelists from CSV...")
  codelist_export <- read.csv(file.path(out, "all_codelists_for_partners.csv"))

  fracture_concept_ids <- codelist_export |> filter(codelist == "frailty_fracture") |> pull(concept_id)
  egfr_codes           <- codelist_export |> filter(codelist == "egfr") |> pull(concept_id)
  creatinine_codes     <- codelist_export |> filter(codelist == "creatinine") |> pull(concept_id)
  ckd_diag_ids         <- codelist_export |> filter(codelist == "ckd_diagnosis") |> pull(concept_id)
  bmi_codes            <- codelist_export |> filter(codelist == "bmi") |> pull(concept_id)
  height_codes         <- codelist_export |> filter(codelist == "height") |> pull(concept_id)
  weight_codes         <- codelist_export |> filter(codelist == "weight") |> pull(concept_id)

  fracture_codelist    <- newCodelist(list("frailty_fracture" = fracture_concept_ids))
  egfr_codelist        <- newCodelist(list("egfr" = egfr_codes))
  creatinine_codelist  <- newCodelist(list("creatinine" = creatinine_codes))
  ckd_diag_codelist    <- newCodelist(list("ckd_diagnosis" = ckd_diag_ids))

  calcium_codes     <- codelist_export |> filter(codelist == "calcium") |> pull(concept_id)
  calcium_codelist  <- newCodelist(list("calcium" = calcium_codes))

  frailty_codes     <- codelist_export |> filter(codelist == "frailty_score") |> pull(concept_id)
  frailty_codelist  <- newCodelist(list("frailty_score" = frailty_codes))

  comorbidity_names <- setdiff(
    unique(codelist_export$codelist),
    c("frailty_fracture", "egfr", "creatinine", "ckd_diagnosis", "bmi", "height", "weight", "calcium", "frailty_score")
  )
  comorbidity_codes <- purrr::map(comorbidity_names, ~ {
    codelist_export |> filter(codelist == .x) |> pull(concept_id)
  }) |> setNames(comorbidity_names) |> newCodelist()

  imd_source_concept_id <- 35812882

  message(glue("Loaded {n_distinct(codelist_export$codelist)} codelists from cache."))
}

# ============================================================================
# 3. DATABASE-WIDE CODE USE DIAGNOSTICS
# ============================================================================
message("\n========== DATABASE-WIDE CODE USE ==========")

message("Checking fracture code use...")
fracture_code_use <- summariseCodeUse(fracture_codelist, cdm = cdm)
tableCodeUse(fracture_code_use, type = "flextable") |>
  flextable::save_as_docx(path = file.path(out, "CodeUse_Fractures.docx"))

message("Checking eGFR code use...")
egfr_code_use <- summariseCodeUse(egfr_codelist, cdm = cdm)
tableCodeUse(egfr_code_use, type = "flextable") |>
  flextable::save_as_docx(path = file.path(out, "CodeUse_eGFR.docx"))

message("Checking creatinine code use...")
creatinine_code_use <- summariseCodeUse(creatinine_codelist, cdm = cdm)
tableCodeUse(creatinine_code_use, type = "flextable") |>
  flextable::save_as_docx(path = file.path(out, "CodeUse_Creatinine.docx"))

message("Checking frailty score code use...")
frailty_code_use <- summariseCodeUse(frailty_codelist, cdm = cdm)
tableCodeUse(frailty_code_use, type = "flextable") |>
  flextable::save_as_docx(path = file.path(out, "CodeUse_Frailty.docx"))

message("========== CODE USE DIAGNOSTICS COMPLETE ==========\n")

# ============================================================================
# 4. COHORT CONSTRUCTION (Fractures-first)
# ============================================================================
message("\n========== COHORT CONSTRUCTION ==========")

# --- 4a. Fracture cohort (N1) -----------------------------------------------
cdm$fractures <- conceptCohort(
  cdm        = cdm,
  conceptSet = fracture_codelist,
  name       = "fractures",
  exit       = "event_start_date"
)

# Apply date range and age filter
cdm$fractures <- cdm$fractures |>
  requireInDateRange(dateRange = c(study_start, study_end)) |>
  requireAge(ageRange = c(min_age, 150), indexDate = "cohort_start_date")

n1_count <- cohortCount(cdm$fractures)
message(glue("N1 (Fractures, age>={min_age}, {study_start} to {study_end}): {n1_count$number_subjects} subjects, {n1_count$number_records} records"))

# --- 4b. Extend cohort_end_date ---------------------------------------------
message("4b. Extending cohort_end_date...")

obs_periods <- cdm$observation_period |>
  select(person_id, observation_period_end_date) |>
  collect() |> patch_int64()

death_dates <- tryCatch(
  cdm$death |> select(person_id, death_date) |> collect() |> patch_int64(),
  error = function(e) tibble(person_id = integer(0), death_date = as.Date(character(0)))
)

fracture_local <- cdm$fractures |> collect() |> patch_int64()

# Keep first fracture per patient
fracture_extended <- fracture_local |>
  group_by(subject_id, cohort_definition_id) |>
  slice_min(cohort_start_date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  left_join(obs_periods, by = c("subject_id" = "person_id")) |>
  left_join(death_dates, by = c("subject_id" = "person_id")) |>
  mutate(
    cohort_end_date = pmin(observation_period_end_date, death_date, study_end, na.rm = TRUE)
  ) |>
  select(-observation_period_end_date, -death_date)

cdm <- insertTable(cdm, name = "fractures", table = fracture_extended)
cdm$fractures <- newCohortTable(cdm$fractures)

n1_subjects <- n_distinct(fracture_extended$subject_id)
message(glue("N1 (after dedup to first fracture): {n1_subjects} unique subjects"))

# ============================================================================
# 5. ATTRITION WATERFALL & ANALYSIS COHORT
# ============================================================================
message("\n========== ATTRITION WATERFALL ==========")

# --- 5a. N2: Have eGFR within ±egfr_window_days of fracture? ----------------
# Pull all eGFR measurements for fracture patients
fracture_persons <- fracture_extended |> distinct(subject_id, cohort_start_date)

all_egfr_raw <- cdm$measurement |>
  filter(measurement_concept_id %in% !!egfr_codes) |>
  select(person_id, measurement_date, value_as_number = value_as_number) |>
  collect() |> patch_int64()

# Join to fracture cohort, filter to ±window
egfr_around_fx <- fracture_persons |>
  inner_join(all_egfr_raw, by = c("subject_id" = "person_id")) |>
  mutate(days_from_fx = as.numeric(measurement_date - cohort_start_date)) |>
  filter(abs(days_from_fx) <= egfr_window_days)

n2_subjects <- n_distinct(egfr_around_fx$subject_id)
message(glue("N2 (eGFR within +/-{egfr_window_days} days of fracture): {n2_subjects} subjects"))

# --- 5b. N3: eGFR < threshold (ANALYSIS COHORT) -----------------------------
# Take the eGFR closest to fracture date per patient
egfr_closest <- egfr_around_fx |>
  group_by(subject_id, cohort_start_date) |>
  slice_min(abs(days_from_fx), n = 1, with_ties = FALSE) |>
  ungroup()

egfr_below_threshold <- egfr_closest |>
  filter(!is.na(value_as_number) & value_as_number < egfr_threshold)

n3_subjects <- n_distinct(egfr_below_threshold$subject_id)
message(glue("N3 (eGFR < {egfr_threshold}): {n3_subjects} subjects (ANALYSIS COHORT)"))

# Also report N with eGFR < 30
egfr_below_30 <- egfr_closest |>
  filter(!is.na(value_as_number) & value_as_number < 30)
n5_subjects <- n_distinct(egfr_below_30$subject_id)
message(glue("N5 (eGFR < 30): {n5_subjects} subjects"))

# Build analysis cohort table (N3)
analysis_cohort_ids <- egfr_below_threshold |> distinct(subject_id)

analysis_extended <- fracture_extended |>
  inner_join(analysis_cohort_ids, by = "subject_id")

cdm <- insertTable(cdm, name = "analysis_cohort", table = analysis_extended)
cdm$analysis_cohort <- newCohortTable(cdm$analysis_cohort)

message(glue("Analysis cohort inserted: {nrow(analysis_extended)} records"))

# --- 5c. N4: Have ≥1 creatinine in 5yr pre-fracture? ------------------------
all_creat_raw <- cdm$measurement |>
  filter(measurement_concept_id %in% !!creatinine_codes) |>
  select(person_id, measurement_date, value_as_number) |>
  collect() |> patch_int64()

analysis_persons <- analysis_extended |> distinct(subject_id, cohort_start_date)

creat_pre_fx <- analysis_persons |>
  inner_join(all_creat_raw, by = c("subject_id" = "person_id")) |>
  filter(
    measurement_date >= (cohort_start_date - lookback_creatinine_days),
    measurement_date <= cohort_start_date
  )

# Count creatinine measurements per patient
creat_count_per_patient <- creat_pre_fx |>
  group_by(subject_id) |>
  summarise(n_creat = n(), median_value = median(value_as_number, na.rm = TRUE), .groups = "drop")

# N4a: ≥1 creatinine in lookback period
n4a_subjects <- n_distinct(creat_pre_fx$subject_id)
n4a_pct <- round(100 * n4a_subjects / n3_subjects, 1)

# N4b: ≥5 creatinine in lookback period
n4b_subjects <- creat_count_per_patient |> filter(n_creat >= 5) |> nrow()
n4b_pct <- round(100 * n4b_subjects / n3_subjects, 1)

creat_summary <- creat_count_per_patient |>
  summarise(
    n_patients_with_creat = n(),
    median_count = median(n_creat),
    iqr_count = paste0(quantile(n_creat, 0.25), " - ", quantile(n_creat, 0.75)),
    median_value = round(median(median_value, na.rm = TRUE), 1),
    iqr_value = paste0(round(quantile(median_value, 0.25, na.rm = TRUE), 1), " - ",
                       round(quantile(median_value, 0.75, na.rm = TRUE), 1))
  )

message(glue("N4a (>=1 creatinine in {lookback_creatinine_days/365}yr pre-fracture): {n4a_subjects} ({n4a_pct}%) of analysis cohort"))
message(glue("N4b (>=5 creatinine in {lookback_creatinine_days/365}yr pre-fracture): {n4b_subjects} ({n4b_pct}%) of analysis cohort"))
message("Creatinine summary (5yr pre-fracture):")
print(creat_summary)

# --- Attrition table ---------------------------------------------------------
attrition_table <- tibble(
  Step = c(
    glue("N1: Frailty fracture (age>={min_age}, {study_start} to {study_end})"),
    glue("N2: eGFR measurement within +/-{egfr_window_days} days"),
    glue("N3: eGFR < {egfr_threshold} (analysis cohort)"),
    glue("N4a: >=1 creatinine in {lookback_creatinine_days/365}yr pre-fracture"),
    glue("N4b: >=5 creatinine in {lookback_creatinine_days/365}yr pre-fracture"),
    "N5: eGFR < 30"
  ),
  N = c(n1_subjects, n2_subjects, n3_subjects, n4a_subjects, n4b_subjects, n5_subjects),
  Pct_of_N1 = round(100 * c(n1_subjects, n2_subjects, n3_subjects, n4a_subjects, n4b_subjects, n5_subjects) / n1_subjects, 1),
  Pct_of_previous = c(
    100,
    round(100 * n2_subjects / n1_subjects, 1),
    round(100 * n3_subjects / n2_subjects, 1),
    round(100 * n4a_subjects / n3_subjects, 1),
    round(100 * n4b_subjects / n3_subjects, 1),
    round(100 * n5_subjects / n3_subjects, 1)
  )
)

message("\n--- Attrition Waterfall ---")
print(attrition_table)

attrition_table |>
  flextable::flextable() |> flextable::autofit() |>
  flextable::save_as_docx(path = file.path(out, "Attrition_Waterfall.docx"))

message("========== COHORT CONSTRUCTION COMPLETE ==========\n")

# ============================================================================
# 6. DEMOGRAPHIC MISSINGNESS (of N3 analysis cohort)
# ============================================================================
message("\n========== DEMOGRAPHIC MISSINGNESS ==========")

# --- Sex ---
cohort_with_sex <- cdm$analysis_cohort |> addSex() |> collect() |> patch_int64()
n_sex_available <- sum(!is.na(cohort_with_sex$sex) & cohort_with_sex$sex != "None")

# --- Ethnicity AND Race (check BOTH fields for valid codes) ---
person_eth_race <- cdm$person |>
  filter(person_id %in% !!analysis_persons$subject_id) |>
  select(person_id, ethnicity_concept_id, race_concept_id) |>
  collect() |> patch_int64()

n_ethnicity_available <- sum(!is.na(person_eth_race$ethnicity_concept_id) &
                               person_eth_race$ethnicity_concept_id != 0)
n_race_available <- sum(!is.na(person_eth_race$race_concept_id) &
                          person_eth_race$race_concept_id != 0)

# Report which field actually has data (UK OMOP often stores ethnicity in race_concept_id)
message(glue("  Ethnicity field (ethnicity_concept_id): {n_ethnicity_available} with valid code"))
message(glue("  Race field (race_concept_id): {n_race_available} with valid code"))

# Use whichever field has more data
n_eth_best <- max(n_ethnicity_available, n_race_available)
eth_source_field <- if (n_race_available > n_ethnicity_available) "race_concept_id" else "ethnicity_concept_id"
message(glue("  Using '{eth_source_field}' as ethnicity source ({n_eth_best} available)"))

# Detailed breakdown of ethnicity/race concept IDs used
eth_breakdown <- if (eth_source_field == "race_concept_id") {
  person_eth_race |> filter(race_concept_id != 0, !is.na(race_concept_id)) |>
    count(race_concept_id, name = "n") |> arrange(desc(n))
} else {
  person_eth_race |> filter(ethnicity_concept_id != 0, !is.na(ethnicity_concept_id)) |>
    count(ethnicity_concept_id, name = "n") |> arrange(desc(n))
}
message("  Ethnicity/Race concept ID distribution:")
print(eth_breakdown)

# --- IMD ---
imd_obs <- cdm$observation |>
  filter(observation_source_concept_id == imd_source_concept_id) |>
  select(person_id, observation_date, imd_quintile = value_as_number) |>
  collect() |> patch_int64()

# IMD: any date within study period (IMD is linked by LSOA, date is not critical)
imd_in_study <- analysis_persons |>
  left_join(imd_obs, by = c("subject_id" = "person_id")) |>
  filter(
    is.na(observation_date) |
    (observation_date >= study_start & observation_date <= study_end)
  ) |>
  group_by(subject_id) |>
  # Take most recent IMD within study period
  slice_max(observation_date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(has_imd = !is.na(imd_quintile))
n_imd_study_period <- sum(imd_in_study$has_imd)

# IMD: any ever (no date restriction)
imd_any_ever <- analysis_persons |>
  left_join(imd_obs, by = c("subject_id" = "person_id")) |>
  group_by(subject_id) |>
  summarise(has_imd = any(!is.na(imd_quintile)), .groups = "drop")
n_imd_any <- sum(imd_any_ever$has_imd)

message(glue("  IMD (any date in study period): {n_imd_study_period} | IMD (any ever): {n_imd_any}"))

# --- BMI ---
bmi_persons <- if (length(bmi_codes) > 0) {
  cdm$measurement |>
    filter(measurement_concept_id %in% !!bmi_codes,
           person_id %in% !!analysis_persons$subject_id) |>
    summarise(n = n_distinct(person_id)) |> collect() |> pull(n)
} else { 0L }

# --- Height ---
height_persons <- if (length(height_codes) > 0) {
  cdm$measurement |>
    filter(measurement_concept_id %in% !!height_codes,
           person_id %in% !!analysis_persons$subject_id) |>
    summarise(n = n_distinct(person_id)) |> collect() |> pull(n)
} else { 0L }

# --- Weight ---
weight_persons <- if (length(weight_codes) > 0) {
  cdm$measurement |>
    filter(measurement_concept_id %in% !!weight_codes,
           person_id %in% !!analysis_persons$subject_id) |>
    summarise(n = n_distinct(person_id)) |> collect() |> pull(n)
} else { 0L }

# --- Rural/Urban classification (via location table LSOA + ONS lookup) --------
# Load ONS Rural-Urban Classification lookup (included in project folder)
ru_lookup <- read.csv(here("rural_urban.csv"), stringsAsFactors = FALSE) |>
  select(LSOA21CD, rural_urban_flag = Rural.Urban.flag)

rural_urban_data <- tryCatch({
  # Join person -> location to get LSOA
  person_location <- cdm$person |>
    filter(person_id %in% !!analysis_persons$subject_id) |>
    select(person_id, location_id) |>
    inner_join(
      cdm$location |> select(location_id, location_source_value),
      by = "location_id"
    ) |>
    collect() |> patch_int64()

  # Join LSOA to ONS rural/urban lookup
  person_location |>
    mutate(has_lsoa = !is.na(location_source_value) & location_source_value != "") |>
    left_join(ru_lookup, by = c("location_source_value" = "LSOA21CD"))
}, error = function(e) {
  message(glue("  WARN: Could not access location table: {e$message}"))
  tibble(person_id = integer(0), has_lsoa = logical(0), rural_urban_flag = character(0))
})

n_rural_urban <- sum(!is.na(rural_urban_data$rural_urban_flag))
message(glue("Rural/Urban classified (via LSOA -> ONS lookup): {n_rural_urban} of {n3_subjects} persons"))

# Distribution
ru_distribution <- rural_urban_data |>
  filter(!is.na(rural_urban_flag)) |>
  count(rural_urban_flag) |>
  mutate(pct = round(100 * n / sum(n), 1))
message("Rural/Urban distribution:")
print(ru_distribution)

# --- Frailty Score (query measurement table for value_as_number) ---
frailty_codes_vec <- as.list(frailty_codelist)[["frailty_score"]]
all_frailty_raw <- cdm$measurement |>
  filter(measurement_concept_id %in% !!frailty_codes_vec,
         person_id %in% !!analysis_persons$subject_id) |>
  select(person_id, measurement_date, value_as_number) |>
  collect() |> patch_int64()

# Join to analysis cohort: latest frailty score up to fracture date
frailty_pre_fx <- analysis_persons |>
  inner_join(all_frailty_raw, by = c("subject_id" = "person_id")) |>
  filter(measurement_date <= cohort_start_date, !is.na(value_as_number))

# Latest frailty score per patient
frailty_latest <- frailty_pre_fx |>
  group_by(subject_id) |>
  slice_max(measurement_date, n = 1, with_ties = FALSE) |>
  ungroup()

n_frailty_available <- n_distinct(frailty_latest$subject_id)
message(glue("Frailty score available (any prior to fracture): {n_frailty_available} of {n3_subjects} ({round(100*n_frailty_available/n3_subjects,1)}%)"))

# Frailty score distribution (counts by score value)
# Rockwood Clinical Frailty Scale: 1 = Very Fit, 2 = Well, 3 = Managing Well,
# 4 = Vulnerable, 5 = Mildly Frail, 6 = Moderately Frail, 7 = Severely Frail,
# 8 = Very Severely Frail, 9 = Terminally Ill.  Lower = less frail.
frailty_score_dist <- frailty_latest |>
  count(value_as_number) |>
  mutate(
    pct = round(100 * n / sum(n), 1),
    # Suppress counts below threshold for disclosure control
    n_display = suppress_count(n),
    pct_display = suppress_pct(n, pct)
  ) |>
  arrange(value_as_number)
message("Frailty score distribution (Rockwood CFS: 1=Very Fit … 9=Terminally Ill):")
print(frailty_score_dist |> select(value_as_number, n_display, pct_display))

# --- Build missingness table ---
missingness_table <- tibble(
  Variable = c("Age", "Sex", "Ethnicity", "IMD (study period)", "IMD (any ever)",
               "Rural/Urban", "BMI", "Height", "Weight", "Frailty Score"),
  N_Total = n3_subjects,
  N_Available = c(
    n3_subjects,  # Age always derivable
    n_sex_available,
    n_eth_best,
    n_imd_study_period,
    n_imd_any,
    n_rural_urban,
    bmi_persons,
    height_persons,
    weight_persons,
    n_frailty_available
  )
) |>
  mutate(
    N_Missing = N_Total - N_Available,
    Pct_Available = round(100 * N_Available / N_Total, 1),
    Pct_Missing = round(100 * N_Missing / N_Total, 1)
  )

message("\n--- Demographic Missingness (Analysis Cohort, N3) ---")
print(missingness_table)

missingness_table |>
  flextable::flextable() |> flextable::autofit() |>
  flextable::save_as_docx(path = file.path(out, "Missingness_Analysis_Cohort.docx"))

message("========== MISSINGNESS COMPLETE ==========\n")

# ============================================================================
# 6b. COMPREHENSIVE TABLE 1 (Clinical & Person Characteristics)
# ============================================================================
message("\n========== TABLE 1: FULL CHARACTERISTICS ==========")

# --- Demographics via PatientProfiles ----------------------------------------
cdm$analysis_cohort <- cdm$analysis_cohort |>
  addDemographics(ageGroup = list(c(50, 64), c(65, 74), c(75, 84), c(85, 150))) |>
  addSex()

# --- Death ---
death_in_cohort <- analysis_persons |>
  left_join(death_dates, by = c("subject_id" = "person_id")) |>
  mutate(died = !is.na(death_date))
n_died <- sum(death_in_cohort$died)
n_died_pct <- round(100 * n_died / n3_subjects, 1)
message(glue("Deaths in analysis cohort: {n_died} ({n_died_pct}%)"))

# --- Comorbidity flags (prior to fracture) ---
cdm$analysis_cohort <- cdm$analysis_cohort |>
  addConceptIntersectFlag(
    conceptSet = comorbidity_codes,
    window     = c(-Inf, 0),
    nameStyle  = "{concept_name}_flag"
  )

# --- CKD diagnosis flag ---
cdm$analysis_cohort <- cdm$analysis_cohort |>
  addConceptIntersectFlag(
    conceptSet = ckd_diag_codelist,
    window     = c(-Inf, 0),
    nameStyle  = "ckd_dx_flag"
  )

# --- Lab-based Hypercalcaemia & Hypocalcaemia (unit-aware) -------------------
# Thresholds (adjusted calcium, per NICE CKD guidelines):
#   Hypercalcaemia: > 2.6 mmol/L  (or > 10.4 mg/dL)
#   Hypocalcaemia:  < 2.2 mmol/L  (or < 8.8 mg/dL)
# The code checks unit_concept_id and standardises to mmol/L before thresholding.
# OMOP standard unit concept IDs: 8753 = mmol/L, 8840 = mg/dL
# Heuristic fallback: if unit unknown, values > 20 assumed mg/dL (normal Ca ~2.2-2.6 mmol/L).
calcium_codes_vec <- as.list(calcium_codelist)[["calcium"]]
hypercalcaemia_lab <- if (length(calcium_codes_vec) > 0) {
  message("Computing lab-based calcium abnormalities (unit-aware)...")
  all_calcium_raw <- cdm$measurement |>
    filter(measurement_concept_id %in% !!calcium_codes_vec,
           person_id %in% !!analysis_persons$subject_id) |>
    select(person_id, measurement_date, value_as_number, unit_concept_id) |>
    collect() |> patch_int64()

  # Report unit distribution
  unit_dist <- all_calcium_raw |>
    count(unit_concept_id) |> arrange(desc(n))
  message("  Calcium unit_concept_id distribution:")
  print(unit_dist)

  # Standardise to mmol/L
  # unit_concept_id: 8753 = mmol/L, 8840 = mg/dL
  all_calcium_std <- all_calcium_raw |>
    mutate(
      ca_mmol = case_when(
        unit_concept_id == 8840 ~ value_as_number / 4,    # mg/dL → mmol/L (approx)
        unit_concept_id == 8753 ~ value_as_number,         # already mmol/L
        value_as_number > 20    ~ value_as_number / 4,     # heuristic: likely mg/dL
        TRUE                    ~ value_as_number           # assume mmol/L
      )
    )

  # Join to analysis cohort, restrict to prior-to-fracture measurements
  calcium_pre_fx <- analysis_persons |>
    inner_join(all_calcium_std, by = c("subject_id" = "person_id")) |>
    filter(measurement_date <= cohort_start_date,
           !is.na(ca_mmol))

  # Flag: any calcium > 2.6 mmol/L prior to fracture (hypercalcaemia)
  hypercalcaemia_patients <- calcium_pre_fx |>
    filter(ca_mmol > 2.6) |>
    distinct(subject_id)

  n_hypercalcaemia_lab <- nrow(hypercalcaemia_patients)
  pct_hypercalcaemia_lab <- round(100 * n_hypercalcaemia_lab / n3_subjects, 1)
  message(glue("  Hypercalcaemia (lab: Ca > 2.6 mmol/L): {n_hypercalcaemia_lab} ({pct_hypercalcaemia_lab}%) of analysis cohort"))
  message(glue("  Total calcium measurements (pre-fracture): {nrow(calcium_pre_fx)} across {n_distinct(calcium_pre_fx$subject_id)} patients"))

  list(n = n_hypercalcaemia_lab, pct = pct_hypercalcaemia_lab,
       n_with_ca = n_distinct(calcium_pre_fx$subject_id),
       n_measurements = nrow(calcium_pre_fx))
} else {
  message("  No calcium concept IDs available — skipping lab-based calcaemia.")
  list(n = 0L, pct = 0, n_with_ca = 0L, n_measurements = 0L)
}

# --- Lab-based Hypocalcaemia (Ca < 2.2 mmol/L / < 8.8 mg/dL) ---
# Uses the same unit-standardised calcium_pre_fx from above.
hypocalcaemia_lab <- if (length(calcium_codes_vec) > 0 && exists("calcium_pre_fx")) {
  message("Computing lab-based hypocalcaemia (Ca < 2.2 mmol/L)...")

  # Flag: any calcium < 2.2 mmol/L prior to fracture
  hypocalcaemia_patients <- calcium_pre_fx |>
    filter(ca_mmol < 2.2) |>
    distinct(subject_id)

  n_hypocalcaemia_lab <- nrow(hypocalcaemia_patients)
  pct_hypocalcaemia_lab <- round(100 * n_hypocalcaemia_lab / n3_subjects, 1)
  message(glue("  Hypocalcaemia (lab: Ca < 2.2 mmol/L): {n_hypocalcaemia_lab} ({pct_hypocalcaemia_lab}%) of analysis cohort"))

  list(n = n_hypocalcaemia_lab, pct = pct_hypocalcaemia_lab)
} else {
  message("  No calcium data available — skipping lab-based hypocalcaemia.")
  list(n = 0L, pct = 0)
}

# --- eGFR timeliness (days from fracture to closest eGFR) ---
# Reference: Levey et al. 2009 CKD-EPI; reporting timing per KDIGO guidance
egfr_timeliness <- egfr_below_threshold |>
  mutate(
    abs_days = abs(days_from_fx),
    timing_cat = case_when(
      abs_days <= 1  ~ "Same day",
      abs_days <= 7  ~ "Within 1 week",
      TRUE           ~ "1 week+"
    )
  )

egfr_timing_table <- egfr_timeliness |>
  count(timing_cat) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  arrange(match(timing_cat, c("Same day", "Within 1 week", "1 week+")))

message("\n--- eGFR Timeliness (days from fracture to closest eGFR) ---")
print(egfr_timing_table)

# --- IMD distribution (quintiles, from study period) ---
imd_distribution <- imd_in_study |>
  filter(has_imd) |>
  count(imd_quintile) |>
  mutate(pct = round(100 * n / sum(n), 1))

# --- Collect full Table 1 data -----------------------------------------------
table1_local <- cdm$analysis_cohort |> collect() |> patch_int64()

# Build Table 1
table1_rows <- list()

# Demographics
table1_rows$n_total <- tibble(Variable = "Total N", Value = as.character(n3_subjects), Pct = "100")

# Age
age_stats <- table1_local |>
  summarise(median = median(age, na.rm = TRUE),
            q25 = quantile(age, 0.25, na.rm = TRUE),
            q75 = quantile(age, 0.75, na.rm = TRUE))
table1_rows$age <- tibble(Variable = "Age, median (IQR)",
                           Value = glue("{age_stats$median} ({age_stats$q25}-{age_stats$q75})"), Pct = "")

# Age groups
age_groups <- table1_local |>
  count(age_group) |> mutate(pct = round(100 * n / sum(n), 1))
table1_rows$age_grp <- age_groups |>
  transmute(Variable = glue("  Age {age_group}"), Value = suppress_count(n), Pct = suppress_pct(n, pct))

# Sex
sex_dist <- table1_local |> count(sex) |> mutate(pct = round(100 * n / sum(n), 1))
table1_rows$sex <- sex_dist |>
  transmute(Variable = glue("  Sex: {sex}"), Value = suppress_count(n), Pct = suppress_pct(n, pct))

# Ethnicity
table1_rows$ethnicity <- tibble(
  Variable = glue("Ethnicity available (from {eth_source_field})"),
  Value = suppress_count(n_eth_best),
  Pct = suppress_pct(n_eth_best, round(100 * n_eth_best / n3_subjects, 1))
)

# Death
table1_rows$death <- tibble(Variable = "Death (any time after fracture)",
                             Value = suppress_count(n_died),
                             Pct = suppress_pct(n_died, n_died_pct))

# IMD
table1_rows$imd <- tibble(Variable = "IMD available (study period)",
                           Value = suppress_count(n_imd_study_period),
                           Pct = suppress_pct(n_imd_study_period, round(100 * n_imd_study_period / n3_subjects, 1)))

# Rural/Urban (with breakdown)
table1_rows$rural_urban <- bind_rows(
  tibble(Variable = "Rural/Urban classified (via LSOA)",
         Value = suppress_count(n_rural_urban),
         Pct = suppress_pct(n_rural_urban, round(100 * n_rural_urban / n3_subjects, 1))),
  ru_distribution |>
    transmute(Variable = glue("  {rural_urban_flag}"), Value = suppress_count(n), Pct = suppress_pct(n, pct))
)

# BMI
table1_rows$bmi <- tibble(Variable = "BMI available",
                           Value = suppress_count(bmi_persons),
                           Pct = suppress_pct(bmi_persons, round(100 * bmi_persons / n3_subjects, 1)))

# CKD diagnosis
ckd_flag_col <- names(table1_local)[str_detect(names(table1_local), "ckd_dx")]
n_ckd_dx <- if (length(ckd_flag_col) > 0) sum(table1_local[[ckd_flag_col[1]]] > 0, na.rm = TRUE) else 0L
table1_rows$ckd_dx <- tibble(Variable = "CKD diagnosis code (any prior)",
                              Value = suppress_count(n_ckd_dx),
                              Pct = suppress_pct(n_ckd_dx, round(100 * n_ckd_dx / n3_subjects, 1)))

# eGFR value
egfr_val_stats <- egfr_below_threshold |>
  summarise(median = round(median(value_as_number, na.rm = TRUE), 1),
            q25 = round(quantile(value_as_number, 0.25, na.rm = TRUE), 1),
            q75 = round(quantile(value_as_number, 0.75, na.rm = TRUE), 1))
table1_rows$egfr_val <- tibble(Variable = "Index eGFR, median (IQR)",
                                Value = glue("{egfr_val_stats$median} ({egfr_val_stats$q25}-{egfr_val_stats$q75})"),
                                Pct = "")

# eGFR < 30 subset
table1_rows$egfr30 <- tibble(Variable = "  eGFR < 30",
                              Value = suppress_count(n5_subjects),
                              Pct = suppress_pct(n5_subjects, round(100 * n5_subjects / n3_subjects, 1)))

# Creatinine
table1_rows$creat_1 <- tibble(
  Variable = glue("Creatinine >=1 ({lookback_creatinine_days/365}yr pre-fracture)"),
  Value = suppress_count(n4a_subjects), Pct = suppress_pct(n4a_subjects, n4a_pct)
)
table1_rows$creat_5 <- tibble(
  Variable = glue("Creatinine >=5 ({lookback_creatinine_days/365}yr pre-fracture)"),
  Value = suppress_count(n4b_subjects), Pct = suppress_pct(n4b_subjects, n4b_pct)
)
table1_rows$creat_count <- tibble(
  Variable = "  Median creatinine count (IQR)",
  Value = glue("{creat_summary$median_count} ({creat_summary$iqr_count})"), Pct = ""
)
table1_rows$creat_val <- tibble(
  Variable = "  Median creatinine value (IQR)",
  Value = glue("{creat_summary$median_value} ({creat_summary$iqr_value})"), Pct = ""
)

# Comorbidities
comorbidity_flag_cols <- names(table1_local)[str_detect(names(table1_local), "_flag$")]
comorbidity_flag_cols <- setdiff(comorbidity_flag_cols, ckd_flag_col)

comorbidity_table1 <- purrr::map_dfr(comorbidity_flag_cols, function(col) {
  n_pos <- sum(table1_local[[col]] > 0, na.rm = TRUE)
  pct_pos <- round(100 * n_pos / n3_subjects, 1)
  clean_name <- str_replace_all(col, "_flag$", "") |> str_replace_all("_", " ") |> str_to_title()
  tibble(Variable = glue("  {clean_name}"), Value = suppress_count(n_pos), Pct = suppress_pct(n_pos, pct_pos))
})

# Combine all
table1_full <- bind_rows(
  table1_rows$n_total,
  tibble(Variable = "--- Demographics ---", Value = "", Pct = ""),
  table1_rows$age,
  table1_rows$age_grp,
  table1_rows$sex,
  table1_rows$ethnicity,
  table1_rows$death,
  table1_rows$imd,
  table1_rows$rural_urban,
  table1_rows$bmi,
  tibble(Variable = "--- Kidney Function ---", Value = "", Pct = ""),
  table1_rows$egfr_val,
  table1_rows$egfr30,
  table1_rows$ckd_dx,
  table1_rows$creat_1,
  table1_rows$creat_5,
  table1_rows$creat_count,
  table1_rows$creat_val,
  tibble(Variable = "--- eGFR Timeliness ---", Value = "", Pct = ""),
  egfr_timing_table |> transmute(Variable = glue("  {timing_cat}"), Value = suppress_count(n), Pct = suppress_pct(n, pct)),
  tibble(Variable = "--- Comorbidities (prior to fracture) ---", Value = "", Pct = ""),
  comorbidity_table1,
  tibble(Variable = "--- Lab-based Comorbidities (measurement threshold) ---", Value = "", Pct = ""),
  tibble(Variable = "  Hypercalcaemia (lab: Ca > 2.6 mmol/L)",
         Value = suppress_count(hypercalcaemia_lab$n),
         Pct = suppress_pct(hypercalcaemia_lab$n, hypercalcaemia_lab$pct)),
  tibble(Variable = "  Hypocalcaemia (lab: Ca < 2.2 mmol/L)",
         Value = suppress_count(hypocalcaemia_lab$n),
         Pct = suppress_pct(hypocalcaemia_lab$n, hypocalcaemia_lab$pct)),
  tibble(Variable = "    Patients with any Ca measurement (pre-fracture)",
         Value = suppress_count(hypercalcaemia_lab$n_with_ca),
         Pct = suppress_pct(hypercalcaemia_lab$n_with_ca, round(100 * hypercalcaemia_lab$n_with_ca / n3_subjects, 1))),
  tibble(Variable = "--- Frailty Score (latest prior to fracture) ---", Value = "", Pct = ""),
  tibble(Variable = "  Frailty score available (Rockwood CFS: 1=Very Fit, 9=Terminally Ill)",
         Value = suppress_count(n_frailty_available),
         Pct = suppress_pct(n_frailty_available, round(100 * n_frailty_available / n3_subjects, 1))),
  frailty_score_dist |>
    transmute(Variable = glue("    Score {value_as_number}"), Value = n_display, Pct = pct_display)
)

message("\n--- TABLE 1 ---")
print(table1_full, n = Inf)

# Save Table 1
table1_full |>
  rename(`N / Value` = Value, `%` = Pct) |>
  flextable::flextable() |> flextable::autofit() |>
  flextable::save_as_docx(path = file.path(out, "Table1_Full_Characteristics.docx"))

message("Table 1 saved to Table1_Full_Characteristics.docx")
message("========== TABLE 1 COMPLETE ==========\n")

# ============================================================================
# 7. CONSOLIDATED FEASIBILITY REPORT (HTML)
# ============================================================================
message("\n========== GENERATING CONSOLIDATED REPORT ==========")

# --- Codelist summary for HTML (top concept IDs per list, max 15 per codelist) ---
codelist_html_rows <- codelist_export |>
  group_by(codelist) |>
  summarise(n_concepts = n(), concept_ids_sample = paste(head(concept_id, 15), collapse = ", "), .groups = "drop") |>
  mutate(concept_ids_sample = if_else(n_concepts > 15, paste0(concept_ids_sample, ", ..."), concept_ids_sample))

report_lines <- c(
  "<!DOCTYPE html><html><head><meta charset='UTF-8'>",
  "<title>CKD-Bone Feasibility Report</title>",
  "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:auto;padding:20px}",
  "table{border-collapse:collapse;width:100%;margin:15px 0}",
  "th,td{border:1px solid #ddd;padding:8px;text-align:left}",
  "th{background-color:#4472C4;color:white}",
  "tr:nth-child(even){background-color:#f2f2f2}",
  ".red{color:#c00;font-weight:bold}.green{color:#060;font-weight:bold}.amber{color:#c60;font-weight:bold}",
  "h1{color:#2F5496}h2{color:#4472C4;border-bottom:2px solid #4472C4;padding-bottom:5px}",
  "h3{color:#333;margin-top:20px}",
  ".note{background:#fff3cd;padding:10px;border-left:4px solid #ffc107;margin:10px 0}",
  ".key-stat{font-size:1.3em;font-weight:bold;color:#2F5496}",
  ".section-header td{background:#2F5496;color:white;font-weight:bold}",
  "code{background:#f4f4f4;padding:2px 5px;border-radius:3px;font-size:0.9em}",
  "</style></head><body>",
  "<h1>CKD-Bone Feasibility Report</h1>",
  glue("<p><strong>Generated:</strong> {Sys.time()} | <strong>Database:</strong> {db_name}</p>"),
  glue("<p><strong>Study period:</strong> {study_start} to {study_end} | <strong>Min age:</strong> {min_age}</p>"),
  glue("<p><strong>eGFR threshold:</strong> &lt; {egfr_threshold} | <strong>eGFR window:</strong> &plusmn;{egfr_window_days} days of fracture</p>"),
  "",
  # --- Section 1: Attrition ---
  "<h2>1. Attrition Waterfall</h2>",
  "<table><tr><th>Step</th><th>N</th><th>% of N1</th><th>% of Previous</th></tr>",
  purrr::pmap_chr(attrition_table, function(Step, N, Pct_of_N1, Pct_of_previous) {
    glue("<tr><td>{Step}</td><td>{N}</td><td>{Pct_of_N1}%</td><td>{Pct_of_previous}%</td></tr>")
  }),
  "</table>",
  "",
  # --- Section 2: Table 1 ---
  "<h2>2. Table 1: Baseline Characteristics (Analysis Cohort, N3)</h2>",
  "<table><tr><th>Variable</th><th>N / Value</th><th>%</th></tr>",
  purrr::pmap_chr(table1_full, function(Variable, Value, Pct) {
    if (str_detect(Variable, "^---")) {
      glue("<tr class='section-header'><td colspan='3'>{str_replace_all(Variable, '---', '')}</td></tr>")
    } else {
      glue("<tr><td>{Variable}</td><td>{Value}</td><td>{Pct}</td></tr>")
    }
  }),
  "</table>",
  "",
  # --- Section 3: Creatinine ---
  "<h2>3. Creatinine Characterisation</h2>",
  glue("<p>Patients with &ge;1 creatinine in {lookback_creatinine_days/365}yr pre-fracture: <span class='key-stat'>{n4a_subjects} ({n4a_pct}%)</span></p>"),
  glue("<p>Patients with &ge;5 creatinine in {lookback_creatinine_days/365}yr pre-fracture: <span class='key-stat'>{n4b_subjects} ({n4b_pct}%)</span></p>"),
  glue("<p>Median creatinine count per patient: <strong>{creat_summary$median_count}</strong> (IQR: {creat_summary$iqr_count})</p>"),
  glue("<p>Median creatinine value: <strong>{creat_summary$median_value}</strong> (IQR: {creat_summary$iqr_value})</p>"),
  "",
  # --- Section 4: Variable Availability ---
  "<h2>4. Variable Availability & Missingness</h2>",
  "<table><tr><th>Variable</th><th>N Total</th><th>N Available</th><th>N Missing</th><th>% Available</th><th>% Missing</th></tr>",
  purrr::pmap_chr(missingness_table, function(Variable, N_Total, N_Available, N_Missing, Pct_Available, Pct_Missing) {
    css_class <- case_when(Pct_Missing < 5 ~ "green", Pct_Missing < 40 ~ "amber", TRUE ~ "red")
    glue("<tr><td>{Variable}</td><td>{N_Total}</td><td>{N_Available}</td><td>{N_Missing}</td><td class='{css_class}'>{Pct_Available}%</td><td>{Pct_Missing}%</td></tr>")
  }),
  "</table>",
  glue("<p><em>Ethnicity source: <code>{eth_source_field}</code> ({n_eth_best} with valid code; {n_ethnicity_available} in ethnicity_concept_id, {n_race_available} in race_concept_id)</em></p>"),
  "",
  # --- Section 5: eGFR ---
  "<h2>5. eGFR Distribution &amp; Timeliness</h2>",
  glue("<p>Patients with eGFR measurement within &plusmn;{egfr_window_days} days: {n2_subjects} of {n1_subjects} fracture patients ({round(100*n2_subjects/n1_subjects,1)}%)</p>"),
  glue("<p>eGFR &lt; {egfr_threshold}: <span class='key-stat'>{n3_subjects}</span> ({round(100*n3_subjects/n2_subjects,1)}% of those measured)</p>"),
  glue("<p>eGFR &lt; 30: <span class='key-stat'>{n5_subjects}</span> ({round(100*n5_subjects/n2_subjects,1)}% of those measured)</p>"),
  "<h3>eGFR Timeliness (proximity to fracture date)</h3>",
  "<table><tr><th>Timing</th><th>N</th><th>%</th></tr>",
  purrr::pmap_chr(egfr_timing_table, function(timing_cat, n, pct) {
    glue("<tr><td>{timing_cat}</td><td>{n}</td><td>{pct}%</td></tr>")
  }),
  "</table>",
  "",
  # --- Section 6: Guidance ---
  "<h2>6. Missingness Guidance</h2>",
  "<div class='note'><ul>",
  "<li><strong>&lt;5% missing:</strong> Complete-case analysis acceptable</li>",
  "<li><strong>5-40% missing:</strong> Multiple imputation (MICE) recommended; include auxiliary variables</li>",
  "<li><strong>&gt;40% missing:</strong> Exclude from primary analysis; report descriptively only</li>",
  "<li><strong>Ethnicity:</strong> Likely &gt;90% missing in UK hospital EHR — exclude from primary models, report descriptively</li>",
  "<li><strong>Frailty:</strong> If unavailable at partner sites, consider Hospital Frailty Risk Score (HFRS) from ICD-10 codes (Gilbert et al. 2018)</li>",
  "<li><strong>IMD:</strong> Linked by LSOA. If partially available, multiple imputation assuming MAR is appropriate</li>",
  "<li><strong>BMI:</strong> If sparse, consider using weight alone or omitting from primary model</li>",
  "</ul></div>",
  "",
  # --- Section 7: Frailty Score ---
  "<h2>7. Frailty Score</h2>",
  "<p><em>Rockwood Clinical Frailty Scale: 1 = Very Fit, 2 = Well, 3 = Managing Well, 4 = Vulnerable, 5 = Mildly Frail, 6 = Moderately Frail, 7 = Severely Frail, 8 = Very Severely Frail, 9 = Terminally Ill. Lower score = less frail.</em></p>",
  glue("<p>Frailty score available for <span class='key-stat'>{n_frailty_available}</span> of {n3_subjects} patients ({round(100*n_frailty_available/n3_subjects,1)}%).</p>"),
  glue("<p>Concept IDs used: {paste(frailty_codes_vec, collapse=', ')} (includes LTHT-specific {frailty_known})</p>"),
  "<h3>Latest Frailty Score Distribution (prior to fracture)</h3>",
  glue("<p><em>Counts &lt;{suppression_threshold} suppressed for disclosure control.</em></p>"),
  "<table><tr><th>Score</th><th>N</th><th>%</th></tr>",
  purrr::pmap_chr(frailty_score_dist, function(value_as_number, n, pct, n_display, pct_display) {
    glue("<tr><td>{value_as_number}</td><td>{n_display}</td><td>{pct_display}</td></tr>")
  }),
  "</table>",
  "",
  # --- Section 8: Codelists ---
  "<h2>8. Codelist Summary</h2>",
  "<p>Full codelists exported to <code>all_codelists_for_partners.csv</code></p>",
  "<table><tr><th>Codelist</th><th>N Concepts</th><th>Sample Concept IDs</th></tr>",
  purrr::pmap_chr(codelist_html_rows, function(codelist, n_concepts, concept_ids_sample) {
    pretty_name <- str_replace_all(codelist, "_", " ") |> str_to_title()
    glue("<tr><td>{pretty_name}</td><td>{n_concepts}</td><td><code>{concept_ids_sample}</code></td></tr>")
  }),
  "</table>",
  "",
  # --- Footer ---
  "<hr>",
  glue("<p><em>Report generated by feasibility2.R | CKD-Bone study | {Sys.time()}</em></p>"),
  "</body></html>"
)

writeLines(report_lines, file.path(out, "Feasibility_Report.html"))
message(glue("Report saved: {file.path(out, 'Feasibility_Report.html')}"))

# Save key tables
missingness_table |>
  flextable::flextable() |> flextable::autofit() |>
  flextable::save_as_docx(path = file.path(out, "Missingness_Analysis_Cohort.docx"))

message("\n========== ALL OUTPUTS COMPLETE ==========")
message(glue("Results directory: {out}"))

# ============================================================================
# 8. DISCONNECT
# ============================================================================
CDMConnector::cdmDisconnect(cdm)
message("Feasibility 2 complete. Database disconnected.")

