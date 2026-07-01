# =============================================================================
# CKD and Bone Disease Feasibility Study
# OMOP CDM Analysis
#
# Description:
#   Investigates the feasibility of a study examining relationships between
#   frailty fractures and chronic kidney disease (CKD) using routinely
#   collected data mapped to the OMOP Common Data Model (CDM).
#
# Prerequisites:
#   See README.md for full list of required R packages and setup instructions.
# =============================================================================

library(DatabaseConnector)
library(SqlRender)
library(dplyr)
library(readr)
library(here)

# =============================================================================
# Configuration
# =============================================================================

# Update connection details to match your database environment.
# Do NOT hard-code credentials here; use environment variables or a secure
# secrets manager instead.
connectionDetails <- createConnectionDetails(
  dbms     = Sys.getenv("OMOP_DBMS", "postgresql"),
  server   = Sys.getenv("OMOP_SERVER"),
  user     = Sys.getenv("OMOP_USER"),
  password = Sys.getenv("OMOP_PASSWORD"),
  port     = as.integer(Sys.getenv("OMOP_PORT", "5432"))
)

cdmDatabaseSchema    <- Sys.getenv("OMOP_CDM_SCHEMA",    "cdm")
resultsDatabaseSchema <- Sys.getenv("OMOP_RESULTS_SCHEMA", "results")
outputFolder         <- here("output")

# =============================================================================
# Helpers
# =============================================================================

#' Execute a parameterised SQL query and return results as a data frame.
#'
#' @param connection  A DatabaseConnector connection object.
#' @param sql         SQL string (SqlRender dialect).
#' @param ...         Named substitution parameters passed to SqlRender.
#' @return            A data frame of query results.
executeQuery <- function(connection, sql, ...) {
  renderedSql <- SqlRender::render(sql, ...)
  translatedSql <- SqlRender::translate(renderedSql,
                                        targetDialect = connectionDetails$dbms)
  DatabaseConnector::querySql(connection, translatedSql)
}

# =============================================================================
# Load Reference Tables
# =============================================================================

ckdConcepts      <- read_csv(here("reference_tables", "ckd_concepts.csv"),
                             show_col_types = FALSE)
fractureConcepts <- read_csv(here("reference_tables", "fracture_concepts.csv"),
                             show_col_types = FALSE)

ckdConceptIds      <- paste(ckdConcepts$concept_id,      collapse = ", ")
fractureConceptIds <- paste(fractureConcepts$concept_id, collapse = ", ")

# =============================================================================
# Database Queries
# =============================================================================

connection <- DatabaseConnector::connect(connectionDetails)
on.exit(DatabaseConnector::disconnect(connection), add = TRUE)

# --- 1. Overall CDM population size -------------------------------------------

sqlPopulation <- "
SELECT COUNT(DISTINCT person_id) AS n_patients
FROM @cdm_schema.person;
"

populationCount <- executeQuery(connection, sqlPopulation,
                                cdm_schema = cdmDatabaseSchema)

# --- 2. CKD cohort ------------------------------------------------------------

sqlCkd <- "
SELECT COUNT(DISTINCT co.person_id) AS n_ckd_patients
FROM @cdm_schema.condition_occurrence co
WHERE co.condition_concept_id IN (@ckd_concept_ids);
"

ckdCount <- executeQuery(connection, sqlCkd,
                         cdm_schema      = cdmDatabaseSchema,
                         ckd_concept_ids = ckdConceptIds)

# --- 3. Fracture cohort -------------------------------------------------------

sqlFracture <- "
SELECT COUNT(DISTINCT co.person_id) AS n_fracture_patients
FROM @cdm_schema.condition_occurrence co
WHERE co.condition_concept_id IN (@fracture_concept_ids);
"

fractureCount <- executeQuery(connection, sqlFracture,
                              cdm_schema          = cdmDatabaseSchema,
                              fracture_concept_ids = fractureConceptIds)

# --- 4. Patients with BOTH CKD and fracture -----------------------------------

sqlOverlap <- "
SELECT COUNT(DISTINCT ckd.person_id) AS n_overlap_patients
FROM (
  SELECT DISTINCT person_id
  FROM @cdm_schema.condition_occurrence
  WHERE condition_concept_id IN (@ckd_concept_ids)
) ckd
INNER JOIN (
  SELECT DISTINCT person_id
  FROM @cdm_schema.condition_occurrence
  WHERE condition_concept_id IN (@fracture_concept_ids)
) frac ON ckd.person_id = frac.person_id;
"

overlapCount <- executeQuery(connection, sqlOverlap,
                             cdm_schema          = cdmDatabaseSchema,
                             ckd_concept_ids     = ckdConceptIds,
                             fracture_concept_ids = fractureConceptIds)

# --- 5. Age and sex breakdown for the overlap cohort -------------------------

sqlDemographics <- "
SELECT
  p.gender_concept_id,
  FLOOR((YEAR(co_min.min_date) - p.year_of_birth) / 10) * 10 AS age_group_start,
  COUNT(DISTINCT p.person_id) AS n_patients
FROM @cdm_schema.person p
INNER JOIN (
  SELECT person_id, MIN(condition_start_date) AS min_date
  FROM @cdm_schema.condition_occurrence
  WHERE condition_concept_id IN (@ckd_concept_ids)
  GROUP BY person_id
) co_min ON p.person_id = co_min.person_id
WHERE p.person_id IN (
  SELECT DISTINCT person_id
  FROM @cdm_schema.condition_occurrence
  WHERE condition_concept_id IN (@fracture_concept_ids)
)
GROUP BY p.gender_concept_id,
         FLOOR((YEAR(co_min.min_date) - p.year_of_birth) / 10) * 10
ORDER BY age_group_start, p.gender_concept_id;
"

demographicsData <- executeQuery(connection, sqlDemographics,
                                 cdm_schema          = cdmDatabaseSchema,
                                 ckd_concept_ids     = ckdConceptIds,
                                 fracture_concept_ids = fractureConceptIds)

# --- 6. Observation period coverage ------------------------------------------

sqlObservation <- "
SELECT
  AVG(DATEDIFF(day, observation_period_start_date,
                    observation_period_end_date)) AS mean_obs_days,
  MIN(DATEDIFF(day, observation_period_start_date,
                    observation_period_end_date)) AS min_obs_days,
  MAX(DATEDIFF(day, observation_period_start_date,
                    observation_period_end_date)) AS max_obs_days
FROM @cdm_schema.observation_period op
WHERE op.person_id IN (
  SELECT DISTINCT person_id
  FROM @cdm_schema.condition_occurrence
  WHERE condition_concept_id IN (@ckd_concept_ids)
)
AND op.person_id IN (
  SELECT DISTINCT person_id
  FROM @cdm_schema.condition_occurrence
  WHERE condition_concept_id IN (@fracture_concept_ids)
);
"

observationData <- executeQuery(connection, sqlObservation,
                                cdm_schema          = cdmDatabaseSchema,
                                ckd_concept_ids     = ckdConceptIds,
                                fracture_concept_ids = fractureConceptIds)

# =============================================================================
# Compile Feasibility Summary
# =============================================================================

feasibilitySummary <- data.frame(
  metric = c(
    "Total patients in CDM",
    "Patients with CKD",
    "Patients with fracture",
    "Patients with both CKD and fracture",
    "Mean observation days (overlap cohort)",
    "Min observation days (overlap cohort)",
    "Max observation days (overlap cohort)"
  ),
  value = c(
    populationCount$N_PATIENTS,
    ckdCount$N_CKD_PATIENTS,
    fractureCount$N_FRACTURE_PATIENTS,
    overlapCount$N_OVERLAP_PATIENTS,
    round(observationData$MEAN_OBS_DAYS, 1),
    observationData$MIN_OBS_DAYS,
    observationData$MAX_OBS_DAYS
  )
)

# =============================================================================
# Write Outputs
# =============================================================================

if (!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

write_csv(feasibilitySummary,
          file.path(outputFolder, "feasibility_summary.csv"))

write_csv(demographicsData,
          file.path(outputFolder, "demographics_by_age_sex.csv"))

message("Feasibility study complete. Results written to: ", outputFolder)
print(feasibilitySummary)
