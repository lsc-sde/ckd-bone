# CKD-Bone: Frailty Fractures in Patients with CKD

Feasibility analysis for studying frailty fractures in patients with chronic kidney disease (CKD), using OMOP CDM and OHDSI R packages.

## Study Design

**Base cohort:** Patients aged ≥50 with a frailty fracture (hip, thoracic/lumbar vertebra, pelvis, wrist, proximal forearm) in the study period.

**Attrition waterfall:**
1. N1: Frailty fracture (age ≥50, study period)
2. N2: eGFR measurement within ±4 weeks of fracture
3. N3: eGFR < 44 → **Analysis cohort**
4. N4: ≥1 creatinine in 5-year pre-fracture lookback
5. N5: eGFR < 30

## Requirements

- R ≥ 4.3
- Access to an OMOP CDM (SQL Server) with write permissions
- `renv` for package management

## Setup

```bash
# 1. Clone the repo and restore packages
renv::restore()

# 2. Copy the config template and fill in your site-specific values
cp config.yml.example config.yml
```

Edit `config.yml` with your database connection details:

```yaml
default:
  database:
    driver: "ODBC Driver 17 for SQL Server"
    server: "YOUR_SERVER_HERE"
    database: "YOUR_DATABASE_HERE"
    cdm_schema: "YOUR_CDM_SCHEMA"
    write_schema: "YOUR_WRITE_SCHEMA"
    write_prefix: "feas2_"
    db_name: "YOUR_DB_DISPLAY_NAME"
```

> **Note:** `config.yml` is in `.gitignore` and will not be committed. Each site maintains their own copy.

Then in R:

```r
# 3. (Optional) Adjust study parameters in feasibility_run.R section 0
# 4. Set generate_codelists = TRUE for first run
# 5. Source or run interactively
source("feasibility_run.R")
```

## Key Files

| File | Description |
|------|-------------|
| `feasibility_run.R` | Main analysis script (run interactively, section by section) |
| `CKD-Bone.Rproj` | RStudio project file (sets working directory, editor preferences) |
| `.Rprofile` | Bootstraps `renv` on session start |
| `rural_urban.csv` | ONS Rural-Urban Classification 2021 lookup (LSOA → Rural/Urban) |
| `config.yml` | Site-specific database configuration (not committed — copy from `config.yml.example`) |
| `config.yml.example` | Template config with placeholder values for partner sites |
| `renv/` | Package management (renv lockfile + activate script) |
| `Results/Feasibility2/` | Output directory (HTML report, DOCX tables, codelist CSV) |

## Configurable Parameters (Section 0)

```r
study_start            <- as.Date("2023-01-01")
study_end              <- as.Date("2026-06-25")
min_age                <- 50
egfr_threshold         <- 44       # eGFR cut-off
egfr_window_days       <- 28       # ±4 weeks around fracture
lookback_creatinine_days <- 1825   # 5 years
generate_codelists     <- TRUE     # Set FALSE after first run (uses cached CSV)
```

## Outputs

- `Feasibility_Report.html` — Self-contained HTML report with all tables and findings
- `Table1_Full_Characteristics.docx` — Full Table 1 (demographics, death, kidney function, comorbidities, eGFR timeliness)
- `Attrition_Waterfall.docx` — N1→N5 waterfall
- `Missingness_Analysis_Cohort.docx` — Variable availability/missingness
- `all_codelists_for_partners.csv` — Exportable codelists for partner OMOP sites
- `CodeUse_*.docx` — Database-wide code usage reports

## Important Notes on Feasibility Code

- **LTHT-specific patches:** The R script includes server-side type patches specific to the LTHT OMOP instance (integer64 handling, date casting via `TRY_CAST`, `measurement_time` coercion). Partner sites may need to adjust or remove these patches (see Section 1, "Server-side type patches for LTHT") depending on their CDM implementation.
- **IMD (Index of Multiple Deprivation):** The feasibility code uses IMD as a concept from the OMOP CDM. Sites must have IMD mapped as a concept in their CDM for deprivation analysis to work, or should be able to adjust the code to join IMD from LSOA (see location note below).
- **Rural/Urban Classification:** The code uses LSOA (Lower Layer Super Output Area) from the `location_source_value` field to map patients to rural/urban categories via the `rural_urban.csv` lookup file.
- **Frailty Score:** Uses LTHT-specific concept ID `40483383` plus keyword-generated concepts for "frailty" and "rockwood frailty score". Partner sites should verify which frailty concept IDs are available in their local vocabulary.

## Partner Site Instructions

1. Clone this repo
2. Run `renv::restore()` to install packages
3. Copy `config.yml.example` to `config.yml` and fill in your site-specific database details:
   - `server` — your SQL Server hostname
   - `database` — your OMOP database name
   - `cdm_schema` — schema containing OMOP CDM tables
   - `write_schema` — schema where you have write permissions (for temp cohort tables)
   - `write_prefix` — prefix for tables written by the script (default `"feas2_"`)
   - `db_name` — a display name for your database (used in reports)
4. Set `generate_codelists = TRUE` (first run) or use the shared `all_codelists_for_partners.csv`
5. Run interactively — pre-flight checks will validate your config and permissions immediately

## OMOP Packages Used

- [CDMConnector](https://darwin-eu.github.io/CDMConnector/)
- [CodelistGenerator](https://darwin-eu.github.io/CodelistGenerator/)
- [CohortConstructor](https://ohdsi.github.io/CohortConstructor/)
- [PatientProfiles](https://darwin-eu.github.io/PatientProfiles/)
- [CohortCharacteristics](https://darwin-eu.github.io/CohortCharacteristics/)
- [PhenotypeR](https://darwin-eu.github.io/PhenotypeR/)
- [OmopSketch](https://darwin-eu.github.io/OmopSketch/)
