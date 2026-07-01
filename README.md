# CKD-Bone: Frailty Fractures in Patients with CKD

Feasibility analysis for Round 2 grant application studying frailty fractures in patients with chronic kidney disease (CKD), using OMOP CDM and OHDSI R packages.

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

```r
# Restore packages via renv
renv::restore()

# Open feasibility_run.R and update:
# 1. Database connection parameters (Server, Database, schema)
# 2. Study parameters (dates, thresholds) in section 0
# 3. Set generate_codelists = TRUE for first run
```

## Key Files

| File | Description |
|------|-------------|
| `feasibility_run.R` | Main analysis script (run interactively, section by section) |
| `CKD-Bone.Rproj` | RStudio project file (sets working directory, editor preferences) |
| `.Rprofile` | Bootstraps `renv` on session start |
| `rural_urban.csv` | ONS Rural-Urban Classification 2021 lookup (LSOA → Rural/Urban) |
| `config.yml` | Project configuration |
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

- **IMD (Index of Multiple Deprivation):** The feasibility code uses IMD as a concept from the OMOP CDM. Sites must have IMD mapped as a concept in their CDM for deprivation analysis to work.
- **Rural/Urban Classification:** The code uses LSOA (Lower Layer Super Output Area) from the `location_source_value` field to map patients to rural/urban categories via the `rural_urban.csv` lookup file.

## Partner Site Instructions

1. Clone this repo
2. `renv::restore()` to install packages
3. Update connection parameters in `feasibility_run.R` section 1
4. Set `generate_codelists = TRUE` (first run) or use the shared `all_codelists_for_partners.csv`
5. Run interactively — pre-flight checks will validate your setup immediately

## OMOP Packages Used

- [CDMConnector](https://darwin-eu.github.io/CDMConnector/)
- [CodelistGenerator](https://darwin-eu.github.io/CodelistGenerator/)
- [CohortConstructor](https://ohdsi.github.io/CohortConstructor/)
- [PatientProfiles](https://darwin-eu.github.io/PatientProfiles/)
- [CohortCharacteristics](https://darwin-eu.github.io/CohortCharacteristics/)
- [PhenotypeR](https://darwin-eu.github.io/PhenotypeR/)
- [OmopSketch](https://darwin-eu.github.io/OmopSketch/)