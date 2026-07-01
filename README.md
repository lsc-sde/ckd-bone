# ckd-bone

Feasibility study investigating the relationships between frailty fractures and chronic kidney disease (CKD) using routinely collected health data mapped to the [OMOP Common Data Model (CDM)](https://ohdsi.github.io/CommonDataModel/).

---

## Repository structure

```
ckd-bone/
├── feasibility.R          # Main analysis script
├── reference_tables/      # OMOP concept lookup tables
│   ├── ckd_concepts.csv
│   └── fracture_concepts.csv
├── output/                # Generated results (excluded from version control)
└── README.md
```

---

## Prerequisites

### R version

R ≥ 4.2.0 is recommended.

### Required R packages

Install all dependencies from CRAN and GitHub before running the script:

```r
install.packages(c("dplyr", "readr", "here"))

# HADES packages (from GitHub via remotes)
install.packages("remotes")
remotes::install_github("OHDSI/SqlRender")
remotes::install_github("OHDSI/DatabaseConnector")
```

### Database driver

`DatabaseConnector` requires a JDBC driver for your database. Download the appropriate driver using:

```r
DatabaseConnector::downloadJdbcDrivers("postgresql")  # or "sql server", "redshift", etc.
```

Refer to the [DatabaseConnector documentation](https://ohdsi.github.io/DatabaseConnector/) for full driver setup instructions.

### Database access

You need read access to an OMOP CDM database. The following environment variables must be set before running the script:

| Variable | Description |
|---|---|
| `OMOP_DBMS` | Database type (`postgresql`, `sql server`, `redshift`, …) |
| `OMOP_SERVER` | Server hostname and database name |
| `OMOP_USER` | Database username |
| `OMOP_PASSWORD` | Database password |
| `OMOP_PORT` | Port (default: `5432`) |
| `OMOP_CDM_SCHEMA` | Schema containing the CDM tables (default: `cdm`) |
| `OMOP_RESULTS_SCHEMA` | Schema for any writeable results tables (default: `results`) |

Set these in your shell or in an `.Renviron` file (do **not** hard-code credentials in the script):

```bash
export OMOP_DBMS="postgresql"
export OMOP_SERVER="your-server/your-database"
export OMOP_USER="your-username"
export OMOP_PASSWORD="your-password"
export OMOP_CDM_SCHEMA="cdm"
export OMOP_RESULTS_SCHEMA="results"
```

---

## Running the feasibility study

1. Clone the repository and open a terminal in the project root.
2. Set the environment variables described above.
3. Run the script from the terminal:

```bash
Rscript feasibility.R
```

or from an interactive R session:

```r
source(here::here("feasibility.R"))
```

Results will be written to the `output/` directory:

| File | Description |
|---|---|
| `output/feasibility_summary.csv` | Top-level counts and observation period statistics |
| `output/demographics_by_age_sex.csv` | Age-band and sex breakdown for the overlap cohort |

---

## Reference tables

The `reference_tables/` directory contains the OMOP standard concept IDs used to define study cohorts.

| File | Contents |
|---|---|
| `ckd_concepts.csv` | SNOMED concepts for CKD stages 1–5 and ESRD |
| `fracture_concepts.csv` | SNOMED concepts for fragility / frailty fractures |

To extend or refine the concept sets, add or remove rows from these CSV files and re-run the script. Concept IDs can be browsed using [ATHENA](https://athena.ohdsi.org/).

---

## Notes

- All patient counts in the output are aggregated; no row-level patient data are written to disk.
- The script uses `SqlRender` to translate SQL into the dialect of your target database, so it should run against any OMOP-compliant data source.
