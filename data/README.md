# Data Directory

This directory contains the input datasets required to reproduce the analyses.

The original surveillance and climate datasets are **not included** in this repository because they are subject to data-sharing restrictions. Users should obtain the data from the appropriate sources and place them in this directory using the filenames expected by the pipeline.

## Required files

| File | Description |
|------|-------------|
| `annual_data.xlsx` | Annual state-level Lassa fever case totals. |
| `Cases_rainfal_data.xlsx` | Weekly state-level Lassa fever case counts and rainfall (mm). Expected columns include `state`, `epi_week`, `year`, `cases`, and `rainfall_mm`. |
| `gadm41_NGA_1.shp` | Nigeria state boundary shapefile (GADM v4.1, level 1). The accompanying `.dbf`, `.shx`, `.prj`, and other associated files must also be present. |

## Directory structure

```
data/
├── README.md
├── raw/
│   ├── annual_data.xlsx
│   ├── Cases_rainfal_data.xlsx
│   └── gadm41_NGA_1.*
└── processed/
```

The pipeline expects the data directory to be defined in `000_setup.R` as:

```r
DATA_DIR <- "data"
```

If your scripts instead reference `data/raw/`, update the `DATA_DIR` variable accordingly before running the pipeline.
