# data/

This folder is not tracked in the repository (see `.gitignore`) because the
underlying surveillance data are not freely redistributable. To run the
pipeline, place the following three files here, using exactly these names:

| File | Description |
|---|---|
| `annual_data.xlsx` | Annual state-level Lassa fever case totals |
| `Cases_rainfal_data.xlsx` | Weekly state-level case counts and rainfall (mm), columns: `state`, `epi_week`, `year`, `cases`, `rainfall_mm` |
| `climate_features.csv` | Weekly NDVI / temperature / soil-moisture features (not currently consumed by any pipeline stage; reserved for future extensions) |
| `gadm41_NGA_1.shp` (+ companion `.dbf`/`.shx`/`.prj`) | Nigeria state-boundary shapefile (GADM v4.1, level 1), used for the choropleth maps in `02_spatial_choropleth_moran.R` |

`00_setup.R` defines the exact paths it expects (`DATA_DIR <- "data"`).

