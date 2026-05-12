# Data

Raw data files are not committed due to file size limits.

## Download instructions

### CAL FIRE DINS — primary dataset (~58MB)
Two separate files for each fire:
- Palisades and Eaton: [https://gis.data.cnra.ca.gov/datasets/CALFIRE-Forestry::dins-2025-palisades-public-view](https://gis.data.cnra.ca.gov/datasets/994d3dc4569640caadbbc3198d5a3da1_0/explore?location=34.543492%2C-117.853526%2C8)

On each page: Download → Spreadsheet (CSV). Place both in `data/raw/`.

### Census TIGER tract boundaries (~32MB)
- https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2023&layergroup=Census+Tracts
- Select California → Download → unzip into `data/raw/tl_2023_06_tract/`

### Census ACS income data
Pulled live from the Census API inside `src/01_xgboost_shap.ipynb`.
Get a free API key at https://api.census.gov/data/key_signup.html
and set it as an environment variable: `export CENSUS_API_KEY='your_key'`
