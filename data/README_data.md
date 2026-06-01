# Data Description

## File in This Folder

`pwt_empirical_panel.xlsx` — the pre-processed estimation panel used in
the empirical application. Four sheets:

| Sheet | Contents |
|-------|----------|
| `README` | Data description and variable definitions |
| `2.EstimationPanel` | Full 4,147-row estimation dataset (load this in R) |
| `3.CountryInfo` | One row per country with static attributes |
| `4.SummaryStats` | Descriptive statistics of all variables |

---

## Sample Dimensions

| | |
|---|---|
| Countries (N) | 143 |
| Estimation years (T) | 29 (1995 to 2023) |
| Observations | 4,147 (fully balanced) |
| Lag year (construction only) | 1994 (not in estimation panel) |

---

## Variable Dictionary

The `2.EstimationPanel` sheet contains the following columns:

| Column | Description | Unit | Source |
|--------|-------------|------|--------|
| `id` | Numeric panel ID (1–143), sorted alphabetically by ISO code | integer | Constructed |
| `countrycode` | 3-letter ISO 3166-1 alpha-3 code | string | PWT 11.0 |
| `country` | Country name | string | PWT 11.0 |
| `year` | Calendar year | integer | PWT 11.0 |
| `ly` | Log real output per capita: `log(rgdpo / pop)` | log(million 2021 USD / million persons) | PWT 11.0 |
| `lylag` | `ly` lagged one year | same as `ly` | Constructed |
| `lk` | Log capital per capita: `log(rnna / pop)` | log(million 2021 USD / million persons) | PWT 11.0 |
| `lh` | Log human capital index: `log(hc)` | log index | PWT 11.0 |
| `larea` | Log land area: `log(land_area_sqkm)` | log(sq km) — **time-invariant** | World Bank |
| `xbar_lk` | Within-country mean of `lk` over 1995–2023 | same as `lk` | Constructed |
| `xbar_lh` | Within-country mean of `lh` over 1995–2023 | same as `lh` | Constructed |
| `land_area_sqkm` | Land area in square kilometres | sq km — **time-invariant** | World Bank |

**PWT 11.0 source variables used:**

| PWT variable | Description |
|---|---|
| `rgdpo` | Output-side real GDP at chained PPPs (million 2021 USD) |
| `pop` | Population (millions) |
| `rnna` | Capital stock at constant national prices (million 2021 USD) |
| `hc` | Human capital index (based on Barro-Lee years of schooling and Mincerian returns) |

---

## Data Sources and Download Instructions

### Penn World Table 11.0

- **URL:** https://www.rug.nl/ggdc/productivity/pwt/
- **File to download:** `pwt110.xlsx` (Excel, version 11.0)
- **Instructions:** Navigate to the URL above, scroll to the latest
  version (11.0), click the Excel download link. Open the file and
  use the `Data` sheet.
- **Reference:** Feenstra, R. C., Inklaar, R., & Timmer, M. P. (2015).
  The next generation of the Penn World Table. *American Economic Review*,
  *105*(10), 3150–3182.

### World Bank Land Area

- **URL:** https://data.worldbank.org/indicator/AG.LND.TOTL.K2
- **Indicator code:** `AG.LND.TOTL.K2` (Land area in square kilometres)
- **Instructions:** Navigate to the URL above, click **Download** and
  select **CSV**. The downloaded ZIP contains a data CSV with year
  columns labelled `YYYY [YRYYYY]`. Rename the file to `land_area.csv`.
- **Treatment:** Land area does not change year-to-year in the World Bank
  series. The first available non-missing annual value is used for each
  country as a time-invariant attribute.

---

## Panel Construction Steps

To reconstruct `pwt_empirical_panel.xlsx` from raw sources:

1. Download `pwt110.xlsx` from the PWT website.
2. Download `land_area.csv` from the World Bank (indicator AG.LND.TOTL.K2).
3. From PWT, compute:
   - `ly  = log(rgdpo / pop)`
   - `lk  = log(rnna / pop)`
   - `lh  = log(hc)`
4. From World Bank, extract one land area value per country code and compute
   `larea = log(land_area_sqkm)`.
5. Merge PWT and land area on `countrycode`. Drop observations missing any
   of `rgdpo`, `pop`, `rnna`, `hc`, or `land_area_sqkm`.
6. Retain years 1994–2023. Keep only countries with a complete 30-year record
   (1994 is the lag year; estimation uses 1995–2023). This yields 143 countries.
7. Sort by `countrycode` and `year`. Assign numeric `id` (1 to 143).
8. Construct `lylag = lag(ly, 1)` within each country.
9. Construct Mundlak averages over the estimation period:
   `xbar_lk = mean(lk)` and `xbar_lh = mean(lh)` within each country.
10. Drop year 1994 (lag year) from the estimation panel.
    Final panel: 143 × 29 = 4,147 observations.

---

## Notes

- The panel is **fully balanced**: every country is observed in every year
  from 1995 to 2023.
- Land area values are constant over time in this dataset. The World Bank
  series does contain small year-to-year revisions (boundary adjustments),
  but for the purpose of this study a single value per country is used.
- The human capital index in PWT 11.0 follows the methodology described in
  Feenstra, Inklaar and Timmer (2015), based on average years of schooling
  from Barro and Lee (2013) and returns to education from Mincerian regressions.
