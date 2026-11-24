# sap_flow

Sap flux density (Granier thermal dissipation probes, 2.5 cm) of the COMPASS synoptic trees,
OWC plots UP (upland) and TR (transition), 2022-2025.
Script: `R_scripts/processing_sapflow.R` 

## Inputs

- **Raw data (not in this repository):** COMPASS L2 v2-1 at `compass_dir`: sap flow
  `sapflow-3.5cm` (mV, 15 min; deep sensors `xxD` not used), weather `wx-rh15` and `wx-tempavg15`
  (used for VPD, kPa, `plantecophys::RHtoVPD`).
- **Design tables:** `Design/sapflow_daily_issues.csv` (bad days from the technician notes),
  `Design/sapflow_manual_checks.csv` (days removed after visual check).
- **Field data in this folder:** `LE_TreePlotInventory_OWC.csv` (tree inventory, DBH per year)
  and `bai_data_OWC.csv` (cores, sapwood depth). Not used by the script; kept for stand
  transpiration scaling. The inventory does not link sensors to trees, so all OWC sap flow trees
  are labelled Quercus.

## Processing steps and criteria

1. **Technician notes:** readings in periods with sensor or power problems are set to NA.
2. **dTmax (zero-flow dT), per sensor and night:**
   - Night window **21:00-06:00**, assigned to the following day. The COMPASS standard is
     23:00-04:00, but with that window 10-22% of readings at 05-06 h and 21-22 h exceeded dTmax,
     because dT keeps rising until dawn (zero flow lasts from ~21 h to ~6 h at OWC).
   - A night reading is accepted only if VPD < 0.2 kPa and it is not a spike (|Z| <= 2.5 within the night).
   - dTmax = maximum accepted reading of the night.
   - **Stability filter:** nights with dTmax more than 2 SD from the 15-day rolling mean are removed.
   - **Gap filling (`dTmax_final`):** nights without a valid dTmax repeat the last valid value
     until the criteria are met again.
3. **Sap flux density (Granier):** F = 0.011899 x (dTmax_final / dT - 1)^1.231 (cm3 cm-2 s-1).
4. **QA/QC cascade (rows removed):**
   1. Bad mV days: >= 10% of the day's readings below the sensor's 1st percentile of positive mV, and any F > 0.005.
   2. Sharp drop: reading < 10% of the median of the previous 2 h.
   3. Sustained fault: sensor daily median < 10% of the plot median (>= 3 sensors) for >= 2 consecutive days.
   4. dT < 0.1 mV, or F not finite.
   5. Days flagged in the visual check (`Design/sapflow_manual_checks.csv`).
   - The MAD outlier step of the original qmd is not applied (it removed real signal).
5. **Hourly mean per sensor.** F set to NA in Nov-Jan (sensors not maintained) and from
   14 Oct 2024 to 31 May 2025 (bad data).

**Open question (check with Nate):** readings with dT > dTmax_final (no flow) give non-finite F
and are currently **removed** in step 4.4 (2.3% of TR and 6.4% of UP valid readings). The
alternative is to set F = 0, which keeps these zeros in hourly and daily means.

## processing_outputs/

| File | Content |
|---|---|
| `sapflow_hourly_sensor.parquet` | datetime (Etc/GMT+5), Site, Plot, Sensor_ID, Species, F (cm3 cm-2 s-1), sapflow_mV, n_readings |

## Data coverage (study window 5 Apr 2024 - 14 Jan 2025)

TR: 6 sensors with data from April/May to 15 Sep-9 Oct 2024. UP: 3 sensors. No sap flow from
mid-Oct 2024 to Jan 2025 (leaf-off, transpiration ~0).
