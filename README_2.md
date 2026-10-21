# Cigarette Smoking Status and All Cause Mortality Among US Adults (NHANES 2009 to 2018)

Survey weighted survival analysis of cigarette smoking status (never, former, current) and all cause mortality among US adults aged 20 years or older, using the National Health and Nutrition Examination Survey (NHANES) 2009 to 2018 linked to the NCHS 2019 Public Use Linked Mortality File.

**Author:** Satwika Muddasani, MPH
**Preprint:** [add medRxiv link once posted]

## Data

All data are public and de identified. The script downloads them automatically from CDC:

* NHANES demographics (DEMO), smoking (SMQ) and body measures (BMX), cycles F to J
* NCHS Public Use Linked Mortality Files, follow up through December 31, 2019

No data files are stored in this repository.

## Methods in brief

* Exposure: smoking status from SMQ020 (100 cigarettes in lifetime) and SMQ040 (smoke now)
* Outcome: all cause mortality; time from MEC exam to death or end of follow up
* Design: SDMVPSU, SDMVSTRA and MEC weights divided by 5 for pooled cycles (survey package)
* Analyses: weighted Kaplan Meier curves, Rao Scott log rank test, survey weighted Cox models
  * Model 1: age, sex
  * Model 2: plus race and ethnicity, education, income to poverty ratio
  * Model 3: plus BMI category
* Proportional hazards checked with Schoenfeld residuals (cox.zph)
* Sensitivity analyses: excluding the first 12 months of follow up; ever vs never smoking
* Replication of an earlier unweighted SAS model on NHANES 2017 to 2018

## How to run

1. Install R (4.2 or newer) and RStudio, or use Posit Cloud.
2. Open `analysis.R` and click **Source**. First run downloads about 60 MB.
3. Results appear in `results/`:
   * `participant_flow.csv`
   * `table1_weighted.csv`
   * `table2_cox_models.csv`
   * `table3_sensitivity.csv`
   * `figure1_km_curves.png`
   * `ph_test.txt`
   * `results_summary.txt`

## Reporting

The manuscript follows the STROBE guidelines for cohort studies.

## License

Code released under the MIT License.
