# Cigarette Smoking Status and All Cause Mortality Among US Adults (NHANES 2009 to 2018)

Survey weighted survival analysis of cigarette smoking status (never, former, current) and all cause mortality among US adults aged 20 years or older, using the National Health and Nutrition Examination Survey (NHANES) 2009 to 2018 linked to the NCHS 2019 Public Use Linked Mortality Files.

**Author:** Satwika Muddasani, MPH ([ORCID 0009-0002-4850-6863](https://orcid.org/0009-0002-4850-6863))

## Data

The data are publicly available from CDC, and the code downloads them directly from the official source:

* NHANES Demographics, Smoking and Cigarette Use, and Body Measures files, cycles 2009 to 2010 through 2017 to 2018
* NCHS 2019 Public Use Linked Mortality Files, with follow up through December 31, 2019

No data files are stored in this repository.

## Methods in brief

* Exposure: smoking status from SMQ020 (smoked 100 cigarettes in lifetime) and SMQ040 (smoke now)
* Outcome: all cause mortality, with follow up from the MEC exam to death or December 31, 2019
* Survey design: strata (SDMVSTRA), clusters (SDMVPSU), and MEC weights divided by 5 for the pooled cycles
* Weighted Kaplan Meier curves and Rao Scott log rank test
* Survey weighted Cox models, all fit in the same complete case sample:
  * Model 1: age and sex
  * Model 2: plus race and ethnicity, education, and income to poverty ratio
  * Model 3: plus BMI category
* Proportional hazards checked with Schoenfeld residuals
* Sensitivity analyses: excluding the first 12 months of follow up, ever vs never smoking, and attained age as the time scale

## Work completed

The R script (`analysis.R`) builds the analytic dataset and runs the full analysis from publicly available source files:

1. **Retrieve and pool source files.** It downloads demographics, smoking, and body-measure XPT files for the five NHANES two-year cycles from 2009–2010 through 2017–2018, along with each cycle's NCHS linked mortality file. It joins them by participant identifier and pools the cycles.
2. **Define the study variables and sample.** Smoking status is classified as never, former, or current from the NHANES cigarette-use questions. The script derives demographic and socioeconomic covariates, BMI categories, follow-up time, and mortality status. The primary analytic sample is adults aged 20 or older who are eligible for mortality linkage and have valid smoking and follow-up data; the fully adjusted models use a common complete-case sample.
3. **Account for the survey design.** Analyses use NHANES examination weights, strata, and primary sampling units, with the two-year MEC weights divided by five for the pooled cycles.
4. **Estimate associations.** Outputs include weighted participant characteristics, Kaplan–Meier curves and a Rao–Scott log-rank test, and nested survey-weighted Cox models adjusting first for age and sex, then race/ethnicity, education and income-to-poverty ratio, and finally BMI.
5. **Assess assumptions and robustness.** The script checks proportional hazards using Schoenfeld residuals and runs prespecified sensitivity analyses: removing the first year of follow-up, comparing ever versus never smoking, and using attained age as the time scale.
6. **Save reproducible outputs.** Tables, a survival figure, model diagnostics, a text summary, and a cleaned one-row-per-participant analysis file are written to `results/`.

This is an observational analysis. The models estimate associations and do not establish that smoking status caused an individual's mortality outcome.

## How to run

1. Install R (version 4.2 or newer) and RStudio, or use Posit Cloud.
2. Open analysis.R and click Source. The first run downloads the NHANES and mortality files from CDC, so it may take 10 to 20 minutes.
3. The results are saved in the results folder:
   * participant_flow.csv: sample size at each step
   * table1_weighted.csv: participant characteristics
   * table2_cox_models.csv: hazard ratios from the Cox models
   * table3_sensitivity.csv: sensitivity analyses
   * figure1_km_curves.png: survival curves
   * ph_test.txt: proportional hazards check
   * results_summary.txt: all key numbers in one file
   * nhanes_smoking_mortality_clean.csv: the cleaned analytic dataset (one row per participant)

## Reporting

The analysis is reported following the STROBE guidelines for observational cohort studies. A manuscript based on this analysis is in preparation for journal submission.

## License

Code released under the MIT License.
