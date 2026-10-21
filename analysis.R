# ==============================================================================
# Cigarette smoking status and all cause mortality among US adults
# NHANES 2009 to 2018 linked to the NCHS 2019 Public Use Linked Mortality File
#
# Author : Satwika Muddasani, MPH
# Purpose: Survey weighted survival analysis (Kaplan Meier and Cox models)
# Run in : RStudio or Posit Cloud (https://posit.cloud). Needs internet access.
# Output : results/ folder (tables, figure, results_summary.txt)
# ==============================================================================

# ---- 0. Packages --------------------------------------------------------------
pkgs <- c("haven", "dplyr", "survey", "survival", "readr")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

options(survey.lonely.psu = "adjust")   # standard NHANES setting
set.seed(2026)
dir.create("data",    showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---- 1. Download helpers ------------------------------------------------------
cycles <- data.frame(
  years  = c("2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"),
  start  = c(2009, 2011, 2013, 2015, 2017),
  suffix = c("F", "G", "H", "I", "J"),
  stringsAsFactors = FALSE
)

# CDC moved NHANES files in 2024; try the new path first, then the old one.
get_xpt <- function(file, start, years) {
  dest <- file.path("data", paste0(file, ".xpt"))
  if (!file.exists(dest)) {
    urls <- c(
      sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.xpt", start, file),
      sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s.XPT", years, file)
    )
    ok <- FALSE
    for (u in urls) {
      ok <- tryCatch({ download.file(u, dest, mode = "wb", quiet = TRUE); TRUE },
                     error = function(e) FALSE, warning = function(w) FALSE)
      if (ok && file.size(dest) > 1000) break
    }
    if (!ok) stop("Could not download ", file, ". Download it by hand into data/.")
  }
  haven::read_xpt(dest)
}

# NCHS linked mortality file (fixed width). Layout from the NCHS R read in program.
get_mort <- function(years) {
  f    <- sprintf("NHANES_%s_MORT_2019_PUBLIC.dat", sub("-", "_", years))
  dest <- file.path("data", f)
  if (!file.exists(dest)) {
    download.file(paste0("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/",
                         "datalinkage/linked_mortality/", f),
                  dest, mode = "wb", quiet = TRUE)
  }
  readr::read_fwf(
    dest,
    col_positions = readr::fwf_cols(
      SEQN = c(1, 6), ELIGSTAT = c(15, 15), MORTSTAT = c(16, 16),
      UCOD_LEADING = c(17, 19), PERMTH_INT = c(43, 45), PERMTH_EXM = c(46, 48)
    ),
    col_types = "iiiiii", na = c("", ".")
  )
}

# ---- 2. Build the pooled analytic file ---------------------------------------
one_cycle <- function(i) {
  s <- cycles$suffix[i]; y <- cycles$years[i]; st <- cycles$start[i]
  demo <- get_xpt(paste0("DEMO_", s), st, y) |>
    select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, DMDEDUC2, INDFMPIR,
           SDMVPSU, SDMVSTRA, WTMEC2YR)
  smq  <- get_xpt(paste0("SMQ_", s), st, y) |> select(SEQN, SMQ020, SMQ040)
  bmx  <- get_xpt(paste0("BMX_", s), st, y) |> select(SEQN, BMXBMI)
  mort <- get_mort(y)
  demo |>
    left_join(smq,  by = "SEQN") |>
    left_join(bmx,  by = "SEQN") |>
    left_join(mort, by = "SEQN") |>
    mutate(cycle = y)
}

raw <- bind_rows(lapply(seq_len(nrow(cycles)), one_cycle))

dat <- raw |>
  mutate(
    # Smoking status. Never = fewer than 100 cigarettes in lifetime.
    smoke3 = case_when(
      SMQ020 == 2                    ~ "Never",
      SMQ020 == 1 & SMQ040 == 3      ~ "Former",
      SMQ020 == 1 & SMQ040 %in% 1:2  ~ "Current",
      TRUE                           ~ NA_character_
    ),
    smoke3 = factor(smoke3, levels = c("Never", "Former", "Current")),
    ever   = factor(ifelse(smoke3 == "Never", "Never", "Ever"),
                    levels = c("Never", "Ever")),
    sex    = factor(RIAGENDR, 1:2, c("Male", "Female")),
    race   = factor(RIDRETH1, 1:5, c("Mexican American", "Other Hispanic",
                                     "Non Hispanic White", "Non Hispanic Black",
                                     "Other or Multiracial")),
    race   = relevel(race, ref = "Non Hispanic White"),
    educ   = case_when(DMDEDUC2 %in% 1:2 ~ "Less than high school",
                       DMDEDUC2 == 3     ~ "High school or GED",
                       DMDEDUC2 %in% 4:5 ~ "Some college or more"),
    educ   = factor(educ, c("Some college or more", "High school or GED",
                            "Less than high school")),
    pir_cat = cut(INDFMPIR, c(-Inf, 1.3, 3.5, Inf), right = FALSE,
                  labels = c("Below 1.3", "1.3 to 3.49", "3.5 or higher")),
    bmi_cat = cut(BMXBMI, c(-Inf, 18.5, 25, 30, Inf), right = FALSE,
                  labels = c("Under 18.5", "18.5 to 24.9", "25 to 29.9", "30 or higher")),
    bmi_cat = relevel(bmi_cat, ref = "18.5 to 24.9"),
    time_m  = PERMTH_EXM,                  # months from MEC exam to death or 12/31/2019
    time_y  = PERMTH_EXM / 12,
    death   = as.integer(MORTSTAT == 1),
    wt5     = WTMEC2YR / 5,                # 5 pooled cycles
    in_analysis = RIDAGEYR >= 20 & ELIGSTAT == 1 & !is.na(smoke3) &
                  !is.na(time_m) & !is.na(death) & WTMEC2YR > 0
  )

# Complete case flag for the fully adjusted model
dat <- dat |>
  mutate(in_analysis = in_analysis %in% TRUE,
         in_full = in_analysis & !is.na(race) & !is.na(educ) &
                   !is.na(pir_cat) & !is.na(bmi_cat))

# Flow of participants (report these in the Methods / Figure 1)
flow <- tibble::tibble(
  step = c("All NHANES participants 2009 to 2018",
           "Aged 20 years or older",
           "Eligible for mortality linkage",
           "Valid smoking status and follow up (analytic sample)",
           "Complete covariates (fully adjusted model)"),
  n = c(nrow(dat),
        sum(dat$RIDAGEYR >= 20, na.rm = TRUE),
        sum(dat$RIDAGEYR >= 20 & dat$ELIGSTAT == 1, na.rm = TRUE),
        sum(dat$in_analysis, na.rm = TRUE),
        sum(dat$in_full, na.rm = TRUE))
)
write.csv(flow, "results/participant_flow.csv", row.names = FALSE)

# ---- 3. Survey design (create on full data, THEN subset) ---------------------
des_all <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt5,
                     nest = TRUE, data = dat |> mutate(wt5 = ifelse(is.na(wt5), 0, wt5)))
des  <- subset(des_all, in_analysis)
desF <- subset(des_all, in_full)

# ---- 4. Table 1: weighted characteristics by smoking status ------------------
cat_row <- function(var) {
  tab <- svytable(as.formula(paste0("~", var, " + smoke3")), des)
  pc  <- round(100 * prop.table(tab, 2), 1)
  ov  <- round(100 * prop.table(svytable(as.formula(paste0("~", var)), des)), 1)
  data.frame(variable = var, level = rownames(pc), Overall = as.numeric(ov),
             Never = pc[, "Never"], Former = pc[, "Former"], Current = pc[, "Current"],
             row.names = NULL)
}
num_row <- function(var) {
  f  <- as.formula(paste0("~", var))
  tb <- svyby(f, ~smoke3, des, svymean, na.rm = TRUE)
  ov <- svymean(f, des, na.rm = TRUE)
  data.frame(variable = var, level = "mean (SE)",
             Overall = sprintf("%.1f (%.2f)", coef(ov), SE(ov)),
             t(sprintf("%.1f (%.2f)", tb[[var]], tb[[paste0("se")]])) |>
               `colnames<-`(levels(dat$smoke3)), check.names = FALSE)
}
n_unw <- table(dat$smoke3[dat$in_analysis])
deaths_unw <- tapply(dat$death[dat$in_analysis], dat$smoke3[dat$in_analysis], sum)
w_pct <- round(100 * coef(svymean(~smoke3, des)), 1)

table1 <- bind_rows(
  data.frame(variable = "N (unweighted)", level = "",
             Overall = as.character(sum(n_unw)), t(as.character(n_unw)) |>
               `colnames<-`(names(n_unw)), check.names = FALSE),
  data.frame(variable = "Weighted percent", level = "", Overall = "100",
             t(as.character(w_pct)) |> `colnames<-`(levels(dat$smoke3)),
             check.names = FALSE),
  data.frame(variable = "Deaths (unweighted)", level = "",
             Overall = as.character(sum(deaths_unw)), t(as.character(deaths_unw)) |>
               `colnames<-`(names(deaths_unw)), check.names = FALSE),
  num_row("RIDAGEYR"),
  cat_row("sex") |> mutate(across(everything(), as.character)),
  cat_row("race") |> mutate(across(everything(), as.character)),
  cat_row("educ") |> mutate(across(everything(), as.character)),
  cat_row("pir_cat") |> mutate(across(everything(), as.character)),
  cat_row("bmi_cat") |> mutate(across(everything(), as.character))
)
write.csv(table1, "results/table1_weighted.csv", row.names = FALSE)

fu <- svyquantile(~time_y, des, 0.5, na.rm = TRUE)
fu_median <- round(coef(fu)[1], 1)
person_years <- round(sum(dat$time_y[dat$in_analysis]))

# ---- 5. Weighted Kaplan Meier curves -----------------------------------------
km <- svykm(Surv(time_y, death) ~ smoke3, design = des)
cols <- c("#1b9e77", "#7570b3", "#d95f02")
png("results/figure1_km_curves.png", width = 2000, height = 1500, res = 300)
par(mar = c(4.5, 4.5, 1, 1))
plot(NA, xlim = c(0, 11.5), ylim = c(0.80, 1), las = 1,
     xlab = "Years of follow up", ylab = "Survival probability (weighted)")
for (i in 1:3) lines(km[[i]]$time, km[[i]]$surv, type = "s", col = cols[i], lwd = 2)
legend("bottomleft", legend = levels(dat$smoke3), col = cols, lwd = 2, bty = "n")
dev.off()

# Weighted log rank style test (Rao Scott)
lr <- svylogrank(Surv(time_y, death) ~ smoke3, design = des)

# ---- 6. Survey weighted Cox models -------------------------------------------
m1 <- svycoxph(Surv(time_y, death) ~ smoke3 + RIDAGEYR + sex, design = des)
m2 <- svycoxph(Surv(time_y, death) ~ smoke3 + RIDAGEYR + sex + race + educ + pir_cat,
               design = desF)
m3 <- svycoxph(Surv(time_y, death) ~ smoke3 + RIDAGEYR + sex + race + educ + pir_cat +
                 bmi_cat, design = desF)

hr_tab <- function(m, label) {
  ci <- exp(confint(m)); est <- exp(coef(m)); p <- summary(m)$coefficients[, ncol(summary(m)$coefficients)]
  r  <- grep("^smoke3", names(est))
  data.frame(model = label, term = sub("smoke3", "", names(est)[r]),
             HR = round(est[r], 2), LCL = round(ci[r, 1], 2), UCL = round(ci[r, 2], 2),
             p = signif(p[r], 3), n = m$n, events = m$nevent, row.names = NULL)
}
cox_tab <- bind_rows(
  hr_tab(m1, "Model 1: age, sex"),
  hr_tab(m2, "Model 2: + race and ethnicity, education, income to poverty ratio"),
  hr_tab(m3, "Model 3: + BMI category")
)
write.csv(cox_tab, "results/table2_cox_models.csv", row.names = FALSE)

# Current vs former smokers (Model 3, reference = Former)
desF2 <- update(desF, smoke3b = relevel(smoke3, ref = "Former"))
m3b <- svycoxph(Surv(time_y, death) ~ smoke3b + RIDAGEYR + sex + race + educ + pir_cat +
                  bmi_cat, design = desF2)
cur_vs_former <- exp(cbind(HR = coef(m3b), confint(m3b)))["smoke3bCurrent", ]

# ---- 7. Proportional hazards check -------------------------------------------
# Weighted coxph with weights normalized to mean 1 (raw survey weights make the
# Schoenfeld test meaningless), robust SEs, same covariates as Model 3.
d_ph <- dat[dat$in_full, ]
d_ph$w <- d_ph$wt5 / mean(d_ph$wt5)
fit_ph <- coxph(Surv(time_y, death) ~ smoke3 + RIDAGEYR + sex + race + educ + pir_cat +
                  bmi_cat, data = d_ph, weights = w, robust = TRUE)
ph <- cox.zph(fit_ph)
capture.output(print(ph), file = "results/ph_test.txt")

# ---- 8. Sensitivity analyses -------------------------------------------------
# (a) Exclude deaths and censoring in the first 12 months (reverse causation)
desS <- subset(desF, time_m > 12)
s1 <- svycoxph(Surv(time_y, death) ~ smoke3 + RIDAGEYR + sex + race + educ + pir_cat +
                 bmi_cat, design = desS)
# (b) Ever vs never smoker, fully adjusted
s2 <- svycoxph(Surv(time_y, death) ~ ever + RIDAGEYR + sex + race + educ + pir_cat +
                 bmi_cat, design = desF)
# (c) Attained age as the time scale (age did not meet the PH assumption)
desA <- update(desF, age_in = RIDAGEYR, age_out = RIDAGEYR + pmax(time_y, 1 / 24))
s3 <- svycoxph(Surv(age_in, age_out, death) ~ smoke3 + sex + race + educ + pir_cat +
                 bmi_cat, design = desA)
sens <- bind_rows(
  hr_tab(s1, "Sensitivity: excluding first 12 months"),
  hr_tab(s3, "Sensitivity: attained age as time scale"),
  data.frame(model = "Sensitivity: ever vs never (Model 3 covariates)", term = "Ever",
             HR = round(exp(coef(s2))["everEver"], 2),
             LCL = round(exp(confint(s2))["everEver", 1], 2),
             UCL = round(exp(confint(s2))["everEver", 2], 2),
             p = signif(summary(s2)$coefficients["everEver", ncol(summary(s2)$coefficients)], 3),
             n = s2$n, events = s2$nevent, row.names = NULL)
)
write.csv(sens, "results/table3_sensitivity.csv", row.names = FALSE)

# ---- 9. Replication of the original SAS analysis (2017 to 2018 only) ---------
# Unweighted, ever vs never, age and sex adjusted, interview based follow up.
rep <- dat |> filter(cycle == "2017-2018", RIDAGEYR >= 20, ELIGSTAT == 1, !is.na(ever))
r1  <- coxph(Surv(PERMTH_INT, death) ~ ever + RIDAGEYR + sex, data = rep)
rep_hr <- exp(cbind(HR = coef(r1), confint(r1)))["everEver", ]

# ---- 10. Plain text summary to copy into the manuscript ----------------------
fmt <- function(r) sprintf("HR %.2f (95%% CI %.2f to %.2f)", r$HR, r$LCL, r$UCL)
sink("results/results_summary.txt")
cat("NHANES 2009 to 2018 smoking and all cause mortality: RESULTS SUMMARY\n")
cat("Generated:", format(Sys.time()), " | R", R.version$major, ".", R.version$minor,
    " | survey", as.character(packageVersion("survey")), "\n\n")
cat("PARTICIPANT FLOW\n"); print(flow); cat("\n")
cat("Analytic sample (unweighted N):", sum(n_unw), "\n")
cat("Deaths:", sum(deaths_unw), "\n")
cat("Person years of follow up:", person_years, "\n")
cat("Weighted median follow up (years):", fu_median, "\n")
cat("Weighted prevalence (%): Never", w_pct[1], "| Former", w_pct[2], "| Current", w_pct[3], "\n\n")
cat("RAO SCOTT LOG RANK TEST\n"); print(lr); cat("\n")
cat("COX MODELS (reference = Never)\n"); print(cox_tab, row.names = FALSE); cat("\n")
for (i in seq_len(nrow(cox_tab)))
  cat(cox_tab$model[i], "|", cox_tab$term[i], ":", fmt(cox_tab[i, ]), "\n")
cat("\nCurrent vs Former (Model 3): HR", round(cur_vs_former[1], 2),
    "(95% CI", round(cur_vs_former[2], 2), "to", round(cur_vs_former[3], 2), ")\n\n")
cat("SENSITIVITY\n"); print(sens, row.names = FALSE); cat("\n")
cat("PROPORTIONAL HAZARDS (global p should be > 0.05)\n"); print(ph); cat("\n")
cat("REPLICATION OF SAS RESULT (2017 to 2018, unweighted, age and sex adjusted)\n")
cat("N =", r1$n, "| deaths =", r1$nevent, "| HR ever vs never",
    sprintf("%.2f (95%% CI %.2f to %.2f)", rep_hr[1], rep_hr[2], rep_hr[3]), "\n")
cat("Original SAS: N = 5,524, HR 1.78 (95% CI 1.25 to 2.53), p = 0.001\n")
sink()

cat("\nDone. Open results/results_summary.txt and copy the numbers into the manuscript.\n")
sessionInfo()
