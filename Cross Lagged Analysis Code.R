# Loneliness x housing instability -- modelling on the balanced HRS panel built using ("Final Longitudinal Dataset.csv").
# Authors: Bijeta Pradhan and Sandile Ndabezitha
# Contents:
#   A. Covariate-adjusted GEE, contemporaneous (housing ~ loneliness + covariates)
#   B. Covariate-adjusted cross-lagged GEE, both directions
#   C. loneliness x retirement-status interaction (+ simple effects by ret)
#   D. Random-Intercept Cross-Lagged Panel Model (RI-CLPM) via lavaan


library(tidyverse)
library(geepack)
library(broom)
library(emmeans)
library(lavaan)

out_dir <- "analysis_output"
dir.create(out_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load and prepare
# -----------------------------------------------------------------------------
df <- readr::read_csv("Final Longitudinal Dataset (full).csv", show_col_types = FALSE)

waves_order <- c(2014, 2016, 2018, 2020, 2022)

df <- df %>%
  filter(!is.na(loneliness), !is.na(housing_stability)) %>%
  mutate(
    # binary outcomes (numeric 0/1) for GEE and RI-CLPM
    lonely_bin   = as.integer(loneliness == 1),
    unstable_bin = as.integer(housing_stability == "unstable housing"),

    # labelled factors for readable model output
    loneliness = factor(if_else(lonely_bin == 1, "Lonely", "Not Lonely"),
                        levels = c("Not Lonely", "Lonely")),
    housing    = factor(if_else(unstable_bin == 1, "unstable housing", "stable housing"),
                        levels = c("stable housing", "unstable housing")),

    # covariates
    log_income     = log(total_income + 1),
    education_cat  = fct_relevel(factor(education_cat), "Less than high school"),
    marital_status = factor(marital_status),
    race4          = fct_relevel(factor(race4), "white"),
    hisp           = factor(hisp),
    ret            = factor(ret),                       # "retired" / "not retired"
    wave_num       = as.integer(factor(wave, levels = waves_order))
  ) %>%
  arrange(person_id, wave_num)                          # GEE needs clusters contiguous

# covariate set used in all adjusted models (edit here to add sex, depression, etc.)
covars <- c("age", "education_cat", "marital_status", "race4", "hisp", "log_income")

# build a formula safely (handles interaction terms like "loneliness*ret")
mk_form <- function(response, terms) {
  as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

# helper: tidy a geeglm into an OR table and save
save_gee <- function(model, file, label) {
  tab <- broom::tidy(model, conf.int = TRUE) %>%
    mutate(
      OR       = exp(estimate),
      OR_low   = exp(conf.low),
      OR_high  = exp(conf.high),
      across(c(estimate, std.error, OR, OR_low, OR_high), ~ round(.x, 3)),
      p.value  = ifelse(p.value < .001, "<0.001", as.character(round(p.value, 3)))
    ) %>%
    select(term, estimate, std.error, p.value, OR, OR_low, OR_high)
  readr::write_csv(tab, file.path(out_dir, file))
  cat("\n====", label, "====\n"); print(tab)
  invisible(tab)
}

# -----------------------------------------------------------------------------
# A. Adjusted GEE -- contemporaneous (same-wave)
#    Outcome: unstable housing ; Predictor: loneliness + covariates
# -----------------------------------------------------------------------------
fA_vars <- c("unstable_bin", "loneliness", covars, "person_id")
dA <- df %>% select(all_of(fA_vars)) %>% drop_na()

mA <- geeglm(
  mk_form("unstable_bin", c("loneliness", covars)),
  id = factor(person_id), data = dA,
  family = binomial, corstr = "exchangeable"
)
save_gee(mA, "A_gee_contemporaneous.csv",
         "A. Adjusted GEE: unstable housing ~ loneliness + covariates (same wave)")

# -----------------------------------------------------------------------------
# B. Adjusted cross-lagged GEE (both directions)
#    Lag predictors within person; adjust for prior value of the OUTCOME.
# -----------------------------------------------------------------------------
df_lag <- df %>%
  group_by(person_id) %>%
  arrange(wave_num, .by_group = TRUE) %>%
  mutate(
    lonely_lag   = lag(lonely_bin),
    unstable_lag = lag(unstable_bin)
  ) %>%
  ungroup()

# B1: later housing instability ~ prior loneliness, adjusting for prior housing
dB1 <- df_lag %>%
  select(unstable_bin, lonely_lag, unstable_lag, all_of(covars), person_id) %>%
  drop_na()
mB1 <- geeglm(
  mk_form("unstable_bin", c("lonely_lag", "unstable_lag", covars)),
  id = factor(person_id), data = dB1, family = binomial, corstr = "exchangeable"
)
save_gee(mB1, "B1_gee_house_from_prior_loneliness.csv",
         "B1. Housing instability ~ prior loneliness (+ prior housing + covariates)")

# B2: later loneliness ~ prior housing instability, adjusting for prior loneliness
dB2 <- df_lag %>%
  select(lonely_bin, unstable_lag, lonely_lag, all_of(covars), person_id) %>%
  drop_na()
mB2 <- geeglm(
  mk_form("lonely_bin", c("unstable_lag", "lonely_lag", covars)),
  id = factor(person_id), data = dB2, family = binomial, corstr = "exchangeable"
)
save_gee(mB2, "B2_gee_loneliness_from_prior_house.csv",
         "B2. Loneliness ~ prior housing instability (+ prior loneliness + covariates)")

# -----------------------------------------------------------------------------
# C. loneliness x retirement-status interaction (contemporaneous housing outcome)
# -----------------------------------------------------------------------------
dC <- df %>%
  select(unstable_bin, loneliness, ret, all_of(covars), person_id) %>%
  drop_na()

mC <- geeglm(
  mk_form("unstable_bin", c("loneliness*ret", covars)),
  id = factor(person_id), data = dC, family = binomial, corstr = "exchangeable"
)
save_gee(mC, "C_gee_loneliness_by_ret_interaction.csv",
         "C. unstable housing ~ loneliness * ret + covariates")

# Simple effect of loneliness WITHIN each retirement group (as odds ratios)
simple_or <- tryCatch({
  emm <- emmeans(mC, ~ loneliness | ret)
  ctr <- summary(pairs(emm, reverse = TRUE), infer = TRUE)  # Lonely - Not Lonely on logit
  ctr %>%
    as.data.frame() %>%
    mutate(
      OR      = round(exp(estimate), 3),
      OR_low  = round(exp(asymp.LCL), 3),
      OR_high = round(exp(asymp.UCL), 3),
      p.value = ifelse(p.value < .001, "<0.001", as.character(round(p.value, 3)))
    ) %>%
    select(ret, contrast, estimate, p.value, OR, OR_low, OR_high)
}, error = function(e) {
  message("emmeans simple-effects step failed: ", conditionMessage(e))
  NULL
})
if (!is.null(simple_or)) {
  readr::write_csv(simple_or, file.path(out_dir, "C_simple_effects_by_ret.csv"))
  cat("\n==== C. Simple effect of loneliness on housing within ret group ====\n")
  print(simple_or)
}

# -----------------------------------------------------------------------------
# D. RI-CLPM (Random-Intercept Cross-Lagged Panel Model) via lavaan
#    Wide format: one row per person, 5 waves of loneliness (L) and housing (H).
#    Primary fit treats the binaries as numeric (common, interpretable
#    approximation). A fully categorical version (ordered + WLSMV) is provided
#    as a commented alternative below.
# -----------------------------------------------------------------------------
df_wide <- df %>%
  transmute(person_id, wave_num, L = lonely_bin, H = unstable_bin) %>%
  pivot_wider(names_from = wave_num, values_from = c(L, H), names_sep = "") %>%
  rename_with(~ gsub("^L", "L", .x)) %>%      # columns: L1..L5, H1..H5
  select(person_id, paste0("L", 1:5), paste0("H", 1:5))

riclpm_model <- '
  # ----- random intercepts (stable between-person differences) -----
  RIlon =~ 1*L1 + 1*L2 + 1*L3 + 1*L4 + 1*L5
  RIhou =~ 1*H1 + 1*H2 + 1*H3 + 1*H4 + 1*H5

  # ----- within-person components -----
  wL1 =~ 1*L1
  wL2 =~ 1*L2
  wL3 =~ 1*L3
  wL4 =~ 1*L4
  wL5 =~ 1*L5
  wH1 =~ 1*H1
  wH2 =~ 1*H2
  wH3 =~ 1*H3
  wH4 =~ 1*H4
  wH5 =~ 1*H5

  # observed residual variances fixed to 0 (variance carried by within factors)
  L1 ~~ 0*L1; L2 ~~ 0*L2; L3 ~~ 0*L3; L4 ~~ 0*L4; L5 ~~ 0*L5
  H1 ~~ 0*H1; H2 ~~ 0*H2; H3 ~~ 0*H3; H4 ~~ 0*H4; H5 ~~ 0*H5

  # ----- autoregressive (a, c) and cross-lagged (b, d) paths -----
  #   b = prior housing -> later loneliness ; d = prior loneliness -> later housing
  #   paths held time-invariant; remove the labels to free them across waves
  wL2 ~ a*wL1 + b*wH1
  wL3 ~ a*wL2 + b*wH2
  wL4 ~ a*wL3 + b*wH3
  wL5 ~ a*wL4 + b*wH4

  wH2 ~ c*wH1 + d*wL1
  wH3 ~ c*wH2 + d*wL2
  wH4 ~ c*wH3 + d*wL3
  wH5 ~ c*wH4 + d*wL4

  # ----- within-wave (residual) covariances -----
  wL1 ~~ wH1
  wL2 ~~ wH2
  wL3 ~~ wH3
  wL4 ~~ wH4
  wL5 ~~ wH5

  # ----- random-intercept (co)variances -----
  RIlon ~~ RIlon
  RIhou ~~ RIhou
  RIlon ~~ RIhou

  # random intercepts uncorrelated with the wave-1 within components
  RIlon ~~ 0*wL1 + 0*wH1
  RIhou ~~ 0*wL1 + 0*wH1
'

fit_riclpm <- tryCatch(
  lavaan(riclpm_model, data = df_wide, missing = "ML",
         meanstructure = TRUE, int.ov.free = TRUE),
  error = function(e) { message("RI-CLPM failed: ", conditionMessage(e)); NULL }
)

if (!is.null(fit_riclpm)) {
  cat("\n==== D. RI-CLPM fit measures ====\n")
  print(fitMeasures(fit_riclpm,
        c("npar","chisq","df","pvalue","cfi","tli","rmsea","srmr")))

  est <- parameterEstimates(fit_riclpm, standardized = TRUE) %>%
    filter(label %in% c("a","b","c","d") | op == "~~") %>%
    select(lhs, op, rhs, label, est, se, z, pvalue, std.all)
  readr::write_csv(est, file.path(out_dir, "D_riclpm_estimates.csv"))

  cat("\n==== D. RI-CLPM key paths (a=lon AR, c=hou AR, b=hou->lon, d=lon->hou) ====\n")
  print(est %>% filter(label %in% c("a","b","c","d")))
}

# ----- Fully categorical RI-CLPM (alternative; binary indicators) -------------
# If you prefer to treat L/H as ordered/binary rather than numeric, declare the
# indicators ordered and use the WLSMV estimator (theta parameterization).
# This is more correct for binary data but is more prone to non-convergence:
#
# fit_riclpm_cat <- lavaan(
#   riclpm_model, data = df_wide,
#   ordered = c(paste0("L", 1:5), paste0("H", 1:5)),
#   parameterization = "theta", estimator = "WLSMV", missing = "pairwise"
# )

cat("\nAll outputs written to: ", normalizePath(out_dir), "\n")
