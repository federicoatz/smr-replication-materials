# Empirical anchoring of the safeguard-power inflation factor (Equation 2 /
# eq:markup in Paper 11 "APriori"/Wizard) against the SCORE replication
# dataset (Tyner et al. 2026, Nature 652:143; OSF osf.io/g5sny).
#
# N_safeguard / N_naive ~= (1 - z_gamma / |t_orig|)^(-2)
#
# t_orig prioritizes a directly reported t or z statistic, or the signed
# square root of a one-df F or chi-squared statistic. The Pearson-r conversion
# is used only as a fallback when no direct route is available.

suppressMessages({
  library(dplyr)
  library(pwr)
})

data_dir <- "../g5sny-osfstorage-archive/extracted"
load(file.path(data_dir, "analyst data.RData"))

alpha        <- .05   # significance standard for N_naive / N_safeguard
target_power <- .80   # power standard for N_naive / N_safeguard
gamma        <- .80   # safeguard confidence level (Wizard default)
z_gamma      <- qnorm(gamma)   # one-sided critical value, ~0.8416

## ---- 1. Claim-level analytic frame -----------------------------------

ro <- repli_outcomes %>%
  filter(!is_covid, repli_version_of_record) %>%
  select(claim_id, repli_type, repli_sample_size_value, repli_conv_r,
         repli_score_criteria_met)

oo <- orig_outcomes %>%
  semi_join(ro, by = "claim_id") %>%
  select(claim_id, orig_conv_r, orig_sample_size_value,
         orig_stat_type, orig_stat_value, orig_stat_dof_1)

d <- oo %>%
  inner_join(ro, by = "claim_id") %>%
  filter(!is.na(orig_conv_r), !is.na(repli_conv_r), !is.na(orig_sample_size_value))

cat("Claims with orig_conv_r, repli_conv_r and orig_sample_size_value all present:",
    nrow(d), "(Nature Fig. 4 / Table 2 report 249 claims with both r's; one",
    "additional claim here lacks orig_sample_size_value and cannot supply t_orig)\n")

## ---- 2. t_orig -----------------------------------------------------------
#
# t_orig = delta_hat / SE(delta_hat) is the decisiveness of the ORIGINAL
# test. Reconstructing it from (orig_conv_r, N) via the simple-correlation
# formula t = r*sqrt(N-2)/sqrt(1-r^2) implicitly assumes iid observations.
# That assumption fails badly for clustered/complex designs with very large
# N (survey panels, country-year data): the r-conversion then predicts a
# huge, spurious t, because it cannot see the design effect that inflated
# the true SE. Direct evidence for this: on the 139 claims with a reported
# t-statistic, agreement with the r-conversion is tight for typical N
# (median |diff| = 0.05) but breaks down completely for the largest-N
# claims, dragging the Pearson correlation on |t| down to 0.25 while the
# rank (Spearman) correlation stays at 0.68 -- a classic outlier signature.
#
# Fix: use the DIRECTLY reported test statistic as t_orig whenever the
# original paper's own inferential test already delivers one (t, z, or an
# F/chi-squared with 1 numerator df, which is algebraically |t| or |z|).
# The r-conversion is used only as a fallback for statistic types that
# provide no such direct route (multi-df F, multi-df chi-squared,
# non-parametric statistics, delta-gamma-squared).

d <- d %>%
  mutate(
    t_from_stat = case_when(
      orig_stat_type %in% c("t", "z")                        ~ orig_stat_value,
      orig_stat_type == "F" & orig_stat_dof_1 == 1            ~ sqrt(orig_stat_value),
      orig_stat_type == "chi_squared" & orig_stat_dof_1 == 1  ~ sqrt(orig_stat_value),
      TRUE                                                    ~ NA_real_
    ),
    t_from_stat = sign(orig_conv_r) * abs(t_from_stat),
    t_from_r    = orig_conv_r * sqrt(orig_sample_size_value - 2) /
                  sqrt(1 - orig_conv_r^2),
    t_orig      = coalesce(t_from_stat, t_from_r),
    t_source    = ifelse(!is.na(t_from_stat), "direct_stat", "r_conversion")
  )

cat("\nt_orig source:\n")
print(table(d$t_source))

n_for_r <- function(r) {
  r <- abs(r)
  if (is.na(r) || r <= 0 || r >= 1) return(NA_real_)
  min_n_power <- tryCatch(
    pwr.r.test(n = 4, r = r, sig.level = alpha,
               alternative = "two.sided")$power,
    error = function(e) NA_real_
  )
  if (!is.na(min_n_power) && min_n_power >= target_power) return(4)
  out <- tryCatch(
    pwr.r.test(r = r, sig.level = alpha, power = target_power,
               alternative = "two.sided")$n,
    error = function(e) NA_real_
  )
  if (is.na(out)) NA_real_ else ceiling(out)
}

d <- d %>%
  mutate(
    abs_t             = abs(t_orig),
    safeguard_defined = abs_t > z_gamma,
    markup_closedform = ifelse(safeguard_defined, (1 - z_gamma / abs_t)^(-2), NA_real_),
    r_sg              = ifelse(safeguard_defined,
                                sign(orig_conv_r) * abs(orig_conv_r) * (1 - z_gamma / abs_t),
                                NA_real_)
  )

d$N_naive             <- vapply(d$orig_conv_r, n_for_r, numeric(1))
d$N_safeguard_exact    <- vapply(d$r_sg, n_for_r, numeric(1))
d$N_safeguard_closedform <- d$N_naive * d$markup_closedform

d <- d %>%
  mutate(
    ratio_obs_safeguard = repli_sample_size_value / N_safeguard_exact,
    met_safeguard        = ratio_obs_safeguard >= 1
  )

## ---- 3. Robustness: where the fallback r-conversion is used -----------
#
# For the 139 claims with orig_stat_type == "t" we have ground truth (the
# reported t itself), so we can bound how badly the r-conversion would have
# done had it been used there instead of the direct statistic.

direct_t <- d %>%
  filter(orig_stat_type == "t") %>%
  mutate(t_abs_diff = abs(t_from_r - orig_stat_value))

cat("\n== Fallback check: r-conversion vs. reported t, n =",
    nrow(direct_t), "==\n")
print(summary(abs(direct_t$t_abs_diff)))
cat("Pearson cor (r-conversion vs reported |t|):",
    round(cor(abs(direct_t$t_from_r), abs(direct_t$orig_stat_value)), 3), "\n")
cat("Spearman cor:",
    round(cor(abs(direct_t$t_from_r), abs(direct_t$orig_stat_value),
              method = "spearman"), 3), "\n")
cat("Claims using the r-conversion fallback in the main analysis:",
    sum(d$t_source == "r_conversion"), "of", nrow(d), "\n")

## ---- 4. Primary sample: new-data replications --------------------------

new_data <- d %>% filter(repli_type == "new data")

cat("\n== New-data replications (n =", nrow(new_data), ") ==\n")
cat("Safeguard bound undefined (|t_orig| <= z_gamma):",
    sum(!new_data$safeguard_defined), "of", nrow(new_data), "\n")

nd_valid <- new_data %>% filter(safeguard_defined, !is.na(N_safeguard_exact))
cat("Median numerically solved markup N_safeguard/N_naive:",
    round(median(nd_valid$N_safeguard_exact / nd_valid$N_naive,
                 na.rm = TRUE), 3), "\n")
cat("Median ratio N_observed/N_safeguard:",
    round(median(nd_valid$ratio_obs_safeguard, na.rm = TRUE), 3), "\n")
cat("Share of replications planned BELOW the safeguard criterion:",
    round(mean(nd_valid$ratio_obs_safeguard < 1, na.rm = TRUE) * 100, 1), "%\n")

## ---- 5. Validation test -------------------------------------------------

cat("\n== Validation: success rate by safeguard-criterion status (new data) ==\n")
val_tab <- table(nd_valid$met_safeguard, nd_valid$repli_score_criteria_met,
                  dnn = c("met_safeguard", "success"))
print(val_tab)

if (all(dim(val_tab) == c(2, 2))) {
  vt <- prop.test(c(val_tab["TRUE", "TRUE"], val_tab["FALSE", "TRUE"]),
                   c(sum(val_tab["TRUE", ]), sum(val_tab["FALSE", ])))
  cat("Success rate | met safeguard:   ",
      round(val_tab["TRUE", "TRUE"] / sum(val_tab["TRUE", ]) * 100, 1), "%\n")
  cat("Success rate | below safeguard: ",
      round(val_tab["FALSE", "TRUE"] / sum(val_tab["FALSE", ]) * 100, 1), "%\n")
  cat("Two-proportion test p-value:", signif(vt$p.value, 3), "\n")
}

## ---- 5b. Does the shortfall concentrate on fragile originals? ---------
#
# The paper's actual claim is not "replications are underpowered on
# average" but "the shortfall, where it exists, tracks how fragile the
# original evidence was" (low |t_orig| -> high required markup). Test that
# directly: correlate |t_orig| (fragility, inverted) and the closed-form
# markup against the realized ratio_obs_safeguard.

cat("\n== Does required markup track the realized shortfall? (new data) ==\n")
cat("cor(abs_t, ratio_obs_safeguard), Spearman:",
    round(cor(nd_valid$abs_t, nd_valid$ratio_obs_safeguard,
              method = "spearman", use = "complete.obs"), 3), "\n")
cat("cor(markup_closedform, ratio_obs_safeguard), Spearman:",
    round(cor(nd_valid$markup_closedform, nd_valid$ratio_obs_safeguard,
              method = "spearman", use = "complete.obs"), 3), "\n")

nd_valid <- nd_valid %>%
  mutate(fragility_tercile = ntile(abs_t, 3))
cat("\nMedian ratio_obs_safeguard by original-evidence fragility tercile",
    "(1 = most fragile / lowest |t_orig|):\n")
print(nd_valid %>% group_by(fragility_tercile) %>%
        summarise(n = n(),
                   median_abs_t = round(median(abs_t), 2),
                   median_markup = round(median(markup_closedform), 2),
                   median_ratio = round(median(ratio_obs_safeguard), 2),
                   share_below_1 = round(mean(ratio_obs_safeguard < 1) * 100, 1)))

## ---- 6. Secondary-data contrast (structural, not chosen N) -------------

secondary_data <- d %>% filter(repli_type == "secondary data") %>%
  filter(safeguard_defined, !is.na(N_safeguard_exact))

cat("\n== Secondary-data replications (n =", nrow(secondary_data), ") ==\n")
cat("Median ratio N_available/N_safeguard:",
    round(median(secondary_data$ratio_obs_safeguard, na.rm = TRUE), 3), "\n")
cat("Share below the safeguard criterion:",
    round(mean(secondary_data$ratio_obs_safeguard < 1, na.rm = TRUE) * 100, 1), "%\n")

## ---- 7. Save claim-level output ----------------------------------------

out_path <- "safeguard_inflation_claim_level.csv"
write.csv(d, out_path, row.names = FALSE)
cat("\nSaved claim-level table to", out_path, "\n")
