# Direct test of whether replication planning responded to original-evidence
# fragility, plus the placebo (secondary-data) comparison, a selection check
# on the binary success reversal, a power/gamma sensitivity grid, and a check
# that the r-conversion fallback for t_orig is not tercile-concentrated.
#
# Builds on compute_safeguard_inflation.R; re-derives the claim-level frame
# so that N_naive and the safeguard markup can be recomputed at arbitrary
# (power, gamma) combinations for the sensitivity grid.

suppressMessages({
  library(dplyr)
  library(pwr)
  library(sandwich)
  library(lmtest)
  library(splines)
})

data_dir <- "../g5sny-osfstorage-archive/extracted"
load(file.path(data_dir, "analyst data.RData"))

alpha_sig <- .05  # significance standard, held fixed throughout

## ---- 1. Claim-level frame (as in compute_safeguard_inflation.R) --------

ro <- repli_outcomes %>%
  filter(!is_covid, repli_version_of_record) %>%
  select(claim_id, paper_id, repli_type, repli_sample_size_value,
         repli_conv_r, repli_score_criteria_met)

oo <- orig_outcomes %>%
  semi_join(ro, by = "claim_id") %>%
  select(claim_id, orig_conv_r, orig_sample_size_value,
         orig_stat_type, orig_stat_value, orig_stat_dof_1)

d <- oo %>%
  inner_join(ro, by = "claim_id") %>%
  filter(!is.na(orig_conv_r), !is.na(repli_conv_r), !is.na(orig_sample_size_value))

d <- d %>%
  mutate(
    t_from_stat = case_when(
      orig_stat_type %in% c("t", "z")                        ~ orig_stat_value,
      orig_stat_type == "F" & orig_stat_dof_1 == 1 & orig_stat_value >= 0 ~
        sqrt(pmax(orig_stat_value, 0)),
      orig_stat_type == "chi_squared" & orig_stat_dof_1 == 1 & orig_stat_value >= 0 ~
        sqrt(pmax(orig_stat_value, 0)),
      TRUE                                                    ~ NA_real_
    ),
    t_from_stat = sign(orig_conv_r) * abs(t_from_stat),
    t_from_r    = orig_conv_r * sqrt(orig_sample_size_value - 2) /
                  sqrt(1 - orig_conv_r^2),
    t_orig      = coalesce(t_from_stat, t_from_r),
    abs_t       = abs(t_orig),
    t_source    = ifelse(!is.na(t_from_stat), "direct_stat", "r_conversion")
  )

## discipline control: the same 6-discipline grouping Tyner et al. use
disc <- paper_metadata %>% select(paper_id, discipline = COS_pub_category) %>%
  distinct(paper_id, .keep_all = TRUE)
d <- d %>% left_join(disc, by = "paper_id")

n_for_r <- function(r, power, alpha = alpha_sig) {
  r <- abs(r)
  if (is.na(r) || r <= 0 || r >= 1) return(NA_real_)
  min_n_power <- tryCatch(
    pwr.r.test(n = 4, r = r, sig.level = alpha,
               alternative = "two.sided")$power,
    error = function(e) NA_real_)
  if (!is.na(min_n_power) && min_n_power >= power) return(4)
  out <- tryCatch(
    pwr.r.test(r = r, sig.level = alpha, power = power,
               alternative = "two.sided")$n,
    error = function(e) NA_real_
  )
  if (is.na(out)) NA_real_ else ceiling(out)
}

build_frame <- function(d, power, gamma) {
  z_gamma <- qnorm(gamma)
  d <- d %>%
    mutate(
      safeguard_defined = abs_t > z_gamma,
      markup_approx = ifelse(safeguard_defined,
                             (1 - z_gamma / abs_t)^(-2), NA_real_),
      r_sg   = ifelse(safeguard_defined,
                       sign(orig_conv_r) * abs(orig_conv_r) * (1 - z_gamma / abs_t),
                       NA_real_)
    )
  d$N_naive <- vapply(d$orig_conv_r, n_for_r, numeric(1), power = power)
  d$N_safeguard <- vapply(d$r_sg, n_for_r, numeric(1), power = power)
  d <- d %>% mutate(
    markup = N_safeguard / N_naive,
    ratio_obs_safeguard = repli_sample_size_value / N_safeguard)
  d
}

d80 <- build_frame(d, power = .80, gamma = .80)

## ---- 2. Main regression: did planning respond to fragility? -----------
## log(N_obs) = a + b1 log(N_naive) + b2 log(markup) + discipline FE

fit_elasticity <- function(dat, label) {
  m <- dat %>%
    filter(safeguard_defined, !is.na(N_naive), N_naive > 0,
           !is.na(repli_sample_size_value), repli_sample_size_value > 0,
           !is.na(discipline)) %>%
    mutate(logN_obs = log(repli_sample_size_value),
           logN_naive = log(N_naive),
           logmarkup = log(markup))
  discipline_counts <- table(m$discipline)
  m$discipline_collapsed <- ifelse(
    m$discipline %in% names(discipline_counts[discipline_counts >= 10]),
    as.character(m$discipline), "Other")
  fit <- lm(logN_obs ~ logN_naive + logmarkup + discipline_collapsed, data = m)
  rob <- coeftest(fit, vcov = vcovHC(fit, type = "HC1"))
  cat("\n==", label, "(n =", nrow(m), ") ==\n")
  print(rob)
  list(fit = fit, rob = rob, n = nrow(m), data = m)
}

cat("\n########## MAIN TEST: new-data replications ##########\n")
new_fit <- fit_elasticity(d80 %>% filter(repli_type == "new data"), "New data")

cat("\n########## PLACEBO: secondary-data replications ##########\n")
sec_fit <- fit_elasticity(d80 %>% filter(repli_type == "secondary data"), "Secondary data")

## Original N is both a determinant of original precision (and hence markup)
## and a plausible direct anchor for replication scale. These models therefore
## estimate a distinct within-original-size association, not a robustness check
## for the portfolio estimand above.
fit_scale_conditioned <- function(dat, label) {
  m <- dat %>%
    filter(safeguard_defined, !is.na(N_naive), N_naive > 0,
           !is.na(orig_sample_size_value), orig_sample_size_value > 0,
           !is.na(repli_sample_size_value), repli_sample_size_value > 0,
           !is.na(discipline)) %>%
    mutate(logN_obs = log(repli_sample_size_value),
           logN_naive = log(N_naive), logmarkup = log(markup),
           logN_orig = log(orig_sample_size_value))
  discipline_counts <- table(m$discipline)
  m$discipline_collapsed <- ifelse(
    m$discipline %in% names(discipline_counts[discipline_counts >= 10]),
    as.character(m$discipline), "Other")
  linear <- lm(logN_obs ~ logN_naive + logmarkup + logN_orig, data = m)
  flexible <- lm(logN_obs ~ logN_naive + logmarkup +
                   ns(logN_orig, df = 3) + discipline_collapsed, data = m)
  cat("\n==", label, "conditioned on original sample size ==\n")
  print(coeftest(linear, vcov = vcovCL(linear, cluster = m$paper_id, type = "HC1")))
  cat("\nFlexible original-size function plus discipline:\n")
  print(coeftest(flexible, vcov = vcovCL(flexible, cluster = m$paper_id, type = "HC1")))
  invisible(list(linear = linear, flexible = flexible, data = m))
}

new_scale_fit <- fit_scale_conditioned(
  d80 %>% filter(repli_type == "new data"), "New data")
sec_scale_fit <- fit_scale_conditioned(
  d80 %>% filter(repli_type == "secondary data"), "Secondary data")

## ---- 3. Placebo tercile-gradient table (comparable to Table 4) --------

tercile_table <- function(dat, label) {
  m <- dat %>% filter(safeguard_defined, !is.na(ratio_obs_safeguard)) %>%
    mutate(frag_tercile = ntile(abs_t, 3))
  cat("\n==", label, "tercile gradient ==\n")
  print(as.data.frame(m %>% group_by(frag_tercile) %>% summarise(
    n = n(),
    t_range = paste0(round(min(abs_t),2), "-", round(max(abs_t),2)),
    median_markup = round(median(markup), 3),
    median_ratio = round(median(ratio_obs_safeguard), 3),
    share_below = round(mean(ratio_obs_safeguard < 1) * 100, 1)
  )))
}

tercile_table(d80 %>% filter(repli_type == "new data"), "New data")
tercile_table(d80 %>% filter(repli_type == "secondary data"), "Secondary data")

## ---- 4. Fallback (r-conversion) concentration by tercile ---------------

cat("\n########## Fallback concentration check ##########\n")
fb <- d80 %>% filter(repli_type == "new data", safeguard_defined) %>%
  mutate(frag_tercile = ntile(abs_t, 3))
print(table(fb$frag_tercile, fb$t_source))

fb_all <- d80 %>% filter(safeguard_defined) %>%
  mutate(frag_tercile = ntile(abs_t, 3))
cat("\nAll claims (new+secondary), fallback by tercile:\n")
print(table(fb_all$frag_tercile, fb_all$t_source))

## ---- 5. Selection check on the binary-success reversal -----------------

cat("\n########## Selection check: is the reversal driven by |t_orig|/N_naive? ##########\n")
sel <- d80 %>% filter(repli_type == "new data", safeguard_defined,
                       !is.na(ratio_obs_safeguard)) %>%
  mutate(log_ratio = log(ratio_obs_safeguard),
         success = as.integer(repli_score_criteria_met))

m1 <- glm(success ~ log_ratio, data = sel, family = binomial)
m2 <- glm(success ~ log_ratio + abs_t, data = sel, family = binomial)
m3 <- glm(success ~ log_ratio + log(N_naive), data = sel, family = binomial)
m4 <- glm(success ~ log_ratio + abs_t + log(N_naive), data = sel,
          family = binomial)

cat("\n-- Model 1: success ~ log(ratio) --\n"); print(summary(m1)$coefficients)
cat("\n-- Model 2: success ~ log(ratio) + t_orig --\n"); print(summary(m2)$coefficients)
cat("\n-- Model 3: success ~ log(ratio) + log(N_naive) --\n"); print(summary(m3)$coefficients)
cat("\n-- Model 4: full controls --\n"); print(summary(m4)$coefficients)

cat("\ncor(log_ratio, abs_t):", round(cor(sel$log_ratio, sel$abs_t), 3), "\n")
cat("cor(log_ratio, log(N_naive)):", round(cor(sel$log_ratio, log(sel$N_naive)), 3), "\n")

## ---- 6. Sensitivity grid: power x gamma ---------------------------------

cat("\n########## Sensitivity grid ##########\n")
grid <- expand.grid(power = c(.80, .90), gamma = c(.60, .80))
grid_results <- lapply(seq_len(nrow(grid)), function(i) {
  p <- grid$power[i]; g <- grid$gamma[i]
  dd <- build_frame(d, power = p, gamma = g)
  nd <- dd %>% filter(repli_type == "new data", safeguard_defined,
                       !is.na(ratio_obs_safeguard), !is.na(discipline))
  m <- nd %>% mutate(logN_obs = log(repli_sample_size_value),
                      logN_naive = log(N_naive), logmarkup = log(markup))
  fit <- tryCatch(lm(logN_obs ~ logN_naive + logmarkup + discipline, data = m),
                   error = function(e) NULL)
  b2 <- if (!is.null(fit)) coef(fit)["logmarkup"] else NA_real_
  data.frame(power = p, gamma = g, n = nrow(nd),
             median_markup = median(nd$markup, na.rm = TRUE),
             median_ratio = median(nd$ratio_obs_safeguard, na.rm = TRUE),
             share_below = mean(nd$ratio_obs_safeguard < 1, na.rm = TRUE) * 100,
             beta2 = unname(b2))
})
grid_df <- bind_rows(grid_results)
print(grid_df)

## ---- 7. Save ------------------------------------------------------------

write.csv(d80, "claim_level_with_discipline.csv", row.names = FALSE)
write.csv(grid_df, "sensitivity_grid.csv", row.names = FALSE)
cat("\nSaved claim_level_with_discipline.csv and sensitivity_grid.csv\n")
