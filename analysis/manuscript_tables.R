# Reproduce every numerical result reported in the submission manuscript.
# Run from analysis/: Rscript manuscript_tables.R

suppressMessages({
  library(dplyr)
  library(pwr)
  library(sandwich)
  library(lmtest)
  library(splines)
})

load("../g5sny-osfstorage-archive/extracted/analyst data.RData")

alpha <- .05
target_power <- .80
gamma <- .80
z_gamma <- qnorm(gamma)

ro <- repli_outcomes %>%
  filter(!is_covid, repli_version_of_record) %>%
  select(claim_id, paper_id, repli_type, repli_sample_size_value,
         repli_conv_r, repli_score_criteria_met)
oo <- orig_outcomes %>%
  semi_join(ro, by = "claim_id") %>%
  select(claim_id, orig_conv_r, orig_sample_size_value,
         orig_stat_type, orig_stat_value, orig_stat_dof_1)
disc <- paper_metadata %>%
  select(paper_id, discipline = COS_pub_category) %>%
  distinct(paper_id, .keep_all = TRUE)

n_for_r <- function(r, power = target_power) {
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
    error = function(e) NA_real_)
  if (is.na(out)) NA_real_ else ceiling(out)
}

d <- oo %>%
  inner_join(ro, by = "claim_id") %>%
  left_join(disc, by = "paper_id") %>%
  filter(!is.na(orig_conv_r), !is.na(repli_conv_r),
         !is.na(orig_sample_size_value)) %>%
  mutate(
    t_from_stat = case_when(
      orig_stat_type %in% c("t", "z") ~ orig_stat_value,
      orig_stat_type == "F" & orig_stat_dof_1 == 1 & orig_stat_value >= 0 ~ sqrt(pmax(orig_stat_value, 0)),
      orig_stat_type == "chi_squared" & orig_stat_dof_1 == 1 & orig_stat_value >= 0 ~ sqrt(pmax(orig_stat_value, 0)),
      TRUE ~ NA_real_),
    t_from_stat = sign(orig_conv_r) * abs(t_from_stat),
    t_from_r = orig_conv_r * sqrt(orig_sample_size_value - 2) /
      sqrt(1 - orig_conv_r^2),
    t_orig = coalesce(t_from_stat, t_from_r),
    abs_t = abs(t_orig),
    t_source = if_else(!is.na(t_from_stat), "Direct", "Fallback"),
    safeguard_defined = abs_t > z_gamma,
    markup_approx = if_else(safeguard_defined,
                            (1 - z_gamma / abs_t)^(-2), NA_real_),
    r_sg = if_else(safeguard_defined,
                   sign(orig_conv_r) * abs(orig_conv_r) *
                     (1 - z_gamma / abs_t), NA_real_))
d$N_naive <- vapply(d$orig_conv_r, n_for_r, numeric(1))
d$N_safeguard <- vapply(d$r_sg, n_for_r, numeric(1))
d <- d %>% mutate(
  markup_exact = N_safeguard / N_naive,
  approximation_error_pct = 100 * (markup_approx - markup_exact) / markup_exact,
  abs_approximation_error_pct = abs(approximation_error_pct),
  ratio = repli_sample_size_value / N_safeguard,
  shortfall = pmax(N_safeguard - repli_sample_size_value, 0),
  surplus = pmax(repli_sample_size_value - N_safeguard, 0))

# Shortfall and surplus must be the positive and negative parts of exactly the
# same gap, computed from the numerically solved integer safeguard target.
stopifnot(all.equal(d$shortfall - d$surplus,
                    d$N_safeguard - d$repli_sample_size_value,
                    check.attributes = FALSE),
          all(is.na(d$shortfall) | is.na(d$surplus) |
                d$shortfall == 0 | d$surplus == 0))

valid <- d %>%
  filter(safeguard_defined, !is.na(N_safeguard), N_safeguard > 0,
         !is.na(repli_sample_size_value), repli_sample_size_value > 0) %>%
  group_by(repli_type) %>%
  mutate(tercile = ntile(abs_t, 3)) %>%
  ungroup()

# Body Table 1: tercile gradient, both designs.
gradient <- valid %>%
  group_by(tercile, repli_type) %>%
  summarise(n = n(), min_t = min(abs_t), max_t = max(abs_t),
            median_markup_exact = median(markup_exact, na.rm = TRUE),
            median_markup_approx = median(markup_approx),
            median_abs_approx_error_pct = median(abs_approximation_error_pct,
                                                  na.rm = TRUE),
            max_abs_approx_error_pct = max(abs_approximation_error_pct,
                                           na.rm = TRUE),
            median_ratio = median(ratio),
            pct_below = 100 * mean(ratio < 1), .groups = "drop")
write.csv(gradient, "table1_gradient.csv", row.names = FALSE)

# Body Table 2: elasticity specifications using the numerical target ratio.
regdata <- valid %>% filter(!is.na(N_naive), N_naive > 0) %>%
  mutate(logN = log(repli_sample_size_value),
         logNnaive = log(N_naive), logmarkup = log(markup_exact),
         logNorig = log(orig_sample_size_value))
collapse_disc <- function(x) {
  counts <- table(x)
  ifelse(x %in% names(counts[counts >= 10]), as.character(x), "Other")
}

cluster_extract <- function(fit, cluster, spec, sample_name) {
  V <- vcovCL(fit, cluster = cluster, type = "HC1")
  ct <- coeftest(fit, vcov. = V)
  b <- coef(fit)["logmarkup"]
  se <- sqrt(V["logmarkup", "logmarkup"])
  G <- length(unique(cluster))
  df <- G - 1
  data.frame(sample = sample_name, specification = spec,
             n = nobs(fit), papers = G,
             beta1 = coef(fit)["logNnaive"], beta2 = b, se_beta2 = se,
             ci95_low_beta2 = b - qt(.975, df = df) * se,
             ci95_high_beta2 = b + qt(.975, df = df) * se,
             r2 = summary(fit)$r.squared,
             t_beta2_eq_1 = (b - 1) / se,
             p_beta2_eq_1 = 2 * pt(-abs((b - 1) / se), df = df))
}

fits <- list()
for (typ in c("new data", "secondary data")) {
  x <- regdata %>% filter(repli_type == typ)
  fits[[length(fits) + 1]] <- cluster_extract(
    lm(logN ~ logNnaive + logmarkup, data = x), x$paper_id,
    "Baseline", typ)
  x$discipline_collapsed <- collapse_disc(x$discipline)
  fits[[length(fits) + 1]] <- cluster_extract(
    lm(logN ~ logNnaive + logmarkup + discipline_collapsed, data = x),
    x$paper_id, "Discipline fixed effects", typ)
  fits[[length(fits) + 1]] <- cluster_extract(
    lm(logN ~ logNnaive + logmarkup + logNorig, data = x), x$paper_id,
    "Original sample-size control", typ)
  fits[[length(fits) + 1]] <- cluster_extract(
    lm(logN ~ logNnaive + logmarkup + ns(logNorig, df = 3) +
         discipline_collapsed, data = x), x$paper_id,
    "Flexible original size + discipline", typ)
  xd <- x %>% filter(t_source == "Direct")
  fits[[length(fits) + 1]] <- cluster_extract(
    lm(logN ~ logNnaive + logmarkup, data = xd), xd$paper_id,
    "Direct-statistic subset", typ)
  xp <- x %>% group_by(paper_id) %>%
    summarise(logN = median(logN), logNnaive = median(logNnaive),
              logmarkup = median(logmarkup), .groups = "drop")
  fits[[length(fits) + 1]] <- cluster_extract(
    lm(logN ~ logNnaive + logmarkup, data = xp), xp$paper_id,
    "Paper-level medians", typ)
  xpo <- x %>% group_by(paper_id) %>%
    summarise(logN = median(logN), logNnaive = median(logNnaive),
              logmarkup = median(logmarkup), logNorig = median(logNorig),
              .groups = "drop")
  fits[[length(fits) + 1]] <- cluster_extract(
    lm(logN ~ logNnaive + logmarkup + logNorig, data = xpo), xpo$paper_id,
    "Paper medians + original size", typ)
}
elasticity <- bind_rows(fits)
write.csv(elasticity, "table2_elasticity.csv", row.names = FALSE)

# A transparent practical-equivalence diagnostic for the original-size-
# conditioned new-data coefficient. TOST requires both one-sided tests to
# reject; failure to reject beta = 1 alone is not evidence of equivalence.
equivalence_band <- c(.80, 1.25)
conditioned <- elasticity %>%
  filter(sample == "new data", specification == "Original sample-size control")
stopifnot(nrow(conditioned) == 1)
df_equiv <- conditioned$papers - 1
t_lower <- (conditioned$beta2 - equivalence_band[1]) / conditioned$se_beta2
t_upper <- (conditioned$beta2 - equivalence_band[2]) / conditioned$se_beta2
equivalence <- data.frame(
  sample = conditioned$sample,
  specification = conditioned$specification,
  lower_bound = equivalence_band[1], upper_bound = equivalence_band[2],
  beta2 = conditioned$beta2, se_beta2 = conditioned$se_beta2,
  df = df_equiv,
  ci90_low = conditioned$beta2 - qt(.95, df_equiv) * conditioned$se_beta2,
  ci90_high = conditioned$beta2 + qt(.95, df_equiv) * conditioned$se_beta2,
  p_lower = pt(t_lower, df_equiv, lower.tail = FALSE),
  p_upper = pt(t_upper, df_equiv),
  p_tost = max(pt(t_lower, df_equiv, lower.tail = FALSE),
               pt(t_upper, df_equiv)))
write.csv(equivalence, "appendix_equivalence.csv", row.names = FALSE)

# Body Table 3: numerical integer resource gaps and illustrative monetization.
resources <- valid %>% filter(repli_type == "new data") %>%
  group_by(tercile) %>%
  summarise(n = n(), claims_below = sum(shortfall > 0),
            missing_observations = sum(shortfall),
            surplus_observations = sum(surplus), .groups = "drop") %>%
  mutate(cost_5 = 5 * missing_observations,
         cost_10 = 10 * missing_observations,
         cost_20 = 20 * missing_observations,
         cost_50 = 50 * missing_observations)
write.csv(resources, "table3_resources.csv", row.names = FALSE)

# Frictionless reallocation benchmark. This is an accounting upper bound, not
# a claim that observations or marginal costs are fungible across studies.
nd_resources <- valid %>% filter(repli_type == "new data")
reallocation <- data.frame(
  claims = nrow(nd_resources),
  claims_meeting_observed = sum(nd_resources$shortfall == 0),
  pct_meeting_observed = 100 * mean(nd_resources$shortfall == 0),
  total_shortfall = sum(nd_resources$shortfall),
  total_surplus = sum(nd_resources$surplus),
  least_fragile_surplus = sum(nd_resources$surplus[nd_resources$tercile == 3]),
  shortfall_share_total_surplus = 100 * sum(nd_resources$shortfall) /
    sum(nd_resources$surplus),
  shortfall_share_least_fragile_surplus = 100 * sum(nd_resources$shortfall) /
    sum(nd_resources$surplus[nd_resources$tercile == 3]),
  claims_meeting_frictionless = if_else(
    sum(nd_resources$surplus) >= sum(nd_resources$shortfall),
    nrow(nd_resources), NA_integer_),
  pct_meeting_frictionless = if_else(
    sum(nd_resources$surplus) >= sum(nd_resources$shortfall), 100, NA_real_),
  remaining_surplus = sum(nd_resources$surplus) - sum(nd_resources$shortfall)
)
write.csv(reallocation, "reallocation_benchmark.csv", row.names = FALSE)

# Exact omitted-variable/FWL decomposition of the baseline-to-original-size
# reversal: beta_restricted = beta_full + gamma_Norig * pi_markup.
fwl_decompose <- function(x, sample_label) {
  restricted <- lm(logN ~ logNnaive + logmarkup, data = x)
  full <- lm(logN ~ logNnaive + logmarkup + logNorig, data = x)
  auxiliary <- lm(logNorig ~ logNnaive + logmarkup, data = x)
  rx <- resid(lm(logmarkup ~ logNnaive, data = x))
  rw <- resid(lm(logNorig ~ logNnaive, data = x))
  gamma <- coef(full)["logNorig"]
  pi_markup <- coef(auxiliary)["logmarkup"]
  data.frame(
    sample = sample_label, n = nrow(x), papers = n_distinct(x$paper_id),
    beta_restricted = coef(restricted)["logmarkup"],
    beta_full = coef(full)["logmarkup"], gamma_original_n = gamma,
    pi_markup = pi_markup, omitted_scale_component = gamma * pi_markup,
    reconstructed_restricted = coef(full)["logmarkup"] + gamma * pi_markup,
    partial_cor_markup_original_n = cor(rx, rw)
  )
}
new_reg <- regdata %>% filter(repli_type == "new data")
fwl <- bind_rows(
  fwl_decompose(new_reg, "All new-data claims"),
  fwl_decompose(new_reg %>% filter(t_source == "Direct"),
                "Direct-statistic subset")
)
stopifnot(max(abs(fwl$beta_restricted - fwl$reconstructed_restricted)) < 1e-10)
write.csv(fwl, "appendix_fwl_decomposition.csv", row.names = FALSE)

# Exploratory discipline heterogeneity. Small categories are pooled exactly as
# in the fixed-effect specification. Slopes are linear combinations from one
# interacted model, so log(N_naive) (and, where applicable, log original N)
# have common coefficients across disciplines.
discipline_slopes <- function(x, add_original_n = FALSE) {
  x$discipline_collapsed <- factor(collapse_disc(x$discipline))
  x$discipline_collapsed <- relevel(x$discipline_collapsed, ref = "business")
  form <- if (add_original_n) {
    logN ~ logNnaive + logNorig + logmarkup * discipline_collapsed
  } else {
    logN ~ logNnaive + logmarkup * discipline_collapsed
  }
  fit <- lm(form, data = x)
  V <- vcovCL(fit, cluster = x$paper_id, type = "HC1")
  interaction_terms <- grep("^logmarkup:discipline_collapsed", names(coef(fit)),
                            value = TRUE)
  b_int <- coef(fit)[interaction_terms]
  V_int <- V[interaction_terms, interaction_terms, drop = FALSE]
  q <- length(interaction_terms)
  wald_f <- as.numeric(t(b_int) %*% solve(V_int, b_int) / q)
  joint_p <- pf(wald_f, df1 = q, df2 = n_distinct(x$paper_id) - 1,
                lower.tail = FALSE)
  levs <- levels(x$discipline_collapsed)
  bind_rows(lapply(levs, function(g) {
    L <- setNames(rep(0, length(coef(fit))), names(coef(fit)))
    L["logmarkup"] <- 1
    int_name <- paste0("logmarkup:discipline_collapsed", g)
    if (int_name %in% names(L)) L[int_name] <- 1
    data.frame(
      specification = if_else(add_original_n, "Original-size conditioned",
                              "Portfolio"),
      discipline = g, n = sum(x$discipline_collapsed == g),
      papers = n_distinct(x$paper_id[x$discipline_collapsed == g]),
      beta2 = sum(L * coef(fit)),
      se = sqrt(as.numeric(t(L) %*% V %*% L)),
      joint_interaction_p = joint_p
    )
  }))
}
discipline_heterogeneity <- bind_rows(
  discipline_slopes(new_reg, FALSE), discipline_slopes(new_reg, TRUE)
)
write.csv(discipline_heterogeneity, "appendix_discipline_heterogeneity.csv",
          row.names = FALSE)

# Winner's-curse stress test. If the reported effect is inflated by q, the
# scenario sets the latent effect and t statistic to 1/(1+q) of their reported
# values, corresponding to a fixed reported standard error. This deliberately
# simple exercise is not an estimated selection correction.
winner_curse <- bind_rows(lapply(c(0, .10, .20), function(q) {
  correction <- 1 / (1 + q)
  x <- d %>% filter(repli_type == "new data") %>% mutate(
    r_adjusted = correction * orig_conv_r,
    t_adjusted = correction * t_orig,
    defined_wc = abs(t_adjusted) > z_gamma,
    r_sg_wc = if_else(defined_wc,
      sign(r_adjusted) * abs(r_adjusted) *
        (1 - z_gamma / abs(t_adjusted)), NA_real_)
  )
  x$N_naive_wc <- vapply(x$r_adjusted, n_for_r, numeric(1))
  x$N_safeguard_wc <- vapply(x$r_sg_wc, n_for_r, numeric(1))
  x <- x %>% mutate(
    markup_wc = N_safeguard_wc / N_naive_wc,
    shortfall_wc = pmax(N_safeguard_wc - repli_sample_size_value, 0)
  )
  data.frame(
    assumed_inflation_pct = 100 * q, correction_factor = correction,
    defined_claims = sum(x$defined_wc),
    undefined_claims = sum(!x$defined_wc),
    median_markup = median(x$markup_wc, na.rm = TRUE),
    median_safeguard_target = median(x$N_safeguard_wc, na.rm = TRUE),
    claims_below = sum(x$shortfall_wc > 0, na.rm = TRUE),
    pct_below = 100 * mean(x$shortfall_wc > 0, na.rm = TRUE),
    total_shortfall = sum(x$shortfall_wc, na.rm = TRUE)
  )
}))
write.csv(winner_curse, "appendix_winners_curse.csv", row.names = FALSE)

# Public SCORE claim used as the worked example in Figure 3.
protocol_example <- valid %>%
  filter(claim_id == "WLpV_single-trace") %>%
  transmute(claim_id, orig_conv_r, orig_sample_size_value, t_orig, r_sg,
            N_naive, N_safeguard, markup = markup_exact,
            N_observed = repli_sample_size_value,
            observed_to_safeguard = ratio, shortfall)
stopifnot(nrow(protocol_example) == 1)
write.csv(protocol_example, "protocol_example.csv", row.names = FALSE)

# Appendix A2: fallback source by tercile for all claims.
fallback <- valid %>% group_by(tercile, t_source) %>%
  summarise(n = n(), .groups = "drop")
write.csv(fallback, "appendix_fallback.csv", row.names = FALSE)

# Appendix A3: descriptive statistics, both designs.
desc <- valid %>% group_by(repli_type) %>%
  summarise(n = n(), papers = n_distinct(paper_id),
            median_abs_t = median(abs_t), iqr_abs_t = IQR(abs_t),
            median_naive = median(N_naive, na.rm = TRUE),
            iqr_naive = IQR(N_naive, na.rm = TRUE),
            median_observed = median(repli_sample_size_value),
            iqr_observed = IQR(repli_sample_size_value),
            median_markup = median(markup_exact, na.rm = TRUE),
            iqr_markup = IQR(markup_exact, na.rm = TRUE),
            median_ratio = median(ratio), iqr_ratio = IQR(ratio),
            pct_below = 100 * mean(ratio < 1), .groups = "drop")
write.csv(desc, "appendix_descriptives.csv", row.names = FALSE)

# Appendix A4: sensitivity, using the same clustered baseline specification.
sens <- list()
for (power in c(.80, .90)) for (g in c(.60, .80)) {
  zg <- qnorm(g)
  dg <- d %>% mutate(
    defined_g = abs_t > zg,
    markup_approx_g = if_else(defined_g, (1 - zg / abs_t)^(-2), NA_real_),
    r_sg_g = if_else(defined_g, sign(orig_conv_r) * abs(orig_conv_r) *
                       (1 - zg / abs_t), NA_real_))
  dg$N_naive_g <- vapply(dg$orig_conv_r, n_for_r, numeric(1), power = power)
  dg$N_sg_g <- vapply(dg$r_sg_g, n_for_r, numeric(1), power = power)
  for (typ in c("new data", "secondary data")) {
    x <- dg %>% filter(repli_type == typ, defined_g, !is.na(N_sg_g),
                       !is.na(N_naive_g), repli_sample_size_value > 0) %>%
      mutate(logN = log(repli_sample_size_value),
             markup_g = N_sg_g / N_naive_g,
             logNnaive = log(N_naive_g), logmarkup = log(markup_g),
             ratio_g = repli_sample_size_value / N_sg_g)
    fit <- lm(logN ~ logNnaive + logmarkup, data = x)
    V <- vcovCL(fit, cluster = x$paper_id, type = "HC1")
    sens[[length(sens) + 1]] <- data.frame(
      sample = typ, power = power, gamma = g, n = nrow(x),
      median_markup = median(x$markup_g), median_ratio = median(x$ratio_g),
      pct_below = 100 * mean(x$ratio_g < 1),
      beta2 = coef(fit)["logmarkup"], se = sqrt(V["logmarkup", "logmarkup"]))
  }
}
write.csv(bind_rows(sens), "appendix_sensitivity.csv", row.names = FALSE)

# Appendix A5: success logits, cluster-robust by paper.
sel <- valid %>% filter(repli_type == "new data") %>%
  mutate(success = as.integer(repli_score_criteria_met), logratio = log(ratio),
         logNnaive = log(N_naive))
logit_specs <- list("Ratio only" = success ~ logratio,
                    "Plus original t" = success ~ logratio + abs_t,
                    "Plus naive target" = success ~ logratio + logNnaive,
                    "Full controls" = success ~ logratio + abs_t + logNnaive)
logits <- lapply(names(logit_specs), function(nm) {
  fit <- glm(logit_specs[[nm]], data = sel, family = binomial)
  used_rows <- as.integer(rownames(model.frame(fit)))
  V <- vcovCL(fit, cluster = sel$paper_id[used_rows], type = "HC1")
  ct <- coeftest(fit, vcov. = V)
  data.frame(model = nm, n = nobs(fit), term = rownames(ct), estimate = ct[, 1],
             se = ct[, 2], p = ct[, 4], row.names = NULL)
}) %>% bind_rows()
write.csv(logits, "appendix_logit.csv", row.names = FALSE)

direct_t <- d %>% filter(orig_stat_type == "t")
cat("analytic claims:", nrow(d), "\n")
cat("direct/fallback:", table(d$t_source), "\n")
cat("fallback diagnostic n/median abs diff/Pearson/Spearman:",
    nrow(direct_t), median(abs(direct_t$t_from_r - direct_t$orig_stat_value)),
    cor(abs(direct_t$t_from_r), abs(direct_t$orig_stat_value)),
    cor(abs(direct_t$t_from_r), abs(direct_t$orig_stat_value), method = "spearman"), "\n")
cat("success overall:", sum(ro$repli_score_criteria_met), nrow(ro),
    mean(ro$repli_score_criteria_met), "\n")
