# Figure 1: added-variable (partial-regression) plot of log(N_observed) on
# log(markup), controlling for log(N_naive), new data vs. secondary data.
# The slope in each panel equals beta_2 from Table 2 (Frisch-Waugh-Lovell:
# residualizing both axes on log(N_naive) recovers the partial slope from the
# full trivariate regression). The dashed slope-1 line is the safeguard
# criterion's prediction (N_observed = N_naive * markup, holding N_naive
# fixed implies a one-to-one line in this residualized space).
#
# Figure 2 (optional): concentration curve of the new-data resource shortfall
# across claims ordered by ascending original-evidence precision (most
# fragile first), against the 45-degree equal-share reference.
# Appendix Figure A1: signed approximation error against original decisiveness.

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

d <- read.csv("claim_level_with_discipline.csv")
exact_targets <- read.csv("safeguard_inflation_claim_level.csv") %>%
  select(claim_id, N_safeguard_exact)
d <- d %>% left_join(exact_targets, by = "claim_id") %>%
  mutate(markup_exact = N_safeguard_exact / N_naive)

nd <- d %>% filter(repli_type == "new data", safeguard_defined, !is.na(N_naive)) %>%
  mutate(logN_obs = log(repli_sample_size_value), logN_naive = log(N_naive),
         logmarkup = log(markup_exact))
sd_ <- d %>% filter(repli_type == "secondary data", safeguard_defined, !is.na(N_naive)) %>%
  mutate(logN_obs = log(repli_sample_size_value), logN_naive = log(N_naive),
         logmarkup = log(markup_exact))

partial_frame <- function(dat, label) {
  rY <- resid(lm(logN_obs ~ logN_naive, data = dat))
  rX <- resid(lm(logmarkup ~ logN_naive, data = dat))
  data.frame(sample = label, resid_markup = rX, resid_Nobs = rY)
}

pf <- bind_rows(
  partial_frame(nd, "New data (n = 123)"),
  partial_frame(sd_, "Secondary data (n = 123)")
)
pf$sample <- factor(pf$sample, levels = c("New data (n = 123)", "Secondary data (n = 123)"))

lim <- 4  # zoom for legibility; 6 secondary-data points fall outside the view

beta_labels <- data.frame(
  sample = factor(c("New data (n = 123)", "Secondary data (n = 123)"),
                  levels = levels(pf$sample)),
  resid_markup = c(-3.7, -3.7),
  resid_Nobs = c(3.55, 3.55),
  label = c("hat(beta)[2] == -0.40", "hat(beta)[2] == -2.80")
)

p1 <- ggplot(pf, aes(x = resid_markup, y = resid_Nobs)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_point(alpha = 0.55, size = 1.6, color = "#2c3e50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#c0392b", linewidth = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "#2980b9", linewidth = 0.8) +
  geom_text(data = beta_labels,
            aes(x = resid_markup, y = resid_Nobs, label = label),
            inherit.aes = FALSE, parse = TRUE, hjust = 0,
            color = "#2980b9", size = 4) +
  coord_fixed(ratio = 1, xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  facet_wrap(~ sample) +
  labs(
    x = expression(log(M[numeric]) ~ "|" ~ log(N[naive]) ~ "  (residualized)"),
    y = expression(log(N[observed]) ~ "|" ~ log(N[naive]) ~ "  (residualized)")
  ) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("../manuscript/fig_elasticity_scatter.pdf", p1, width = 8.5, height = 4.4)
cat("Saved fig_elasticity_scatter.pdf\n")

## ---- Figure 2 (optional): concentration curve --------------------------

nd2 <- d %>%
  filter(repli_type == "new data", safeguard_defined,
         !is.na(N_safeguard_exact)) %>%
  mutate(shortfall = pmax(N_safeguard_exact - repli_sample_size_value, 0)) %>%
  arrange(abs_t) %>%  # most fragile (lowest |t_orig|) first
  mutate(
    claim_rank = row_number(),
    cum_claims_pct = 100 * claim_rank / n(),
    cum_shortfall_pct = 100 * cumsum(shortfall) / sum(shortfall)
  )

# prepend the (0,0) origin point for a proper curve
curve_df <- bind_rows(
  data.frame(cum_claims_pct = 0, cum_shortfall_pct = 0),
  nd2 %>% select(cum_claims_pct, cum_shortfall_pct)
)

p2 <- ggplot(curve_df, aes(x = cum_claims_pct, y = cum_shortfall_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.5) +
  geom_line(color = "#c0392b", linewidth = 0.9) +
  annotate("point", x = 100 * 41 / 123, y = 90.2, color = "#2980b9", size = 2.2) +
  annotate("text", x = 100 * 41 / 123 + 3, y = 90.2 - 6,
           label = "most-fragile tercile\ncarries 90.2% of shortfall",
           hjust = 0, size = 3.3, color = "#2980b9") +
  scale_x_continuous(limits = c(0, 100), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0.01, 0)) +
  coord_fixed(ratio = 1) +
  labs(
    x = "Claims, cumulative % (ordered by ascending\noriginal-evidence precision, most fragile first)",
    y = "Resource shortfall, cumulative % of total"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("../manuscript/fig_shortfall_concentration.pdf", p2, width = 6, height = 6.3)
cat("Saved fig_shortfall_concentration.pdf\n")

## ---- Appendix Figure A1: approximation error ---------------------------

approx_df <- d %>%
  filter(safeguard_defined, !is.na(markup_exact), markup_exact > 0,
         !is.na(markup_approx), !is.na(abs_t)) %>%
  mutate(error_pct = 100 * (markup_approx - markup_exact) / markup_exact)
stopifnot(nrow(approx_df) == 246)
worst <- approx_df %>% slice_min(error_pct, n = 1, with_ties = FALSE)
pA1 <- ggplot(approx_df, aes(x = abs_t, y = error_pct)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.45) +
  geom_vline(xintercept = qnorm(.80), linetype = "dotted",
             color = "grey65", linewidth = 0.45) +
  geom_vline(xintercept = 1.96, linetype = "dashed",
             color = "#c0392b", linewidth = 0.55) +
  geom_point(alpha = .55, size = 1.45, color = "#2c3e50") +
  geom_point(data = worst, size = 2.4, color = "#c0392b") +
  annotate("text", x = 2.1, y = -14.8,
           label = "largest understatement: 16.23%",
           hjust = 0, size = 3.2, color = "#c0392b") +
  scale_x_log10(breaks = c(1, 1.96, 3, 5, 10, 30, 100),
                labels = c("1", "1.96", "3", "5", "10", "30", "100")) +
  labs(x = expression("Original decisiveness " * "|" * t[o] * "|" *
                        " (log scale)"),
       y = "Signed approximation error (%)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave("../manuscript/fig_approximation_error.pdf", pA1,
       width = 6.5, height = 3.7)
cat("Saved fig_approximation_error.pdf\n")

figure_slopes <- c(
  "new data" = unname(coef(lm(resid_Nobs ~ resid_markup, data = pf %>% filter(sample == "New data (n = 123)")))[2]),
  "secondary data" = unname(coef(lm(resid_Nobs ~ resid_markup, data = pf %>% filter(sample == "Secondary data (n = 123)")))[2])
)
canonical <- read.csv("table2_elasticity.csv") %>% filter(specification == "Baseline")
canonical_slopes <- setNames(canonical$beta2, canonical$sample)
stopifnot(all(abs(figure_slopes[names(canonical_slopes)] - canonical_slopes) < 1e-10))
cat("\nSlopes verified against Table 2 beta_2:\n")
print(round(figure_slopes, 3))
