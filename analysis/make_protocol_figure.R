#!/usr/bin/env Rscript

# Vector flowchart for the prospective planning protocol in the manuscript.
# Uses only base R's grid package so that it can be regenerated without
# additional graphics dependencies.

library(grid)

out <- file.path("..", "manuscript", "fig_planning_protocol.pdf")
example_file <- "protocol_example.csv"
if (!file.exists(example_file)) {
  stop("Run manuscript_tables.R first to create protocol_example.csv")
}
example <- read.csv(example_file)
stopifnot(nrow(example) == 1)

navy <- "#17365D"
blue_fill <- "#EAF1F8"
gold_fill <- "#FFF3D6"
green_fill <- "#E8F3EE"
grey_fill <- "#F1F2F3"
text_col <- "#20262E"

pdf(out, width = 9.0, height = 7.25, family = "Times", useDingbats = FALSE)
grid.newpage()

draw_box <- function(x, y, w, h, number, title, subtitle,
                     fill = blue_fill, title_size = 11.0,
                     subtitle_size = 9.4, top_row = FALSE) {
  grid.roundrect(
    x = unit(x, "npc"), y = unit(y, "npc"),
    width = unit(w, "npc"), height = unit(h, "npc"),
    r = unit(0.06, "snpc"),
    gp = gpar(fill = fill, col = navy, lwd = 1.15)
  )
  badge_x <- x - w / 2 + 0.034
  grid.circle(
    x = unit(badge_x, "npc"), y = unit(y + h / 2 - 0.034, "npc"),
    r = unit(0.022, "npc"),
    gp = gpar(fill = navy, col = navy)
  )
  grid.text(
    number, x = unit(badge_x, "npc"), y = unit(y + h / 2 - 0.034, "npc"),
    gp = gpar(col = "white", fontsize = 8.2, fontface = "bold")
  )
  grid.text(
    title, x = unit(x, "npc"),
    y = unit(y + if (top_row) -0.005 else 0.034, "npc"),
    gp = gpar(col = text_col, fontsize = title_size, fontface = "bold")
  )
  grid.text(
    subtitle, x = unit(x, "npc"),
    y = unit(y - if (top_row) 0.045 else 0.035, "npc"),
    gp = gpar(col = text_col, fontsize = subtitle_size, lineheight = 1.05)
  )
}

draw_arrow <- function(x0, y0, x1, y1) {
  grid.lines(
    x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
    arrow = arrow(type = "closed", length = unit(0.105, "inches")),
    gp = gpar(col = navy, lwd = 1.25)
  )
}

# Calculation and diagnosis
xs <- c(0.13, 0.375, 0.62, 0.865)
draw_box(xs[1], 0.875, 0.21, 0.16, "1", "Conventional target",
         expression(italic(N)[naive] ~ "from the reported effect"),
         subtitle_size = 8.9, top_row = TRUE)
draw_box(xs[2], 0.875, 0.21, 0.16, "2", "Original fragility",
         expression("Reconstruct" ~ italic(t)[o] ~ "and its provenance"),
         subtitle_size = 8.9, top_row = TRUE)
draw_box(xs[3], 0.875, 0.21, 0.16, "3", "Safeguard target",
         expression("Solve numerical" ~ italic(N)[safeguard] ~ "at stated" ~ gamma),
         subtitle_size = 8.6, top_row = TRUE)
draw_box(xs[4], 0.875, 0.21, 0.16, "4", "Price uncertainty",
         expression(italic(M) == italic(N)[safeguard] / italic(N)[naive]),
         subtitle_size = 8.9, top_row = TRUE)

for (i in 1:3) draw_arrow(xs[i] + 0.108, 0.875, xs[i + 1] - 0.108, 0.875)

# The substantive choice, preregistration, and interpretation
draw_box(
  0.50, 0.635, 0.86, 0.145, "5", "Choose the final sample under explicit constraints",
  "Budget and data  |  design fidelity  |  SESOI or prediction  |  evidential variety",
  fill = gold_fill, title_size = 11.5, subtitle_size = 9.5
)
draw_arrow(0.865, 0.792, 0.865, 0.742)
grid.lines(
  x = unit(c(0.865, 0.50), "npc"), y = unit(c(0.742, 0.742), "npc"),
  gp = gpar(col = navy, lwd = 1.25)
)
draw_arrow(0.50, 0.742, 0.50, 0.708)

draw_box(
  0.50, 0.445, 0.86, 0.125, "6", "Register the allocation decision",
  expression(italic(N)[naive] * ", " * italic(N)[safeguard] * ", planned " * italic(N) *
               ", ratio, and justification"),
  fill = green_fill, title_size = 11.3, subtitle_size = 9.5
)
draw_arrow(0.50, 0.560, 0.50, 0.510)

draw_box(
  0.50, 0.285, 0.86, 0.105, "7", "Keep planning separate from replication success",
  "Compliance with an uncertainty-aware rule is not an outcome criterion",
  fill = grey_fill, title_size = 11.3, subtitle_size = 9.4
)
draw_arrow(0.50, 0.380, 0.50, 0.340)

# A retrospective walk-through using one public SCORE claim. The missing
# allocation reason is itself the prospective protocol's motivating datum.
example_line_1 <- sprintf(
  "Reported r = %.3f, original N = %.0f, t = %.3f  ->  naive N = %.0f",
  example$orig_conv_r, example$orig_sample_size_value, example$t_orig,
  example$N_naive
)
example_line_2 <- sprintf(
  "Safeguarded r = %.3f  ->  safeguard N = %.0f, markup = %.2f",
  example$r_sg, example$N_safeguard, example$markup
)
example_line_3 <- sprintf(
  "Observed N = %.0f, observed/target = %.2f, shortfall = %.0f",
  example$N_observed, example$observed_to_safeguard, example$shortfall
)
grid.roundrect(
  x = unit(0.50, "npc"), y = unit(0.095, "npc"),
  width = unit(0.86, "npc"), height = unit(0.145, "npc"),
  r = unit(0.04, "snpc"),
  gp = gpar(fill = "#F7F4EC", col = navy, lwd = 1.15)
)
grid.text(
  "Worked SCORE example: WLpV_single-trace (retrospective)",
  x = unit(0.50, "npc"), y = unit(0.135, "npc"),
  gp = gpar(col = text_col, fontsize = 10.3, fontface = "bold")
)
grid.text(
  paste(example_line_1, example_line_2, example_line_3,
        "Allocation reason not recorded", sep = "\n"),
  x = unit(0.50, "npc"), y = unit(0.073, "npc"),
  gp = gpar(col = text_col, fontsize = 9.3, lineheight = 1.02)
)

dev.off()
message("Wrote ", normalizePath(out))
