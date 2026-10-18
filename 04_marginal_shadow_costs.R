# =============================================================================
# Marginal shadow costs                                   (Supplementary Code 4)
#
# Economic evaluation of site-specific weed management under imperfect sensing:
# a multi-year bio-economic simulation study.
#
# Figure of the marginal shadow cost of one additional weed plant (ET strategy,
# 1 x 1 m grid cells, all fields pooled) and a summary table by density bin.
#
# The marginal shadow cost is computed per grid cell by the simulation model
# (Supplementary Code 1, SHADOW_COST_MODE = 'fixed_management'): the change in
# the discounted ten-year weed liability caused by one additional weed added to
# the year-0 population before the herbicide application, with the cell's
# treatment decisions in all years held fixed. Marg_SC_D (dicots) and Marg_SC_M
# (monocots) are ten-year NPVs in EUR per weed. They are never negative; the
# share of values below zero is still reported as a check.
#
# INPUT
#   The per-cell CSV files written by the simulation model to
#   <results folder>/cells/ (ENABLE_CELL_EXPORT = True), one per field:
#   <project>_SCENARIO_1_ECON.csv. Attribute tables of the 1 x 1 m ET grid
#   layers exported from QGIS can be used as well (set FILE_PATTERN to
#   "\\.csv$"); they need the columns listed in 'needed' below.
#
# OUTPUT (in OUTPUT_FOLDER)
#   Fig_marginal_shadow_costs.png / .pdf and Marginal_shadow_costs_summary.csv
#
# REQUIREMENTS
#   R >= 4.2 with dplyr, ggplot2, purrr, readr and tidyr.
# =============================================================================

library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(tidyr)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

INPUT_FOLDER  <- "C:/SSWM_simulation/results/cells"   # per-cell files of the model
OUTPUT_FOLDER <- "C:/SSWM_simulation/results"
FILE_PATTERN  <- "_SCENARIO_1_ECON\\.csv$"             # ET runs (1 x 1 m by default)

EXPECTED_N_FILES <- 19           # set to NA to skip the check

SUMMARY_LINE  <- "upper_envelope" # "upper_envelope" (as reported), "median_smooth" or "none"
ENVELOPE_Q    <- 0.95             # upper_envelope: binned quantile the curve follows
                                  #   (0.95 = robust upper envelope; 1 = binned maxima)
SHOW_MEDIAN_TOO <- FALSE          # add the smoothed median as a dashed line
LOESS_SPAN    <- 0.6             # smoothness of the median line (0.4-0.9)
MIN_BIN_N     <- 20              # bins with fewer cells are not summarised
MIN_CELL_AREA <- 0.5             # m2; smaller edge cells are excluded (in a cell of a
                                 # few cm2, one extra weed is an unrealistic density jump)

# density range shown (plants m-2) and bin width per panel
max_x_spray   <- 60
max_x_no_d    <- 10
max_x_no_m    <- 5
bin_w_spray   <- 1.0
bin_w_no_d    <- 0.25
bin_w_no_m    <- 0.10

# y-axis: central share of the points kept in view (outliers still count in
# the summary line and table, they are only outside the plotted range)
Y_QUANTILES <- c(0.005, 0.995)

FIG_SCALE <- 1.0                 # raise to enlarge all text
sz <- function(x) x * FIG_SCALE


# =============================================================================
# 1. IMPORT
# =============================================================================

file_list <- list.files(INPUT_FOLDER, pattern = FILE_PATTERN, full.names = TRUE)
if (length(file_list) == 0) stop("No files matching ", FILE_PATTERN, " in: ", INPUT_FOLDER)
if (!is.na(EXPECTED_N_FILES) && length(file_list) != EXPECTED_N_FILES)
  warning(sprintf("Found %d CSV files, expected %d.", length(file_list), EXPECTED_N_FILES))

needed <- c("Spray_Dec", "P_Dicot", "P_Monocot", "TrueDenD", "TrueDenM",
            "Marg_SC_D", "Marg_SC_M", "Area_m2")

read_one <- function(f) {
  d <- read_csv(f, show_col_types = FALSE, progress = FALSE)
  miss <- setdiff(needed, names(d))
  if (length(miss) > 0) stop(basename(f), " is missing: ", paste(miss, collapse = ", "))
  d %>% select(all_of(needed)) %>%
    mutate(Field = sub("_SCENARIO_1_ECON\\.csv$|_shadow\\.csv$|\\.csv$", "", basename(f)))
}
df <- map_dfr(file_list, read_one)
message(sprintf("Loaded %d cells from %d fields.", nrow(df), n_distinct(df$Field)))

n_sliver <- sum(df$Area_m2 < MIN_CELL_AREA)
df <- df %>% filter(Area_m2 >= MIN_CELL_AREA)
message(sprintf("Excluded %d edge cells smaller than %.2f m2.", n_sliver, MIN_CELL_AREA))


# =============================================================================
# 2. LONG FORMAT, ONE ROW PER CELL AND WEED CLASS
# =============================================================================

panel_levels <- c("Dicots in sprayed cells", "Monocots in sprayed cells",
                  "Dicots in non-sprayed cells", "Monocots in non-sprayed cells")

long <- bind_rows(
  df %>% transmute(Field, Spray_Dec, Class = "Dicots",   Detected = P_Dicot,
                   Dens = TrueDenD, SC = Marg_SC_D),
  df %>% transmute(Field, Spray_Dec, Class = "Monocots", Detected = P_Monocot,
                   Dens = TrueDenM, SC = Marg_SC_M)
) %>%
  mutate(Sprayed = Spray_Dec == "Spray",
         Panel = paste(Class, ifelse(Sprayed, "in sprayed cells", "in non-sprayed cells")),
         Panel = factor(Panel, levels = panel_levels),
         x_max = case_when(Sprayed ~ max_x_spray,
                           Class == "Dicots" ~ max_x_no_d,
                           TRUE ~ max_x_no_m),
         bin_w = case_when(Sprayed ~ bin_w_spray,
                           Class == "Dicots" ~ bin_w_no_d,
                           TRUE ~ bin_w_no_m)) %>%
  # sprayed panels: the class must have been detected (spraying follows detections)
  filter(is.finite(SC), Dens > 0, Dens <= x_max, !Sprayed | Detected > 0)

counts <- long %>% group_by(Panel) %>%
  summarise(n = n(), share_negative = mean(SC < 0), .groups = "drop")
message("Cells per panel:")
for (i in seq_len(nrow(counts)))
  message(sprintf("  %-30s n = %6d   below zero: %4.1f%%", counts$Panel[i],
                  counts$n[i], 100 * counts$share_negative[i]))


# =============================================================================
# 3. BINNED SUMMARIES AND SUMMARY LINE
# =============================================================================

bins <- long %>%
  mutate(x_bin = (floor(Dens / bin_w) + 0.5) * bin_w) %>%
  group_by(Panel, x_bin) %>%
  summarise(n = n(),
            y_med = median(SC), y_q25 = quantile(SC, 0.25),
            y_q75 = quantile(SC, 0.75), y_max = max(SC),
            y_upq = quantile(SC, ENVELOPE_Q), .groups = "drop") %>%
  filter(n >= MIN_BIN_N)

smooth_panel <- function(b, span) {
  if (nrow(b) < 5) return(tibble(x_bin = b$x_bin, y_fit = b$y_med,
                                 lo = b$y_q25, hi = b$y_q75))
  grid_x <- seq(min(b$x_bin), max(b$x_bin), length.out = 300)
  fit <- function(y) {
    m <- tryCatch(loess(y ~ x_bin, data = transform(b, y = y), weights = sqrt(n),
                        span = span, degree = 1), error = function(e) NULL)
    if (is.null(m)) approx(b$x_bin, y, grid_x)$y else predict(m, data.frame(x_bin = grid_x))
  }
  tibble(x_bin = grid_x, y_fit = fit(b$y_med), lo = fit(b$y_q25), hi = fit(b$y_q75))
}

envelope_panel <- function(b, q_col) {
  # Smooth, monotone decreasing curve y = a + b1/(x+c) + b2/(x+c)^2 fitted to
  # the binned upper quantile, starting at the bin where that quantile peaks.
  # The offset c is chosen by grid search (least squares), so the curvature
  # adapts to each panel instead of being fixed.
  b <- b %>% mutate(y = .data[[q_col]])
  b <- b %>% filter(x_bin >= x_bin[which.max(y)])
  grid_x <- seq(min(b$x_bin), max(b$x_bin), length.out = 300)
  if (nrow(b) < 4) return(tibble(x_bin = grid_x, y_fit = approx(b$x_bin, b$y, grid_x)$y))
  design <- function(x, cc) data.frame(X1 = 1 / (x + cc), X2 = 1 / (x + cc)^2)
  best <- NULL; best_sse <- Inf
  for (cc in exp(seq(log(0.01), log(50), length.out = 120))) {
    d <- cbind(design(b$x_bin, cc), y = b$y)
    m <- tryCatch(lm(y ~ X1 + X2, data = d), error = function(e) NULL)
    if (is.null(m)) next
    sse <- sum(residuals(m)^2)
    if (sse < best_sse) { best_sse <- sse; best <- list(m = m, cc = cc) }
  }
  y_fit <- predict(best$m, newdata = design(grid_x, best$cc))
  y_fit <- cummin(pmax(y_fit, 0))          # enforce: never rises, never below zero
  tibble(x_bin = grid_x, y_fit = y_fit)
}

line_df <- NULL; med_df <- NULL
if (SUMMARY_LINE == "median_smooth" || SHOW_MEDIAN_TOO) {
  med_df <- bins %>% group_by(Panel) %>%
    group_modify(~ smooth_panel(.x, LOESS_SPAN)) %>% ungroup()
}
if (SUMMARY_LINE == "median_smooth") {
  line_df <- med_df
} else if (SUMMARY_LINE %in% c("upper_envelope", "envelope")) {
  q_col <- if (ENVELOPE_Q >= 1) "y_max" else "y_upq"
  line_df <- bins %>% group_by(Panel) %>%
    group_modify(~ envelope_panel(.x, q_col)) %>% ungroup()
}


# =============================================================================
# 4. FIGURE (one ggplot, four facets)
# =============================================================================

y_lims <- long %>% group_by(Panel) %>%
  summarise(lo = quantile(SC, Y_QUANTILES[1]), hi = quantile(SC, Y_QUANTILES[2]),
            x_max = first(x_max), .groups = "drop") %>%
  mutate(lo = pmin(lo, 0), hi = pmax(hi, 0))
# the axis must also reach the top of the summary line, and the points up to
# that height are then shown too (otherwise the line starts above the cloud)
if (!is.null(line_df)) {
  y_lims <- y_lims %>%
    left_join(line_df %>% group_by(Panel) %>% summarise(line_hi = max(y_fit, na.rm = TRUE)),
              by = "Panel") %>%
    mutate(hi = pmax(hi, coalesce(line_hi, hi))) %>% select(-line_hi)
}
y_lims <- y_lims %>% mutate(pad = 0.05 * (hi - lo), lo = lo - pad, hi = hi + pad)

# invisible points that set each facet's axis ranges
range_df <- bind_rows(
  y_lims %>% transmute(Panel, Dens = 0,     SC = lo),
  y_lims %>% transmute(Panel, Dens = x_max, SC = hi))

in_view <- long %>% left_join(y_lims %>% select(Panel, lo, hi), by = "Panel") %>%
  filter(SC >= lo, SC <= hi)

label_df <- counts %>%
  mutate(lab = ifelse(share_negative > 0.0005,
                      sprintf("n = %s\n%.1f%% of cells below zero",
                              format(n, big.mark = ","), 100 * share_negative),
                      sprintf("n = %s", format(n, big.mark = ","))))

p <- ggplot(in_view, aes(Dens, SC)) +
  geom_blank(data = range_df) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  geom_point(alpha = 0.10, size = 0.9, colour = "#555555", stroke = 0)

if (SUMMARY_LINE == "median_smooth") {
  p <- p +
    geom_ribbon(data = line_df, inherit.aes = FALSE,
                aes(x = x_bin, ymin = lo, ymax = hi), fill = "black", alpha = 0.15) +
    geom_line(data = line_df, inherit.aes = FALSE,
              aes(x = x_bin, y = y_fit), colour = "black", linewidth = 1.1)
} else if (SUMMARY_LINE %in% c("upper_envelope", "envelope")) {
  if (SHOW_MEDIAN_TOO)
    p <- p + geom_line(data = med_df, inherit.aes = FALSE,
                       aes(x = x_bin, y = y_fit), colour = "black",
                       linewidth = 0.8, linetype = "dashed")
  p <- p + geom_line(data = line_df, inherit.aes = FALSE,
                     aes(x = x_bin, y = y_fit), colour = "black", linewidth = 1.1)
}

env_txt <- if (ENVELOPE_Q >= 1) "binned maxima" else
  sprintf("the %s percentile of each density bin", paste0(round(100 * ENVELOPE_Q), "th"))
line_caption <- switch(SUMMARY_LINE,
  median_smooth  = "Line: smoothed median of density bins; ribbon: interquartile range.",
  upper_envelope = paste0("Line: upper envelope, a smooth decreasing curve fitted to ", env_txt,
                          if (SHOW_MEDIAN_TOO) "; dashed line: smoothed median." else "."),
  envelope       = paste0("Line: upper envelope fitted to ", env_txt, "."),
  none           = "")

p <- p +
  geom_text(data = label_df, inherit.aes = FALSE,
            aes(x = Inf, y = Inf, label = lab), hjust = 1.05, vjust = 1.3,
            size = sz(3.4), colour = "grey30", lineheight = 0.95) +
  facet_wrap(~ Panel, ncol = 2, scales = "free") +
  scale_x_continuous(n.breaks = 7, expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(n.breaks = 7, labels = function(v) format(v, scientific = FALSE,
                                                               drop0trailing = TRUE)) +
  labs(x = expression(bold("Estimated true pre-treatment density (plants "*m^-2*")")),
       y = "Marginal shadow cost (\u20ac per weed, 10-year NPV)",
       caption = paste(strwrap(paste(
                 "ET strategy, 1 x 1 m grid cells pooled across fields;",
                 "future management of each cell held fixed.", line_caption),
                 width = 120), collapse = "\n")) +
  theme_bw(base_size = sz(13)) +
  theme(strip.text       = element_text(face = "bold", size = sz(13)),
        strip.background = element_rect(fill = "grey90"),
        axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(hjust = 0, size = sz(10)))

dir.create(OUTPUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUTPUT_FOLDER, "Fig_marginal_shadow_costs.png"), p,
       width = 12, height = 10, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_FOLDER, "Fig_marginal_shadow_costs.pdf"), p,
       width = 12, height = 10, bg = "white", device = cairo_pdf)   # cairo: renders the euro sign
print(p)


# =============================================================================
# 5. SUMMARY TABLE FOR THE TEXT
# =============================================================================

summary_tab <- bins %>%
  left_join(counts, by = "Panel") %>%
  transmute(Panel, Density_bin_mid = round(x_bin, 3), Cells_in_bin = n.x,
            Median_EUR = signif(y_med, 3), Q25_EUR = signif(y_q25, 3),
            Q75_EUR = signif(y_q75, 3), Upper_Q_EUR = signif(y_upq, 3),
            Max_EUR = signif(y_max, 3),
            Panel_share_below_zero = round(share_negative, 3))
write_csv(summary_tab, file.path(OUTPUT_FOLDER, "Marginal_shadow_costs_summary.csv"))

message("\nMedian marginal shadow cost at the lowest and highest summarised density:")
summary_tab %>% group_by(Panel) %>%
  summarise(first = sprintf("%.5f at %.2f", Median_EUR[1], Density_bin_mid[1]),
            last  = sprintf("%.5f at %.2f", Median_EUR[n()], Density_bin_mid[n()]),
            .groups = "drop") %>%
  { for (i in seq_len(nrow(.))) message(sprintf("  %-30s %s  ->  %s",
                                                .$Panel[i], .$first[i], .$last[i])) }

message("\nFigure and summary written to: ", OUTPUT_FOLDER)
