# =============================================================================
# Main result analysis                                    (Supplementary Code 2)
#
# Economic evaluation of site-specific weed management under imperfect sensing:
# a multi-year bio-economic simulation study.
#
# Reads the CSV files written by the simulation model (Supplementary Code 1),
# adds the technology and information cost scenarios (Table 3) and produces the
# tables and figures of the Results section and the Supplementary Material:
# herbicide savings and the linear mixed-effects model, the partial budget
# decomposition, break-even technology cost and break-even recall, the
# sensitivity analyses, and the efficacy of the reduced background dose.
#
# INPUT (in csv_folder; written by the simulation model)
#   results_thesis.csv                 one row per field x resolution x strategy
#   results_thesis_recall_curve.csv    break-even recall analysis
#   results_thesis_dd_kill_sweep.csv   efficacy of the background dose
#
# OUTPUT (in csv_folder)
#   CSV tables, the individual figures (PDF/PNG) and Results_Analysis.pdf with
#   all figures.
#
# REQUIREMENTS
#   R >= 4.2 with tidyverse, lme4, lmerTest, emmeans, multcomp, multcompView
#   (letters) and pbkrtest (Kenward-Roger degrees of freedom).
#   Tested with R 4.3.1.
# =============================================================================
library(tidyverse)
library(ggplot2)

# mixed-model packages. lmerTest is optional (it only adds
# Satterthwaite p-values); the script works with plain lme4 if it is absent.
have_lme4    <- requireNamespace("lme4",    quietly = TRUE)
have_emmeans <- requireNamespace("emmeans", quietly = TRUE)
have_multcomp<- requireNamespace("multcomp",quietly = TRUE)
have_lmerTest<- requireNamespace("lmerTest",quietly = TRUE)
if (have_lme4)     library(lme4)
if (have_lmerTest) library(lmerTest)
if (have_emmeans)  library(emmeans)


# =============================================================================
# USER CONFIGURATION
# =============================================================================

csv_folder <- "C:/SSWM_simulation/results"   # folder with the model output

csv_files <- c(
  "results_thesis.csv"
)

# Further files written by the simulation model.
recall_sweep_files <- c(          # validation design (spray map fixed); not used
                                  # when the recall curve file is present
  "results_thesis_recall_sweep.csv"
)
recall_curve_files <- c(          # recall curve (population fixed, sensor varied)
  "results_thesis_recall_curve.csv"
)
dd_sweep_files <- c(
  "results_thesis_dd_kill_sweep.csv"
)

# Recall values used for the +/-20% recall sensitivity. They must be present in
# the curve file; the nearest available values are used otherwise.
RECALL_SENS_LOW  <- 0.328
RECALL_SENS_HIGH <- 0.492

# Field-level band drawn around the mean recall curve: "iqr" (25th-75th
# percentile across fields), "sd" (mean +/- 1 SD) or "none".
CURVE_BAND <- "iqr"
# Upper end of the recall axis in the curve figures. Recall 1.0 is a limiting
# case (no escapes at all) where the curve jumps, so the axis ends at 0.99.
CURVE_X_MAX <- 0.99
CURVE_X_MIN <- 0.41   # lower end of the axis (the measured recall)
BE_RECALL_MAX <- 0.99 # highest recall searched for a break-even
# strategies drawn in Fig. 3 (break-even vs recall)
CURVE_STRATEGIES <- c("ET strategy", "PB strategy", "DD strategy")
# draw the field-level recall figure (strategies x resolutions)?
MAKE_FIELD_CURVE_FIGURE <- FALSE
CURVE_X_BREAKS <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99)

# Expected design, used by the completeness check below. 19 fields x
# 3 resolutions x 4 strategies = 228 rows in the main file.
EXPECTED_N_FIELDS <- 19

pdf_output <- file.path(csv_folder, "Results_Analysis.pdf")

# Technology cost scenario used for the main-text decomposition figure.
# Scenario 2 is the least expensive.
main_text_cost_scenario <- "Costs scenario 2"

# Background-dose efficacy at which the DD strategy is reported. The model
# simulates DD at DD_BASE_KILL = 0.85 (full dose 0.95, background-dose years
# charged 50% of the herbicide cost), so the main results file already holds
# the reported specification. The script compares DD_REPORT_K with the
# DD_Base_Kill recorded in the data:
#   equal        -> verify: the DD rows must reproduce the matching row of the
#                   background-efficacy sweep (the script stops otherwise).
#   different    -> substitute the four decomposition columns from the sweep
#                   (the sensitivity columns then belong to another case).
#   NA           -> use the simulated values without a check.
DD_REPORT_K <- 0.85


# =============================================================================
# TECHNOLOGY COST SCENARIOS (EUR ha-1, annualised)
# Applied once in R. Python output contains zero technology cost.
# =============================================================================

# single source of truth for resolution ordering
RES_LEVELS <- c("1x1m", "50x50cm", "25x25cm")

tech_cost_table <- tribble(
  ~Tech_Cost_Label,   ~Resolution,  ~Tech_Cost_EUR,
  "Costs scenario 1", "1x1m",        50,
  "Costs scenario 1", "50x50cm",     50,
  "Costs scenario 1", "25x25cm",     57,
  "Costs scenario 2", "1x1m",        29,
  "Costs scenario 2", "50x50cm",     29,
  "Costs scenario 2", "25x25cm",     32,
  "Costs scenario 3", "1x1m",       130,
  "Costs scenario 3", "50x50cm",    130,
  "Costs scenario 3", "25x25cm",    150,
  "Costs scenario 4", "1x1m",        50,
  "Costs scenario 4", "50x50cm",     50,
  "Costs scenario 4", "25x25cm",     55,
  "Costs scenario 5", "1x1m",       151,
  "Costs scenario 5", "50x50cm",    151,
  "Costs scenario 5", "25x25cm",    158
) %>%
  # Resolution must be a FACTOR with the correct level order. tribble()
  # creates it as character; left_join(factor, character) coerces the result to
  # character and ggplot then sorts alphabetically, putting 25x25cm before
  # 50x50cm on every axis.
  mutate(Resolution = factor(Resolution, levels = RES_LEVELS))

# Forward-looking ET (SCENARIO_4_DYNET). FALSE: excluded from the figures and
# tables; it is compared with the ET in SuppTable_ET_vs_DynET.csv (Section 3.6).
# TRUE: treated as a fourth strategy throughout.
INCLUDE_DYNET <- FALSE

# Greyscale fills, distinguishable in print.
pub_colors <- c(
  "ET strategy"    = "#f0f0f0",
  "PB strategy"    = "#bdbdbd",
  "DD strategy"    = "#252525",
  "DynET strategy" = "#737373"
)

# greyscale fills for the decomposition components, print-safe.
component_colors <- c(
  "Herbicide saved"     = "#e8e8e8",
  "Year-0 yield loss"   = "#a8a8a8",
  "Future weed penalty" = "#5c5c5c",
  "Technology cost"     = "#1a1a1a"
)


# =============================================================================
# DATA LOADING AND CLEANING
# =============================================================================

TEXT_COLUMNS <- c("Project", "Scenario", "Future_Detection", "ET_Future_Rule",
                  "DD_Future_Cost", "Survivor_Cond", "Alloc_Rule")

read_csv_file <- function(filename) {
  path <- file.path(csv_folder, filename)
  if (!file.exists(path)) { warning(paste("File not found:", path)); return(NULL) }
  # Configuration columns are text, all others numeric. Only columns that
  # exist in the file are declared.
  hdr <- names(read_csv(path, n_max = 0, col_types = cols(.default = "c"),
                        show_col_types = FALSE))
  txt <- intersect(TEXT_COLUMNS, hdr)
  ct  <- do.call(cols, c(list(.default = col_double()),
                         setNames(rep(list(col_character()), length(txt)), txt)))
  read_csv(path, col_types = ct) %>%
    mutate(Source_File = filename)
}

# Shared labelling so the main file and the sweep files are treated identically.
add_labels <- function(df) {
  df %>%
    mutate(
      Resolution = case_when(
        str_detect(Project, "-25$") ~ "25x25cm",
        str_detect(Project, "-50$") ~ "50x50cm",
        TRUE                        ~ "1x1m"
      ),
      Resolution = factor(Resolution, levels = RES_LEVELS),
      Field_Name = str_remove(Project, "-25$|-50$") %>%
        str_remove("[_-]econ[_-](model|evaluation)$"),
      Scenario_Label = case_when(
        Scenario == "SCENARIO_1_ECON"   ~ "ET strategy",
        Scenario == "SCENARIO_2_SPOT"   ~ "PB strategy",
        Scenario == "SCENARIO_3_SECURE" ~ "DD strategy",
        Scenario == "SCENARIO_4_DYNET"  ~ "DynET strategy",
        TRUE ~ Scenario
      ),
      Scenario_Label = factor(
        Scenario_Label,
        levels = c("ET strategy", "PB strategy", "DD strategy", "DynET strategy")
      )
    ) %>%
    droplevels()
}

message("Reading simulation output files...")
raw_data <- map_df(csv_files, read_csv_file)
if (nrow(raw_data) == 0) stop("No data loaded. Check csv_folder and csv_files.")
message(sprintf("Loaded %d rows from %d file(s).", nrow(raw_data), length(csv_files)))

clean_data <- add_labels(raw_data)

# Report what was found, then handle DynET according to the toggle.
message("Strategies present: ",
        paste(levels(droplevels(clean_data$Scenario_Label)), collapse = ", "))

has_dynet <- "DynET strategy" %in% levels(droplevels(clean_data$Scenario_Label))

if (has_dynet) {
  # ET vs forward-looking ET (Section 3.6), written regardless of the toggle:
  # does accounting for the future weed penalty change the decision?
  dynet_cmp <- clean_data %>%
    filter(Scenario_Label %in% c("ET strategy", "DynET strategy")) %>%
    select(Field_Name, Resolution, Scenario_Label,
           Sprayed_Pct, Vol_Red_Pct, SSWM_Liab_Base,
           any_of(c("D_Net", "Tech_Breakeven"))) %>%
    pivot_wider(names_from = Scenario_Label,
                values_from = -c(Field_Name, Resolution),
                names_sep = "_") %>%
    mutate(
      Diff_Sprayed = `Sprayed_Pct_DynET strategy` - `Sprayed_Pct_ET strategy`,
      Diff_Liab    = `SSWM_Liab_Base_DynET strategy` - `SSWM_Liab_Base_ET strategy`
    )
  write_csv(dynet_cmp, file.path(csv_folder, "SuppTable_ET_vs_DynET.csv"))
  
  n_identical <- sum(abs(dynet_cmp$Diff_Sprayed) < 0.01 &
                       abs(dynet_cmp$Diff_Liab) < 0.01, na.rm = TRUE)
  message(sprintf(
    "ET vs DynET: identical in %d of %d field x resolution combinations.",
    n_identical, nrow(dynet_cmp)))
  if (n_identical == nrow(dynet_cmp))
    message("  -> DynET never changed a decision. Worth one sentence in the ",
            "discussion; no need to add a strategy to every figure.")
  
  if (!INCLUDE_DYNET) {
    clean_data <- clean_data %>%
      filter(Scenario_Label != "DynET strategy") %>%
      droplevels()
    message("DynET excluded from main figures (INCLUDE_DYNET = FALSE).")
  }
}

read_sweep <- function(files, value_col) {
  paths <- file.path(csv_folder, files)
  ok <- file.exists(paths)
  if (!any(ok)) return(NULL)
  if (!all(ok)) warning("Missing sweep file(s): ",
                        paste(basename(paths[!ok]), collapse = ", "))
  map_df(paths[ok], ~ read_csv(.x, col_types = cols(.default = "d",
                                                    Project = "c",
                                                    Scenario = "c"))) %>%
    add_labels()   # sweep files carry no configuration columns
}



# CONFIGURATION GUARD
# Every row records the model settings that produced it. If a results file mixes
# runs from different settings the analysis is meaningless, so stop rather than
# produce a plausible-looking but invalid table.
# DD_Future_Cost and DD_Base_Kill only enter the calculation when a cell
# receives the reduced background dose, which never happens under ET, PB or the
# forward-looking ET, so they are checked among the DD rows only. Global columns
# must agree everywhere, and DD rows must not come from validation mode ('v13').
cfg_global <- intersect(c("Future_Detection", "ET_Future_Rule", "Recall",
                          "Annuity_Factor", "Survivor_Cond", "Alloc_Rule"),
                        names(clean_data))
cfg_dd     <- c("DD_Future_Cost", "DD_Base_Kill")
if (length(cfg_global) >= 4 && all(cfg_dd %in% names(clean_data))) {
  cfg <- clean_data %>% distinct(across(all_of(cfg_global)))
  if (nrow(cfg) > 1) {
    print(as.data.frame(cfg))
    stop("The results file mixes runs made with different model settings ",
         "(see the table above). Re-run the affected fields, or split the file.")
  }
  message("Model configuration for this run:")
  for (cc in cfg_global) message(sprintf("   %-18s %s", cc, cfg[[cc]][1]))
  if (!is.na(cfg$Future_Detection[1]) && cfg$Future_Detection[1] != "probabilistic")
    warning("Future_Detection is '", cfg$Future_Detection[1], "', not ",
            "'probabilistic'. These are validation-mode results.")
  # the recall curve is only meaningful with conditional survivors
  if ("Survivor_Cond" %in% names(cfg)) {
    if (!identical(as.character(cfg$Survivor_Cond[1]), "thinned"))
      warning("Survivor_Cond is '", cfg$Survivor_Cond[1], "'. The recall curve ",
              "does not converge at r = 1 under 'v16'; re-run with 'thinned'.")
  } else {
    warning("No 'Survivor_Cond' column: output of an earlier model version. ",
            "The recall curve file (if any) should not be used with these rows.")
  }
  # allocation of isolated undetected weeds
  if ("Alloc_Rule" %in% names(cfg)) {
    if (!identical(as.character(cfg$Alloc_Rule[1]), "area"))
      warning("Alloc_Rule is '", cfg$Alloc_Rule[1], "'. The reported results use ",
              "'area' (isolated escapes spread over the field area).")
  } else {
    warning("No 'Alloc_Rule' column: output of an earlier model version ",
            "(isolated escapes packed into zero-detection cells).")
  }
  
  dd_cfg <- clean_data %>%
    filter(Scenario_Label == "DD strategy") %>%
    distinct(across(all_of(cfg_dd)))
  if (nrow(dd_cfg) > 1) {
    print(as.data.frame(dd_cfg))
    stop("The dual-dose rows mix settings (see above). This usually means an ",
         "old DD row survived the substitution. Delete the rows whose ",
         "DD_Future_Cost is 'v13' and re-paste, or re-run those fields.")
  }
  if (nrow(dd_cfg) == 1) {
    message("Dual-dose configuration:")
    for (cc in cfg_dd) message(sprintf("   %-18s %s", cc, dd_cfg[[cc]][1]))
    if (identical(as.character(dd_cfg$DD_Future_Cost[1]), "v13"))
      stop("The dual-dose rows were produced with DD_FUTURE_COST_MODE = 'v13', ",
           "which charges the full herbicide cost in background-dose years and ",
           "is a validation setting only. Re-run SCENARIO_ID = 3 with ",
           "DD_FUTURE_COST_MODE = 'dose_scaled'.")
    DD_SIMULATED_K <- as.numeric(dd_cfg$DD_Base_Kill[1])
  } else {
    DD_SIMULATED_K <- NA_real_
  }
} else {
  warning("No configuration columns found - output of an earlier model version.")
  DD_SIMULATED_K <- NA_real_
}


# COMPLETENESS AND DUPLICATE CHECK
#
# The main file is built by pasting re-run rows over old ones, so the two things
# that can silently go wrong are a row deleted and not replaced, and an old row
# left behind. Both produce means computed over the wrong number of fields
# rather than an error, so they are checked explicitly.
combo_check <- clean_data %>%
  count(Field_Name, Resolution, Scenario_Label, name = "n_rows") %>%
  complete(Field_Name, Resolution, Scenario_Label, fill = list(n_rows = 0L))

missing_combos <- combo_check %>% filter(n_rows == 0)
dupe_combos    <- combo_check %>% filter(n_rows > 1)

n_fields_found <- n_distinct(clean_data$Field_Name)
if (n_fields_found != EXPECTED_N_FIELDS)
  warning(sprintf("Found %d fields, expected %d.", n_fields_found,
                  EXPECTED_N_FIELDS))

if (nrow(dupe_combos) > 0) {
  print(as.data.frame(dupe_combos))
  stop("Duplicate rows for the combinations above. Keep exactly one row per ",
       "field x resolution x strategy.")
}
if (nrow(missing_combos) > 0) {
  print(as.data.frame(missing_combos))
  stop("No row for the combinations above. Every mean would be computed over ",
       "a different number of fields, so the strategies would not be ",
       "comparable. Re-run the missing combination(s) and paste the row in.")
}
message(sprintf("Completeness check passed: %d fields x %d resolutions x %d strategies.",
                n_fields_found, n_distinct(clean_data$Resolution),
                n_distinct(clean_data$Scenario_Label)))

# VERIFY OR SUBSTITUTE THE DUAL-DOSE ROWS
#
# The main file is simulated at the reported background efficacy, so the
# normal path is verification: the DD rows must reproduce the
# matching row of the background-efficacy sweep. Substitution is retained only
# for the case where a different efficacy is to be reported without re-running.
verify_dd <- function(df, k) {
  sw <- read_sweep(dd_sweep_files)
  if (is.null(sw)) {
    warning("DD sweep files not found - cannot verify the dual-dose rows.")
    return(invisible(NULL))
  }
  sel <- sw %>%
    filter(abs(K_Background - k) < 1e-9) %>%
    select(Field_Name, Resolution, sw_herb = D_Herb, sw_yield = D_Yield,
           sw_liab = D_Liab, sw_net = Net)
  cmp <- df %>%
    filter(Scenario_Label == "DD strategy") %>%
    select(Field_Name, Resolution, D_Herb_Saved, D_Yield_Yr0, D_Liab, D_Net) %>%
    inner_join(sel, by = c("Field_Name", "Resolution"))
  if (nrow(cmp) == 0) {
    warning("No matching sweep rows at k = ", k, " - dual-dose not verified.")
    return(invisible(NULL))
  }
  worst <- max(abs(cmp$D_Herb_Saved - cmp$sw_herb),
               abs(cmp$D_Yield_Yr0  - cmp$sw_yield),
               abs(cmp$D_Liab       - cmp$sw_liab),
               abs(cmp$D_Net        - cmp$sw_net), na.rm = TRUE)
  message(sprintf(
    "Dual-dose verification (main rows vs sweep at k = %.2f, %d rows): max deviation %.4f EUR/ha.",
    k, nrow(cmp), worst))
  # The main file rounds to 2 dp, the sweep writes 4 dp, so a deviation up to
  # about 0.005 is rounding and anything larger is a real mismatch.
  if (worst > 0.02)
    stop("The dual-dose rows do not reproduce the background-efficacy sweep. ",
         "The main file and the sweep file are probably from different runs.")
  invisible(NULL)
}

substitute_dd <- function(df, k) {
  if (is.na(k)) {
    message("Dual-dose reported exactly as simulated (DD_REPORT_K = NA).")
    return(df)
  }
  if (!is.na(DD_SIMULATED_K) && abs(k - DD_SIMULATED_K) < 1e-9) {
    message(sprintf(
      "Dual-dose reported at the simulated background efficacy (%.2f); verifying instead of substituting.",
      k))
    verify_dd(df, k)
    return(df)
  }
  warning("DD_REPORT_K (", k, ") differs from the simulated DD_Base_Kill (",
          DD_SIMULATED_K, "). Only the four decomposition columns can be ",
          "substituted, so the DD sensitivity columns will belong to the ",
          "simulated specification, not the reported one.")
  sw <- read_sweep(dd_sweep_files)
  if (is.null(sw)) {
    warning("DD sweep files not found - keeping the simulated 0.95 values.")
    return(df)
  }
  sel <- sw %>%
    filter(abs(K_Background - k) < 1e-9) %>%
    select(Field_Name, Resolution, D_Herb, D_Yield, D_Liab, Net)
  if (nrow(sel) == 0) {
    warning("Background efficacy ", k, " not present in the sweep - ",
            "keeping 0.95. Available: ",
            paste(sort(unique(sw$K_Background)), collapse = ", "))
    return(df)
  }
  n_before <- sum(df$Scenario_Label == "DD strategy")
  # rename the incoming sweep columns BEFORE joining. Joining a column
  # called D_Liab onto a frame that already has one produces D_Liab.x / D_Liab.y
  # and a later rename() then collides.
  sel <- sel %>%
    rename(sw_herb = D_Herb, sw_yield = D_Yield,
           sw_liab = D_Liab, sw_net = Net)
  out <- df %>%
    left_join(sel, by = c("Field_Name", "Resolution")) %>%
    mutate(
      is_dd = Scenario_Label == "DD strategy" & !is.na(sw_net),
      D_Herb_Saved   = ifelse(is_dd, sw_herb,  D_Herb_Saved),
      D_Yield_Yr0    = ifelse(is_dd, sw_yield, D_Yield_Yr0),
      D_Liab         = ifelse(is_dd, sw_liab,  D_Liab),
      D_Net          = ifelse(is_dd, sw_net,   D_Net),
      Tech_Breakeven = ifelse(is_dd, sw_net,   Tech_Breakeven),
      SSWM_Liab_Base = ifelse(is_dd, -sw_liab, SSWM_Liab_Base),
      SSWM_Net_Base  = ifelse(is_dd, BC_Net + sw_net, SSWM_Net_Base)
    ) %>%
    select(-sw_herb, -sw_yield, -sw_liab, -sw_net, -is_dd)
  message(sprintf("Dual-dose rows substituted from the sweep at background efficacy %.2f (%d rows).",
                  k, n_before))
  out
}

clean_data <- substitute_dd(clean_data, DD_REPORT_K)

# Detect whether the decomposition columns are present.
has_decomp <- all(c("D_Herb_Saved", "D_Yield_Yr0", "D_Liab", "D_Net",
                    "Tech_Breakeven") %in% names(clean_data))
has_descr  <- all(c("N_Cells", "Det_Dicot_Per_Ha", "Det_Mono_Per_Ha",
                    "Pct_Cells_Infested") %in% names(clean_data))
if (!has_decomp) warning("Decomposition columns absent — the decomposition figure and the ",
                         "break-even outputs will be skipped. Re-run the ",
                         "simulation model.")


# =============================================================================
# RECALL SWEEP (validation design, used only if the recall curve is missing)
#
# The sweep file is in long format, one row per (field, resolution, strategy,
# recall). 'Net' is the difference vs broadcast BEFORE technology cost, so it is
# directly comparable to D_Net; technology cost is subtracted below.
# =============================================================================

recall_sweep <- read_sweep(recall_sweep_files)

# apply the same DynET filter to the sweep so it matches clean_data
if (!is.null(recall_sweep) && !INCLUDE_DYNET) {
  recall_sweep <- recall_sweep %>%
    filter(Scenario_Label != "DynET strategy") %>% droplevels()
}

# Guard: at a background efficacy below 0.95 the DD strategy responds to recall
# (a missed weed in a background-dose cell survives at 1 - DD_BASE_KILL), so a
# flat DD response indicates a run with DD_BASE_KILL = 0.95.
if (!is.null(recall_sweep)) {
  dd_range <- recall_sweep %>%
    filter(Scenario_Label == "DD strategy") %>%
    group_by(Field_Name, Resolution) %>%
    summarise(rng = max(Net) - min(Net), .groups = "drop")
  if (nrow(dd_range) > 0) {
    message(sprintf(
      "Dual-dose recall sensitivity: mean range %.2f EUR/ha, max %.2f EUR/ha across fields.",
      mean(dd_range$rng, na.rm = TRUE), max(dd_range$rng, na.rm = TRUE)))
    if (max(dd_range$rng, na.rm = TRUE) < 0.01)
      warning("The dual-dose rows do not respond to recall at all. That is the ",
              "signature of a run made with DD_BASE_KILL = 0.95; check the ",
              "DD_BASE_KILL used for the sweep file.")
  }
}

use_old_sweep <- !is.null(recall_sweep) && nrow(recall_sweep) > 0

if (use_old_sweep) {
  message(sprintf("Validation recall sweep loaded (fallback only): %d rows, recall values = %s",
                  nrow(recall_sweep),
                  paste(sort(unique(recall_sweep$Recall)), collapse = ", ")))
}


# =============================================================================
# RECALL CURVE: sensor of recall r' on the same population
#
# One row per (field, resolution, strategy, recall). 'Net' is the difference vs
# broadcast BEFORE technology cost, i.e. the break-even technology cost at that
# recall. The curve at the measured recall must reproduce the baseline
# decomposition of the main file, which is checked here.
# =============================================================================

recall_curve <- read_sweep(recall_curve_files)
if (!is.null(recall_curve) && !INCLUDE_DYNET) {
  recall_curve <- recall_curve %>%
    filter(Scenario_Label != "DynET strategy") %>% droplevels()
}
use_curve <- !is.null(recall_curve) && nrow(recall_curve) > 0

if (use_curve) {
  curve_recalls <- sort(unique(recall_curve$Recall))
  message(sprintf("Recall CURVE loaded: %d rows, %d recall values (%.3f to %.3f).",
                  nrow(recall_curve), length(curve_recalls),
                  min(curve_recalls), max(curve_recalls)))

  # completeness: every field x resolution x strategy must have the full grid
  curve_n <- recall_curve %>%
    count(Field_Name, Resolution, Scenario_Label, name = "n_rec")
  if (any(curve_n$n_rec != length(curve_recalls))) {
    print(as.data.frame(curve_n %>% filter(n_rec != length(curve_recalls))))
    stop("The recall curve file has an incomplete recall grid for the rows above ",
         "(duplicate or missing runs). Keep one run per combination.")
  }
  curve_missing <- clean_data %>%
    distinct(Field_Name, Resolution, Scenario_Label) %>%
    anti_join(curve_n, by = c("Field_Name", "Resolution", "Scenario_Label"))
  if (nrow(curve_missing) > 0) {
    print(as.data.frame(curve_missing))
    stop("No recall-curve rows for the combinations above.")
  }

  # validation against the main file at the measured recall
  measured_recall <- if ("Recall" %in% names(clean_data)) clean_data$Recall[1] else 0.41
  r_base <- curve_recalls[which.min(abs(curve_recalls - measured_recall))]
  if (abs(r_base - measured_recall) > 1e-6)
    stop(sprintf("The curve grid does not contain the measured recall (%.3f).",
                 measured_recall))
  curve_chk <- recall_curve %>%
    filter(abs(Recall - r_base) < 1e-9) %>%
    select(Field_Name, Resolution, Scenario_Label,
           cv_herb = D_Herb, cv_yield = D_Yield, cv_liab = D_Liab, cv_net = Net) %>%
    inner_join(clean_data %>%
                 select(Field_Name, Resolution, Scenario_Label,
                        D_Herb_Saved, D_Yield_Yr0, D_Liab, D_Net),
               by = c("Field_Name", "Resolution", "Scenario_Label"))
  curve_dev <- with(curve_chk, max(abs(cv_herb - D_Herb_Saved), abs(cv_yield - D_Yield_Yr0),
                                   abs(cv_liab - D_Liab), abs(cv_net - D_Net)))
  message(sprintf("Recall CURVE validation (curve@%.3f vs main file, %d rows): max deviation %.4f EUR/ha",
                  r_base, nrow(curve_chk), curve_dev))
  if (curve_dev > 0.02)
    stop("The recall curve at the measured recall does not reproduce the main ",
         "file. The two files are from different runs or different settings.")

  # nearest grid values for the +/-20% sensitivity
  r_lo <- curve_recalls[which.min(abs(curve_recalls - RECALL_SENS_LOW))]
  r_hi <- curve_recalls[which.min(abs(curve_recalls - RECALL_SENS_HIGH))]
  if (abs(r_lo - RECALL_SENS_LOW) > 1e-6 || abs(r_hi - RECALL_SENS_HIGH) > 1e-6)
    message(sprintf("  NOTE: sensitivity points not on the grid; using %.3f and %.3f instead.",
                    r_lo, r_hi))
} else {
  if (use_old_sweep) {
    warning("Recall curve file not found. The +/-20% recall sensitivity falls ",
            "back to the validation design (spray map fixed, population varied). ",
            "No break-even recall can be computed. Run the model with ",
            "ENABLE_RECALL_CURVE = True.")
  } else {
    warning("Neither recall file found. Falling back to the liability-only ",
            "recall columns of the main file, which understate the influence of ",
            "recall because Year-0 revenue was held fixed.")
  }
}
use_new_recall <- use_curve || use_old_sweep


# =============================================================================
# TECHNOLOGY COST APPLICATION AND DIFF CALCULATIONS
#
# DIFF = Python_Net_Column - Tech_Cost_EUR - BC_Net  (single subtraction)
# =============================================================================

message("Applying technology cost scenarios...")

analysis_data <- clean_data %>%
  left_join(tech_cost_table, by = "Resolution", relationship = "many-to-many") %>%
  mutate(
    Tech_Cost_Label = factor(Tech_Cost_Label, levels = c(
      "Costs scenario 1", "Costs scenario 2", "Costs scenario 3",
      "Costs scenario 4", "Costs scenario 5"
    )),
    DIFF_Total      = SSWM_Net_Base   - Tech_Cost_EUR - BC_Net,
    DIFF_Weed_Opt   = Weed_Low_Net    - Tech_Cost_EUR - BC_Net,
    DIFF_Weed_Pess  = Weed_High_Net   - Tech_Cost_EUR - BC_Net,
    DIFF_Split_9010 = Split_90_10_Net - Tech_Cost_EUR - BC_Net,
    DIFF_Split_5050 = Split_50_50_Net - Tech_Cost_EUR - BC_Net,
    DIFF_Kill_Low   = Kill_Low_Net    - Tech_Cost_EUR - BC_Net,
    DIFF_Kill_High  = Kill_High_Net   - Tech_Cost_EUR - BC_Net,
    # liability-only recall columns of the main file, for comparison only
    DIFF_Recall_Low_OLD  = Recall_Low_Net  - Tech_Cost_EUR - BC_Net,
    DIFF_Recall_Base_OLD = Recall_Base_Net - Tech_Cost_EUR - BC_Net,
    DIFF_Recall_High_OLD = Recall_High_Net - Tech_Cost_EUR - BC_Net
  ) %>%
  # guarantee the factor survives the join
  mutate(Resolution = factor(as.character(Resolution), levels = RES_LEVELS))

# Validation: Recall_Base must reproduce the baseline
# The dual-dose rows are no longer excluded. The legacy
# Recall_*_Net columns are now written at the same background efficacy as the
# baseline, so DD must satisfy this identity like every other strategy. It is
# skipped only if DD values were substituted from a different specification.
dd_substituted <- !is.na(DD_REPORT_K) && !is.na(DD_SIMULATED_K) &&
  abs(DD_REPORT_K - DD_SIMULATED_K) > 1e-9
val_diff <- analysis_data %>%
  filter(!dd_substituted | Scenario_Label != "DD strategy") %>%
  summarise(max_diff = max(abs(DIFF_Recall_Base_OLD - DIFF_Total), na.rm = TRUE))
message(sprintf("Recall baseline validation — max deviation: %.4f EUR/ha",
                val_diff$max_diff))
if (val_diff$max_diff > 0.01)
  warning("Recall_Base deviates from DIFF_Total — check simulation output.")

# attach the corrected recall bounds
# From the curve (population fixed, sensor varied) when available; the
# +/-20% points are r_lo / r_hi, NOT the ends of the grid (the curve runs to 1).
if (use_new_recall) {
  recall_src <- if (use_curve) recall_curve else recall_sweep
  if (!use_curve) {
    r_lo <- min(recall_src$Recall); r_hi <- max(recall_src$Recall)
  }
  recall_wide <- recall_src %>%
    left_join(tech_cost_table, by = "Resolution", relationship = "many-to-many") %>%
    mutate(Tech_Cost_Label = factor(Tech_Cost_Label,
                                    levels = levels(analysis_data$Tech_Cost_Label)),
           DIFF = Net - Tech_Cost_EUR) %>%
    group_by(Field_Name, Resolution, Scenario_Label, Tech_Cost_Label) %>%
    summarise(
      DIFF_Recall_Low  = DIFF[which.min(abs(Recall - r_lo))],   # pessimistic sensor
      DIFF_Recall_Base = DIFF[which.min(abs(Recall - 0.41))],
      DIFF_Recall_High = DIFF[which.min(abs(Recall - r_hi))],   # optimistic sensor
      .groups = "drop"
    )
  
  analysis_data <- analysis_data %>%
    left_join(recall_wide,
              by = c("Field_Name", "Resolution", "Scenario_Label", "Tech_Cost_Label"))
  
  chk <- analysis_data %>%
    summarise(m = max(abs(DIFF_Recall_Base - DIFF_Total), na.rm = TRUE))
  message(sprintf("Recall SWEEP validation (sweep@0.41 vs baseline): %.4f EUR/ha",
                  chk$m))
  if (isTRUE(chk$m > 0.01))
    warning("Recall sweep at 0.41 does not reproduce the baseline. ",
            "Check that ANNUITY_FACTOR is pinned and the sweep ran on the same data.")
} else {
  analysis_data <- analysis_data %>%
    mutate(DIFF_Recall_Low  = DIFF_Recall_Low_OLD,
           DIFF_Recall_Base = DIFF_Recall_Base_OLD,
           DIFF_Recall_High = DIFF_Recall_High_OLD)
}


# =============================================================================
# STATISTICAL ANALYSIS — LINEAR MIXED MODEL
#
# Each of the 19 fields is simulated under every strategy x resolution
# combination, so observations are paired within field. A one-way ANOVA treats
# them as independent and understates the residual structure. Field enters as a
# random intercept; letters compare resolutions within each strategy.
#
# Vol_Red_Pct is a bounded percentage; a linear model is used as an approximation
# because values sit well inside (0, 100).
# =============================================================================

message("Fitting mixed model for herbicide savings...")

get_mixed_letters <- function(data) {
  if (!(have_lme4 && have_emmeans && have_multcomp)) {
    warning("lme4/emmeans/multcomp not available — falling back to per-strategy ",
            "aov(). Install them before producing the final figure.")
    return(
      data %>%
        group_by(Scenario_Label) %>%
        do({
          m <- aov(Vol_Red_Pct ~ Resolution, data = .)
          tk <- TukeyHSD(m)$Resolution
          tibble(Resolution = levels(droplevels(.$Resolution)), groups = "")
        }) %>% ungroup()
    )
  }
  
  fit <- lme4::lmer(Vol_Red_Pct ~ Resolution * Scenario_Label + (1 | Field_Name),
                    data = data)
  print(summary(fit))
  if (have_lmerTest) {
    cat("\nType III ANOVA (Satterthwaite):\n")
    print(anova(fit))
  }
  
  emm <- emmeans::emmeans(fit, ~ Resolution | Scenario_Label)
  # decreasing = TRUE assigns "a" to the highest mean.
  # adjust = "sidak" is set explicitly: emmeans silently substitutes it for
  # "tukey" here (tukey applies to a single set of pairwise comparisons only).
  cld <- multcomp::cld(emm, Letters = letters, adjust = "sidak",
                       sort = TRUE, decreasing = TRUE)
  
  as.data.frame(cld) %>%
    transmute(Scenario_Label,
              Resolution = factor(as.character(Resolution), levels = RES_LEVELS),
              groups = trimws(as.character(.group)))
}

vol_letters <- get_mixed_letters(clean_data)
message("\nCompact letter display (resolutions compared within strategy):")
print(as.data.frame(vol_letters), row.names = FALSE)
write_csv(vol_letters, file.path(csv_folder, "Stats_mixed_model_letters.csv"))


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

message("Computing summary statistics...")

summary_stats <- analysis_data %>%
  group_by(Tech_Cost_Label, Scenario_Label, Resolution) %>%
  summarise(
    Econ_Mean        = mean(DIFF_Total,       na.rm = TRUE),
    Econ_SD          = sd(DIFF_Total,         na.rm = TRUE),
    Econ_Min         = min(DIFF_Total,        na.rm = TRUE),
    Econ_Max         = max(DIFF_Total,        na.rm = TRUE),
    N_Fields         = n(),
    N_Positive       = sum(DIFF_Total > 0,    na.rm = TRUE),
    Econ_Weed_Opt    = mean(DIFF_Weed_Opt,    na.rm = TRUE),
    Econ_Weed_Pess   = mean(DIFF_Weed_Pess,   na.rm = TRUE),
    Econ_Split_9010  = mean(DIFF_Split_9010,  na.rm = TRUE),
    Econ_Split_5050  = mean(DIFF_Split_5050,  na.rm = TRUE),
    Econ_Kill_Low    = mean(DIFF_Kill_Low,    na.rm = TRUE),
    Econ_Kill_High   = mean(DIFF_Kill_High,   na.rm = TRUE),
    Econ_Recall_Low  = mean(DIFF_Recall_Low,  na.rm = TRUE),
    Econ_Recall_High = mean(DIFF_Recall_High, na.rm = TRUE),
    Vol_Mean = mean(Vol_Red_Pct, na.rm = TRUE),
    Vol_SD   = sd(Vol_Red_Pct,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # SENSITIVITY BOUNDS MUST BE ORDERED *AND* MUST CONTAIN THE BASELINE.
  #
  # Two problems appeared once the model was recalibrated:
  #
  #  1. The weed-dynamics bounds were passed to the figure as raw
  #     (pessimistic, optimistic) rather than (min, max). Reducing the weed
  #     parameters now produces the WORST outcome, so the pair became inverted
  #     and the error bar was drawn upside down.
  #
  #  2. The response to a parameter is not always monotone, so the baseline can
  #     lie OUTSIDE the interval spanned by the low and high settings. The bar
  #     (the baseline) then sits outside its own error bar. This happens for the
  #     weed-dynamics dimension under ET at 1 x 1 m and for the false-negative
  #     allocation at 50 x 50 cm.
  #
  # Both are handled by defining the interval as the full range across ALL
  # simulated settings for that dimension, INCLUDING the baseline. The bar is
  # then inside its interval by construction, and the interval means what the
  # caption says: the range of outcomes across the settings examined.
  mutate(
    Econ_Weed_Lo   = pmin(Econ_Weed_Opt, Econ_Weed_Pess, Econ_Mean),
    Econ_Weed_Hi   = pmax(Econ_Weed_Opt, Econ_Weed_Pess, Econ_Mean),
    Econ_Split_Lo  = pmin(Econ_Split_5050, Econ_Split_9010, Econ_Mean),
    Econ_Split_Hi  = pmax(Econ_Split_5050, Econ_Split_9010, Econ_Mean),
    Econ_Kill_Lo   = pmin(Econ_Kill_Low,   Econ_Kill_High,  Econ_Mean),
    Econ_Kill_Hi   = pmax(Econ_Kill_Low,   Econ_Kill_High,  Econ_Mean),
    Econ_Recall_Lo = pmin(Econ_Recall_Low, Econ_Recall_High, Econ_Mean),
    Econ_Recall_Hi = pmax(Econ_Recall_Low, Econ_Recall_High, Econ_Mean),
    # flag the non-monotone cases so they can be explained in the text
    NonMono_Weed  = (Econ_Mean < pmin(Econ_Weed_Opt, Econ_Weed_Pess)) |
      (Econ_Mean > pmax(Econ_Weed_Opt, Econ_Weed_Pess)),
    NonMono_Split = (Econ_Mean < pmin(Econ_Split_5050, Econ_Split_9010)) |
      (Econ_Mean > pmax(Econ_Split_5050, Econ_Split_9010))
  ) %>%
  left_join(vol_letters, by = c("Scenario_Label", "Resolution")) %>%
  rename(Vol_Group = groups)

dat <- summary_stats


# =============================================================================
# TABLE 5 — EXPECTED FUTURE WEED PENALTY: WEED DYNAMICS SENSITIVITY
# =============================================================================

# The FWP sensitivity tables include all three strategies. For the DD strategy
# the herbicide-efficacy dimension varies the full-rate efficacy while the
# background dose is held at its own value, so it measures the gap between the
# two doses; the background dose itself is varied in the separate sweep. The
# tables report -SSWM_Liab_Base: where the DD strategy reduces the liability
# relative to broadcast, the value is positive (a benefit).
fwp_sens_data <- clean_data

message("Generating Table 5...")

# COLUMN NAMES CHANGED FROM optimistic/pessimistic TO low/high.
#
# Those labels are now actively misleading. Reducing the weed population
# parameters by 20% makes broadcast application MORE effective at suppressing
# the population, while cells that site-specific application never detects
# continue to escape. The GAP between the two therefore WIDENS, so the
# "optimistic" weed scenario produces the LARGEST penalty for SSWM. The response
# is also non-monotonic at 1 x 1 m, because the population model is close to
# bistable in this parameter region.
#
# The columns are named for what was varied, not for the direction of the
# economic outcome.
table5 <- fwp_sens_data %>%
  group_by(Scenario_Label, Resolution) %>%
  summarise(
    FWP_Params_Low_20pct  = round(-mean(Weed_Low_Liab,  na.rm = TRUE), 2),
    FWP_Params_Baseline   = round(-mean(SSWM_Liab_Base, na.rm = TRUE), 2),
    FWP_Params_High_20pct = round(-mean(Weed_High_Liab, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(Resolution, Scenario_Label) %>%
  mutate(Monotonic = ifelse(
    (FWP_Params_Low_20pct <= FWP_Params_Baseline &
       FWP_Params_Baseline <= FWP_Params_High_20pct) |
      (FWP_Params_Low_20pct >= FWP_Params_Baseline &
         FWP_Params_Baseline >= FWP_Params_High_20pct), "yes", "NO"))

print(as.data.frame(table5), row.names = FALSE)
write_csv(table5, file.path(csv_folder, "Table5_FWP_weed_dynamics.csv"))


# =============================================================================
# TABLE 6 / FIGURE 5 — PARTIAL BUDGET DECOMPOSITION
#
# The four components of the difference vs broadcast:
#   (a) herbicide cost saved      (positive)
#   (b) Year-0 yield loss         (negative)
#   (c) future weed penalty       (negative)
#   (d) technology cost           (negative)
#   (e) net                       (a+b+c+d)
#
# This replaces the assertion that losses are "driven by escapes rather than
# technology cost" with the arithmetic that demonstrates it.
# =============================================================================

if (has_decomp) {
  message("Generating decomposition table and Figure 5...")
  
  decomp <- analysis_data %>%
    group_by(Tech_Cost_Label, Scenario_Label, Resolution) %>%
    summarise(
      `Herbicide saved`     = mean(D_Herb_Saved,  na.rm = TRUE),
      `Year-0 yield loss`   = mean(D_Yield_Yr0,   na.rm = TRUE),
      `Future weed penalty` = mean(D_Liab,        na.rm = TRUE),
      `Technology cost`     = -mean(Tech_Cost_EUR, na.rm = TRUE),
      Net                   = mean(D_Net - Tech_Cost_EUR, na.rm = TRUE),
      .groups = "drop"
    )
  
  decomp_out <- decomp %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    arrange(Tech_Cost_Label, Resolution, Scenario_Label)
  write_csv(decomp_out, file.path(csv_folder, "Table6_partial_budget_decomposition.csv"))
  message("\nPartial budget decomposition (", main_text_cost_scenario, "):")
  print(as.data.frame(decomp_out %>% filter(Tech_Cost_Label == main_text_cost_scenario)),
        row.names = FALSE)
  
  # COMPACT MAIN-TEXT TABLE.
  # The full table has 45 rows (3 strategies x 3 resolutions x 5 cost scenarios),
  # which is too long for a results section. This 9-row version fixes the cost
  # scenario and shows only the four components plus net; the full 45-row table
  # goes to supplementary.
  decomp_compact <- decomp %>%
    filter(Tech_Cost_Label == main_text_cost_scenario) %>%
    select(-Tech_Cost_Label) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    arrange(Resolution, Scenario_Label)
  write_csv(decomp_compact,
            file.path(csv_folder, "Table6_decomposition_MAINTEXT_9row.csv"))
  
  decomp_long <- decomp %>%
    pivot_longer(cols = c(`Herbicide saved`, `Year-0 yield loss`,
                          `Future weed penalty`, `Technology cost`),
                 names_to = "Component", values_to = "Value") %>%
    mutate(Component = factor(Component, levels = names(component_colors)))
}


# =============================================================================
# BREAK-EVEN TECHNOLOGY COST
#
# The maximum annualised technology + information cost at which SSWM matches
# broadcast. Needs no technology cost assumption, so it is reported once per
# strategy x resolution and compared against the five scenarios of Table 4.
# =============================================================================

if (has_decomp) {
  message("Generating break-even technology cost table...")
  
  breakeven <- clean_data %>%
    group_by(Scenario_Label, Resolution) %>%
    summarise(
      BE_Mean   = round(mean(Tech_Breakeven, na.rm = TRUE), 2),
      BE_SD     = round(sd(Tech_Breakeven,   na.rm = TRUE), 2),
      BE_Min    = round(min(Tech_Breakeven,  na.rm = TRUE), 2),
      BE_Median = round(median(Tech_Breakeven, na.rm = TRUE), 2),
      BE_Max    = round(max(Tech_Breakeven,  na.rm = TRUE), 2),
      N_Fields_Positive = sum(Tech_Breakeven > 0, na.rm = TRUE),
      N_Fields  = n(),
      .groups = "drop"
    ) %>%
    left_join(
      tech_cost_table %>%
        group_by(Resolution) %>%
        summarise(Cheapest_Tech_Scenario = min(Tech_Cost_EUR),
                  Dearest_Tech_Scenario  = max(Tech_Cost_EUR), .groups = "drop"),
      by = "Resolution"
    ) %>%
    mutate(Viable_vs_Cheapest = ifelse(BE_Mean >= Cheapest_Tech_Scenario,
                                       "yes", "no")) %>%
    arrange(Resolution, Scenario_Label)
  
  message("\nBreak-even technology cost (EUR/ha):")
  print(as.data.frame(breakeven), row.names = FALSE)
  write_csv(breakeven, file.path(csv_folder, "Table7_breakeven_tech_cost.csv"))
}


# =============================================================================
# BREAK-EVEN RECALL
#
# The recall at which the net return difference vs broadcast (before or after
# technology cost) crosses zero, read off the recall curve by linear
# interpolation between grid points. Computed (a) on the MEAN curve across
# fields, which matches how the rest of the paper reports means, and (b) per
# field, to show the distribution. A "Zero technology cost" scenario is added so
# the free-technology case (the natural companion of a negative break-even
# cost) is reported alongside the five cost scenarios.
#
# Conventions in the per-field summary:
#   BE < lowest grid value  -> the field is viable even at the pessimistic
#                              sensor; reported as the lowest grid value.
#   no crossing up to 1.0   -> NA (not viable at any recall <= 1).
# =============================================================================

crossing_recall <- function(r, d) {
  o <- order(r); r <- r[o]; d <- d[o]
  ok <- is.finite(d); r <- r[ok]; d <- d[ok]
  if (length(d) == 0) return(NA_real_)
  if (d[1] >= 0) return(r[1])
  up <- which(d[-1] >= 0 & d[-length(d)] < 0)
  if (length(up) == 0) return(NA_real_)
  i <- up[1]
  r[i] + (0 - d[i]) * (r[i + 1] - r[i]) / (d[i + 1] - d[i])
}

cost_table_ext <- bind_rows(
  tibble(Tech_Cost_Label = "Zero technology cost",
         Resolution = factor(RES_LEVELS, levels = RES_LEVELS), Tech_Cost_EUR = 0),
  tech_cost_table
) %>%
  mutate(Tech_Cost_Label = factor(Tech_Cost_Label,
                                  levels = c("Zero technology cost",
                                             unique(as.character(tech_cost_table$Tech_Cost_Label)))))

if (use_curve) {
  message("Computing break-even recall...")

  curve_costed <- recall_curve %>%
    left_join(cost_table_ext, by = "Resolution", relationship = "many-to-many") %>%
    mutate(DIFF = Net - Tech_Cost_EUR)

  # (a) mean curve across fields
  curve_mean <- curve_costed %>%
    group_by(Scenario_Label, Resolution, Tech_Cost_Label, Tech_Cost_EUR, Recall) %>%
    summarise(Mean_DIFF = mean(DIFF, na.rm = TRUE),
              SD_DIFF   = sd(DIFF, na.rm = TRUE),
              Q25_DIFF  = quantile(DIFF, 0.25, na.rm = TRUE),
              Q75_DIFF  = quantile(DIFF, 0.75, na.rm = TRUE),
              Min_DIFF  = min(DIFF, na.rm = TRUE),
              Max_DIFF  = max(DIFF, na.rm = TRUE),
              N_Positive = sum(DIFF > 0, na.rm = TRUE),
              N_Fields  = n(),
              .groups = "drop")

  r_min <- min(curve_recalls)
  be_mean <- curve_mean %>%
    group_by(Scenario_Label, Resolution, Tech_Cost_Label, Tech_Cost_EUR) %>%
    summarise(BE_Recall_MeanCurve = crossing_recall(
                Recall[Recall <= BE_RECALL_MAX + 1e-9],
                Mean_DIFF[Recall <= BE_RECALL_MAX + 1e-9]),
              DIFF_at_measured    = Mean_DIFF[which.min(abs(Recall - r_base))],
              DIFF_at_1.0         = Mean_DIFF[which.max(Recall)],
              .groups = "drop") %>%
    mutate(BE_Status = case_when(
      is.na(BE_Recall_MeanCurve)              ~ sprintf("no break-even at any recall <= %.2f", BE_RECALL_MAX),
      BE_Recall_MeanCurve <= r_min + 1e-9     ~ sprintf("viable at every recall examined (>= %.3f)", r_min),
      BE_Recall_MeanCurve <= r_base + 1e-9    ~ "viable at the measured recall",
      TRUE                                    ~ "break-even above the measured recall"))

  # (b) per field
  be_field <- curve_costed %>%
    group_by(Field_Name, Scenario_Label, Resolution, Tech_Cost_Label, Tech_Cost_EUR) %>%
    summarise(BE_Recall = crossing_recall(Recall[Recall <= BE_RECALL_MAX + 1e-9],
                                          DIFF[Recall <= BE_RECALL_MAX + 1e-9]),
              DIFF_at_measured = DIFF[which.min(abs(Recall - r_base))],
              .groups = "drop")
  write_csv(be_field %>% arrange(Tech_Cost_Label, Resolution, Scenario_Label, Field_Name),
            file.path(csv_folder, "SuppTable_breakeven_recall_by_field.csv"))

  be_field_sum <- be_field %>%
    group_by(Scenario_Label, Resolution, Tech_Cost_Label, Tech_Cost_EUR) %>%
    summarise(
      N_Fields              = n(),
      N_Viable_At_Measured  = sum(DIFF_at_measured >= 0, na.rm = TRUE),
      N_With_BE_Below_1     = sum(!is.na(BE_Recall)),
      BE_Recall_Field_Median = median(BE_Recall, na.rm = TRUE),
      BE_Recall_Field_Min    = suppressWarnings(min(BE_Recall, na.rm = TRUE)),
      BE_Recall_Field_Max    = suppressWarnings(max(BE_Recall, na.rm = TRUE)),
      .groups = "drop") %>%
    mutate(across(c(BE_Recall_Field_Min, BE_Recall_Field_Max),
                  ~ ifelse(is.finite(.x), .x, NA_real_)))

  breakeven_recall <- be_mean %>%
    left_join(be_field_sum,
              by = c("Scenario_Label", "Resolution", "Tech_Cost_Label", "Tech_Cost_EUR")) %>%
    mutate(across(c(BE_Recall_MeanCurve, BE_Recall_Field_Median,
                    BE_Recall_Field_Min, BE_Recall_Field_Max), ~ round(.x, 3)),
           across(c(DIFF_at_measured, DIFF_at_1.0), ~ round(.x, 2))) %>%
    select(Tech_Cost_Label, Tech_Cost_EUR, Resolution, Scenario_Label,
           BE_Recall_MeanCurve, BE_Status, DIFF_at_measured, DIFF_at_1.0,
           N_Fields, N_Viable_At_Measured, N_With_BE_Below_1,
           BE_Recall_Field_Median, BE_Recall_Field_Min, BE_Recall_Field_Max) %>%
    arrange(Tech_Cost_Label, Resolution, Scenario_Label)

  message(sprintf("\nBreak-even recall (mean curve; NA = no break-even at any recall <= %.2f):",
                  BE_RECALL_MAX))
  print(as.data.frame(breakeven_recall %>%
                        filter(Tech_Cost_Label %in% c("Zero technology cost",
                                                      main_text_cost_scenario)) %>%
                        select(Tech_Cost_Label, Resolution, Scenario_Label,
                               BE_Recall_MeanCurve, BE_Status, BE_Recall_Field_Median,
                               N_With_BE_Below_1, N_Viable_At_Measured)),
        row.names = FALSE)
  write_csv(breakeven_recall, file.path(csv_folder, "Table8_breakeven_recall.csv"))

  # means along the curve, for the supplementary (herbicide saving, yield loss,
  # penalty and net return at every recall value, cost-free)
  curve_means_table <- recall_curve %>%
    group_by(Scenario_Label, Resolution, Recall) %>%
    summarise(Mean_Herb_Saved = round(mean(D_Herb), 2),
              Mean_Yield_Yr0  = round(mean(D_Yield), 2),
              Mean_FWP        = round(mean(D_Liab), 2),
              Mean_BE_Tech    = round(mean(Net), 2),
              SD_BE_Tech      = round(sd(Net), 2),
              Mean_Sprayed_Pct = round(mean(Sprayed_Pct), 1),
              Mean_Vol_Red_Pct = round(mean(Vol_Red_Pct), 1),
              .groups = "drop") %>%
    arrange(Resolution, Scenario_Label, Recall)
  write_csv(curve_means_table, file.path(csv_folder, "SuppTable_recall_curve_means.csv"))
} else {
  message("No recall curve - break-even recall skipped.")
}


# =============================================================================
# BETWEEN-FIELD DISTRIBUTION
#
# Counts how many fields were profitable, in addition to the means.
# =============================================================================

message("Generating between-field distribution table...")

field_spread <- analysis_data %>%
  group_by(Tech_Cost_Label, Scenario_Label, Resolution) %>%
  summarise(
    N_Fields    = n(),
    Net_Min     = round(min(DIFF_Total,    na.rm = TRUE), 2),
    Net_Q1      = round(quantile(DIFF_Total, 0.25, na.rm = TRUE), 2),
    Net_Median  = round(median(DIFF_Total, na.rm = TRUE), 2),
    Net_Q3      = round(quantile(DIFF_Total, 0.75, na.rm = TRUE), 2),
    Net_Max     = round(max(DIFF_Total,    na.rm = TRUE), 2),
    N_Positive  = sum(DIFF_Total > 0, na.rm = TRUE),
    Pct_Positive= round(100 * mean(DIFF_Total > 0, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(Tech_Cost_Label, Resolution, Scenario_Label)

write_csv(field_spread, file.path(csv_folder, "SuppTable_field_spread.csv"))
message("\nFields with positive net return (", main_text_cost_scenario, "):")
print(as.data.frame(field_spread %>%
                      filter(Tech_Cost_Label == main_text_cost_scenario) %>%
                      select(Scenario_Label, Resolution, N_Positive, N_Fields,
                             Net_Min, Net_Median, Net_Max)),
      row.names = FALSE)


# =============================================================================
# SUPPLEMENTARY TABLE — PER-FIELD DESCRIPTIVE STATISTICS
#
# Weed pressure and grid structure for each of the 19 fields. Values are
# identical across strategies, so one row per field x resolution.
# =============================================================================

if (has_descr) {
  message("Generating per-field descriptive table...")
  
  field_descr <- clean_data %>%
    distinct(Field_Name, Resolution, .keep_all = TRUE) %>%
    transmute(
      Field = Field_Name,
      Resolution,
      Area_ha            = round(Area_Ha, 2),
      N_grid_cells       = N_Cells,
      Dicots_per_ha      = round(Det_Dicot_Per_Ha, 1),
      Monocots_per_ha    = round(Det_Mono_Per_Ha, 1),
      Total_weeds_per_ha = round(Det_Dicot_Per_Ha + Det_Mono_Per_Ha, 1),
      Pct_dicot          = round(100 * Det_Dicot_Per_Ha /
                                   pmax(1e-9, Det_Dicot_Per_Ha + Det_Mono_Per_Ha), 1),
      Pct_cells_infested = round(Pct_Cells_Infested, 2)
    ) %>%
    arrange(Field, Resolution)
  
  write_csv(field_descr, file.path(csv_folder, "SuppTable_field_descriptives.csv"))
  message(sprintf("Per-field table: %d fields x %d resolutions = %d rows",
                  n_distinct(field_descr$Field),
                  n_distinct(field_descr$Resolution),
                  nrow(field_descr)))
  print(head(as.data.frame(field_descr), 6), row.names = FALSE)
  
  # Compact version: detections do not depend on resolution, so one row per
  # field is enough.
  field_descr_compact <- field_descr %>%
    group_by(Field) %>%
    summarise(
      Area_ha            = first(Area_ha),
      Dicots_per_ha      = first(Dicots_per_ha),
      Monocots_per_ha    = first(Monocots_per_ha),
      Total_weeds_per_ha = first(Total_weeds_per_ha),
      Pct_dicot          = first(Pct_dicot),
      Infested_1m        = Pct_cells_infested[Resolution == "1x1m"][1],
      Infested_50cm      = Pct_cells_infested[Resolution == "50x50cm"][1],
      Infested_25cm      = Pct_cells_infested[Resolution == "25x25cm"][1],
      .groups = "drop"
    )
  write_csv(field_descr_compact,
            file.path(csv_folder, "SuppTable_field_descriptives_compact.csv"))
}


# =============================================================================
# SUPPLEMENTARY TABLES 1-3 — SENSITIVITY (FWP basis)
# =============================================================================

supp_table1 <- fwp_sens_data %>%
  group_by(Scenario_Label, Resolution) %>%
  summarise(
    FWP_FN_9010 = round(-mean(Split_90_10_Liab, na.rm = TRUE), 2),
    FWP_FN_7030 = round(-mean(SSWM_Liab_Base,   na.rm = TRUE), 2),
    FWP_FN_5050 = round(-mean(Split_50_50_Liab, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>% arrange(Resolution, Scenario_Label)
write_csv(supp_table1, file.path(csv_folder, "SuppTable1_FN_allocation.csv"))

supp_table2 <- fwp_sens_data %>%
  group_by(Scenario_Label, Resolution) %>%
  summarise(
    FWP_Kill_Low  = round(-mean(Kill_Low_Liab,  na.rm = TRUE), 2),
    FWP_Kill_Base = round(-mean(Kill_Base_Liab, na.rm = TRUE), 2),
    FWP_Kill_High = round(-mean(Kill_High_Liab, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>% arrange(Resolution, Scenario_Label)
write_csv(supp_table2, file.path(csv_folder, "SuppTable2_herbicide_efficacy.csv"))

# Recall sensitivity table: net return differences from the recall curve.
# Column headers name the recall values.
if (use_new_recall) {
  rec_vals <- sort(unique(recall_src$Recall))
  supp_table3 <- recall_src %>%
    filter(abs(Recall - r_lo) < 1e-9 | abs(Recall - 0.41) < 1e-9 | abs(Recall - r_hi) < 1e-9) %>%
    left_join(tech_cost_table, by = "Resolution", relationship = "many-to-many") %>%
    filter(Tech_Cost_Label == main_text_cost_scenario) %>%
    mutate(DIFF = Net - Tech_Cost_EUR) %>%
    group_by(Scenario_Label, Resolution, Recall) %>%
    summarise(Mean_Net = round(mean(DIFF, na.rm = TRUE), 2), .groups = "drop") %>%
    pivot_wider(names_from = Recall, values_from = Mean_Net,
                names_prefix = "Recall_") %>%
    arrange(Resolution, Scenario_Label)
  
  message("\nSupplementary Table 3 (CORRECTED): net return by sensor recall, ",
          main_text_cost_scenario)
  print(as.data.frame(supp_table3), row.names = FALSE)
  write_csv(supp_table3, file.path(csv_folder, "SuppTable3_sensor_recall_CORRECTED.csv"))
  
  # Side-by-side comparison of the old and corrected recall ranges. Use this to
  # decide how to rewrite the sensitivity discussion.
  recall_compare <- analysis_data %>%
    filter(Tech_Cost_Label == main_text_cost_scenario) %>%
    group_by(Scenario_Label, Resolution) %>%
    summarise(
      Range_OLD_liab_only = round(mean(abs(DIFF_Recall_High_OLD -
                                             DIFF_Recall_Low_OLD)), 2),
      Range_NEW_with_rev  = round(mean(abs(DIFF_Recall_High - DIFF_Recall_Low)), 2),
      .groups = "drop"
    ) %>% arrange(Resolution, Scenario_Label)
  # flag non-monotonic responses so they are explained, not discovered
  mono <- recall_src %>%
    left_join(tech_cost_table, by = "Resolution", relationship = "many-to-many") %>%
    filter(Tech_Cost_Label == main_text_cost_scenario) %>%
    group_by(Scenario_Label, Resolution, Recall) %>%
    summarise(Net = mean(Net - Tech_Cost_EUR), .groups = "drop") %>%
    arrange(Scenario_Label, Resolution, Recall) %>%
    group_by(Scenario_Label, Resolution) %>%
    summarise(monotone = all(diff(Net) >= 0) || all(diff(Net) <= 0),
              .groups = "drop") %>%
    filter(!monotone)
  if (nrow(mono) > 0) {
    message("\nNON-MONOTONIC response to recall in these cases:")
    print(as.data.frame(mono), row.names = FALSE)
    message("  For the dual-dose strategy this is expected: a better sensor moves")
    message("  cells from the cheap background dose to the full dose (less herbicide")
    message("  saved) while reducing escapes.")
  }
  
  message("\nRecall sensitivity range: liability-only columns vs recall curve:")
  print(as.data.frame(recall_compare), row.names = FALSE)
  write_csv(recall_compare, file.path(csv_folder, "SuppTable3b_recall_old_vs_new.csv"))
} else {
  supp_table3 <- clean_data %>%
    group_by(Scenario_Label, Resolution) %>%
    summarise(
      FWP_Recall_0328 = round(-mean(Recall_Low_Liab,  na.rm = TRUE), 2),
      FWP_Recall_041  = round(-mean(Recall_Base_Liab, na.rm = TRUE), 2),
      FWP_Recall_0492 = round(-mean(Recall_High_Liab, na.rm = TRUE), 2),
      .groups = "drop"
    ) %>% arrange(Resolution, Scenario_Label)
  write_csv(supp_table3, file.path(csv_folder, "SuppTable3_sensor_recall.csv"))
}


# =============================================================================
# DUAL-DOSE BACKGROUND EFFICACY SWEEP
#
# Net return of the DD strategy as a function of the efficacy of the reduced
# background dose (Section 3.5).
# =============================================================================

dd_sweep <- read_sweep(dd_sweep_files)

if (!is.null(dd_sweep) && nrow(dd_sweep) > 0) {
  message("Generating dual-dose background efficacy table...")
  
  dd_table <- dd_sweep %>%
    left_join(tech_cost_table, by = "Resolution", relationship = "many-to-many") %>%
    mutate(Tech_Cost_Label = factor(Tech_Cost_Label,
                                    levels = levels(analysis_data$Tech_Cost_Label)),
           DIFF = Net - Tech_Cost_EUR) %>%
    group_by(Tech_Cost_Label, Resolution, K_Background) %>%
    summarise(
      Mean_Net    = round(mean(DIFF, na.rm = TRUE), 2),
      Mean_Yield  = round(mean(D_Yield, na.rm = TRUE), 2),
      Mean_Liab   = round(mean(D_Liab,  na.rm = TRUE), 2),
      N_Positive  = sum(DIFF > 0, na.rm = TRUE),
      N_Fields    = n(),
      .groups = "drop"
    ) %>% arrange(Tech_Cost_Label, Resolution, desc(K_Background))
  
  message("\nDual-dose background efficacy (", main_text_cost_scenario, "):")
  print(as.data.frame(dd_table %>% filter(Tech_Cost_Label == main_text_cost_scenario)),
        row.names = FALSE)
  write_csv(dd_table, file.path(csv_folder, "SuppTable_DD_background_efficacy.csv"))
} else {
  message("No dual-dose sweep files found — skipping DD efficacy outputs.")
}


# =============================================================================
# TORNADO DATA — all four sensitivity dimensions on one scale
#
# Answers "which parameter matters most" directly, instead of leaving the reader
# to compare ranges across four separate supplementary figures.
# =============================================================================

message("Building tornado comparison...")

# The DD strategy is included. For DD the herbicide-efficacy dimension measures
# the gap between the full and the background dose rather than herbicide
# efficacy as such.
#
# The tornado is grouped by strategy, so including DD adds a panel rather than
# diluting the other two. Set TORNADO_INCLUDE_DD to FALSE to restore the
# two-strategy version.
TORNADO_INCLUDE_DD <- TRUE

# Full spread across the three settings of a dimension: the low setting, the
# baseline, and the high setting. Means are taken first (across fields), then
# the span, so the result is the spread of the MEAN net return — the same
# quantity as before, just measured over three points instead of two.
span3 <- function(lo, base, hi) {
  v <- c(mean(lo,   na.rm = TRUE),
         mean(base, na.rm = TRUE),
         mean(hi,   na.rm = TRUE))
  max(v) - min(v)
}

# TRUE when the baseline falls outside the interval spanned by the low and high
# settings, i.e. when both perturbations move the result in the same direction.
nonmono3 <- function(lo, base, hi) {
  l <- mean(lo, na.rm = TRUE); b <- mean(base, na.rm = TRUE); h <- mean(hi, na.rm = TRUE)
  isTRUE(b < min(l, h) - 1e-9) || isTRUE(b > max(l, h) + 1e-9)
}

# DIFF_Total is the baseline for every dimension. (DIFF_Recall_Base is the
# baseline recorded in the recall sweep and is identical to DIFF_Total; the
# script validates this earlier, so either can be used here.)
tornado <- analysis_data %>%
  filter(TORNADO_INCLUDE_DD | Scenario_Label != "DD strategy") %>%
  group_by(Tech_Cost_Label, Scenario_Label, Resolution) %>%
  summarise(
    Baseline = mean(DIFF_Total, na.rm = TRUE),
    # labels wrapped onto two lines so they still fit once the text is
    # enlarged. The full descriptions belong in the caption, not the axis.
    `Weed population\ndynamics (±20%)` =
      span3(DIFF_Weed_Opt,   DIFF_Total, DIFF_Weed_Pess),
    `False-negative\nallocation (50–90%)` =
      span3(DIFF_Split_9010, DIFF_Total, DIFF_Split_5050),
    `Herbicide efficacy\n(k = 0.85–0.99)` =
      span3(DIFF_Kill_Low,   DIFF_Total, DIFF_Kill_High),
    `Sensor recall\n(±20%)` =
      span3(DIFF_Recall_Low, DIFF_Total, DIFF_Recall_High),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with(c("Weed", "False", "Herbicide", "Sensor")),
               names_to = "Dimension", values_to = "Range")
# NOTE: mutate(Range = abs(Range)) is no longer needed — span3() is a max minus
# a min and is therefore never negative.

# Report which bars come from a non-monotone response.
tornado_nonmono <- analysis_data %>%
  filter(TORNADO_INCLUDE_DD | Scenario_Label != "DD strategy") %>%
  group_by(Tech_Cost_Label, Scenario_Label, Resolution) %>%
  summarise(
    `Weed population\ndynamics (±20%)` =
      nonmono3(DIFF_Weed_Opt,   DIFF_Total, DIFF_Weed_Pess),
    `False-negative\nallocation (50–90%)` =
      nonmono3(DIFF_Split_9010, DIFF_Total, DIFF_Split_5050),
    `Herbicide efficacy\n(k = 0.85–0.99)` =
      nonmono3(DIFF_Kill_Low,   DIFF_Total, DIFF_Kill_High),
    `Sensor recall\n(±20%)` =
      nonmono3(DIFF_Recall_Low, DIFF_Total, DIFF_Recall_High),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with(c("Weed", "False", "Herbicide", "Sensor")),
               names_to = "Dimension", values_to = "NonMonotone") %>%
  filter(NonMonotone, Tech_Cost_Label == main_text_cost_scenario)

tornado_out <- tornado %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  arrange(Tech_Cost_Label, Resolution, Scenario_Label, desc(Range))
write_csv(tornado_out, file.path(csv_folder, "Table_sensitivity_ranking.csv"))

if (nrow(tornado_nonmono) > 0) {
  message("\nNON-MONOTONE dimensions in the tornado (baseline outside the ",
          "low-high interval; the bar spans all three settings):")
  print(as.data.frame(tornado_nonmono %>%
                        mutate(Dimension = gsub("\n", " ", Dimension)) %>%
                        select(Scenario_Label, Resolution, Dimension)),
        row.names = FALSE)
}


# =============================================================================
# FIGURES
# =============================================================================

message("Generating figures...")

# The euro sign is not in the default PDF device's font encoding, which produces
# "conversion failure in mbcsToSbcs" warnings and substitutes dots for the glyph.
# cairo_pdf handles UTF-8 correctly; fall back to pdf() if cairo is unavailable.
if (capabilities("cairo")) {
  cairo_pdf(pdf_output, width = 16, height = 10, onefile = TRUE)
} else {
  warning("cairo not available — the euro symbol may not render. ",
          "Install a cairo-enabled R build, or replace \u20ac with 'EUR'.")
  pdf(pdf_output, width = 16, height = 10)
}

# ADAPTIVE AXIS BREAKS
#
# A single fixed step does not suit every figure: the axis ranges here run from
# about 60 EUR/ha (Fig. 8) to about 500 EUR/ha (Fig. 4b), so a 20-unit step would
# give 3 labels on one and 25 on the other. This picks a "nice" step from a fixed
# set so that each panel ends up with roughly `target` labels, which keeps values
# readable off the axis without crowding.
#
# Fig. 5b (all five cost scenarios) is deliberately left on a coarse 100-unit
# step: with 15 panels there is no room for dense labelling.
nice_breaks <- function(target = 10) {
  function(lims) {
    lo <- min(lims, na.rm = TRUE)
    hi <- max(lims, na.rm = TRUE)
    span <- hi - lo
    if (!is.finite(span) || span <= 0) return(pretty(lims))
    cand <- c(1, 2, 2.5, 5, 10, 20, 25, 50, 100, 200, 250, 500)
    step <- cand[which.min(abs(span / cand - target))]
    seq(floor(lo / step) * step, ceiling(hi / step) * step, by = step)
  }
}

# =============================================================================
# FIGURE TEXT SIZING — ALL SIZES CONTROLLED FROM HERE
#
# Text legibility depends on the RATIO of font size to canvas size, not on the
# font size alone. These figures are drawn on a 16 x 10 inch canvas and then
# scaled to about 6.5 inches when placed on an A4 page, i.e. reduced by a factor
# of roughly 2.5. A size-14 label therefore prints at about 5.7 pt, which is why
# the multi-panel figures were hard to read.
#
# TO MAKE ALL TEXT LARGER OR SMALLER, CHANGE FIG_SCALE ONLY.
#   FIG_SCALE <- 1.0   sizes below, printing at roughly 9-11 pt on an A4 page
#   FIG_SCALE <- 1.2   about 20% larger
#   FIG_SCALE <- 0.85  about 15% smaller
#
# The dense multi-panel figures (Fig. 5b, Fig. 7) are additionally exported at a
# larger canvas at the end of the script, so their text stays legible without
# crowding the panels.
# =============================================================================
FIG_SCALE <- 1.0

sz <- function(x) round(x * FIG_SCALE)

base_theme <- theme_bw() +
  theme(
    legend.position    = "bottom",
    axis.title.x       = element_text(size = sz(26), face = "bold", margin = margin(t = 15)),
    axis.title.y       = element_text(size = sz(26), face = "bold", margin = margin(r = 15)),
    axis.text.y        = element_text(size = sz(24)),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = sz(22), face = "bold"),
    strip.text         = element_text(size = sz(22), face = "bold"),
    legend.title       = element_text(size = sz(24), face = "bold"),
    legend.text        = element_text(size = sz(24)),
    plot.caption       = element_text(size = sz(15), hjust = 0, margin = margin(t = 10)),
    panel.grid.minor.y = element_blank()
  )

# Previously about 0.55 of base_theme, which printed at roughly 5-6 pt.
# Now 0.85 of base, printing at roughly 8-9 pt. Panels are slightly tighter but
# every label is readable at final page size.
compact_theme <- base_theme +
  theme(axis.title.x = element_text(size = sz(22), face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = sz(22), face = "bold", margin = margin(r = 10)),
        axis.text.y  = element_text(size = sz(19)),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = sz(18), face = "bold"),
        strip.text   = element_text(size = sz(19), face = "bold"),
        legend.title = element_text(size = sz(20), face = "bold"),
        legend.text  = element_text(size = sz(19)),
        plot.caption = element_text(size = sz(14), hjust = 0))

supp_caption_theme <- base_theme +
  theme(plot.caption = element_text(size = sz(15), hjust = 0, margin = margin(t = 10)))


# --- Figure 2: Herbicide savings (mixed-model letters) ---

vol_dat <- dat %>% filter(Tech_Cost_Label == main_text_cost_scenario)

fig2 <- ggplot(vol_dat, aes(x = Resolution, y = Vol_Mean, fill = Scenario_Label)) +
  geom_bar(stat = "identity", position = position_dodge(0.9),
           color = "black", alpha = 0.9) +
  # clamp the lower bound at 0. Negative herbicide savings are not
  # meaningful, and scale_y_continuous(limits=) would otherwise DROP the whole
  # error bar rather than clip it, leaving a floating dash above the bar.
  geom_errorbar(aes(ymin = pmax(0, Vol_Mean - Vol_SD), ymax = Vol_Mean + Vol_SD),
                position = position_dodge(0.9), width = 0.25, linewidth = 0.7) +
  geom_text(aes(label = Vol_Group, y = 105),
            position = position_dodge(0.9), vjust = 0,
            size = 8, fontface = "bold", color = "black") +
  scale_fill_manual(values = pub_colors) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, by = 10)) +
  labs(y = "Herbicide savings (%)", x = "Resolution", fill = "Decision Strategy") +
  base_theme


# --- Figure 4: Net return difference — weed dynamics sensitivity ---

fig4 <- ggplot(dat, aes(x = Resolution, y = Econ_Mean, fill = Scenario_Label)) +
  geom_errorbar(aes(ymin = Econ_Weed_Pess, ymax = Econ_Weed_Opt),
                position = position_dodge(0.9), width = 0.5,
                linewidth = 1.1, color = "black") +
  geom_bar(stat = "identity", position = position_dodge(0.9),
           color = "black", alpha = 0.78) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  facet_wrap(~Tech_Cost_Label, ncol = 5) +
  scale_fill_manual(values = pub_colors) +
  scale_y_continuous(breaks = nice_breaks(10)) +   # readable axis
  labs(y = "Difference in Net Return (\u20ac/ha)",
       x = "Resolution", fill = "Decision Strategy") +
  base_theme


# --- Figure 4b: between-field distribution ---

fig4b <- ggplot(analysis_data,
                aes(x = Resolution, y = DIFF_Total, fill = Scenario_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_boxplot(position = position_dodge(0.9), outlier.size = 0.8,
               alpha = 0.78, color = "black", linewidth = 0.4) +
  facet_wrap(~Tech_Cost_Label, ncol = 5) +
  scale_fill_manual(values = pub_colors) +
  scale_y_continuous(breaks = nice_breaks(10)) +   # readable axis
  labs(y = "Difference in Net Return (\u20ac/ha)", x = "Resolution",
       fill = "Decision Strategy",
       caption = paste("Between-field distribution across the 19 simulated fields.",
                       "Boxes show the interquartile range, whiskers 1.5 x IQR.",
                       "The dashed line marks parity with broadcast application.")) +
  compact_theme


# --- Figure 5: partial budget decomposition ---

if (has_decomp) {
  
  # Main-text version: one cost scenario, faceted by resolution
  decomp_main <- decomp_long %>% filter(Tech_Cost_Label == main_text_cost_scenario)
  net_main    <- decomp       %>% filter(Tech_Cost_Label == main_text_cost_scenario)
  
  fig5 <- ggplot(decomp_main, aes(x = Scenario_Label, y = Value, fill = Component)) +
    geom_col(position = "stack", color = "black", linewidth = 0.3, width = 0.7) +
    geom_point(data = net_main, aes(x = Scenario_Label, y = Net),
               inherit.aes = FALSE, shape = 23, size = 5,
               fill = "white", color = "black", stroke = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) +
    facet_wrap(~Resolution, nrow = 1) +
    scale_fill_manual(values = component_colors) +
    scale_y_continuous(breaks = nice_breaks(10)) +   # readable axis
    labs(y = "Contribution to net return difference (\u20ac/ha)",
         x = "Decision strategy", fill = "Budget component",
         caption = paste0("Partial budget decomposition versus broadcast application (",
                          main_text_cost_scenario,
                          "). Positive segments are savings, negative segments costs. ",
                          "The white diamond marks the net difference.")) +
    compact_theme
  
  # Complete version: all five cost scenarios
  fig5_full <- ggplot(decomp_long, aes(x = Scenario_Label, y = Value, fill = Component)) +
    geom_col(position = "stack", color = "black", linewidth = 0.25, width = 0.7) +
    geom_point(data = decomp, aes(x = Scenario_Label, y = Net),
               inherit.aes = FALSE, shape = 23, size = 3,
               fill = "white", color = "black", stroke = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7) +
    facet_grid(Resolution ~ Tech_Cost_Label) +
    scale_fill_manual(values = component_colors) +
    # 15 panels, so keep the coarse 100-unit step but add a minor gridline
    # midway so intermediate values can still be read off.
    scale_y_continuous(breaks = seq(-1000, 1000, by = 100),
                       minor_breaks = seq(-1000, 1000, by = 50)) +
    # 15 panels: trim the x labels slightly so they do not collide. This
    # figure is also exported separately at a larger canvas (see end of script).
    theme(axis.text.x = element_text(size = sz(15))) +
    labs(y = "Contribution to net return difference (\u20ac/ha)",
         x = "Decision strategy", fill = "Budget component",
         caption = paste("Partial budget decomposition across all five technology",
                         "cost scenarios. The white diamond marks the net difference.")) +
    compact_theme
  
  
  # --- Figure 6: break-even technology cost ---
  
  be_plot_dat <- clean_data %>%
    select(Field_Name, Scenario_Label, Resolution, Tech_Breakeven)
  
  tech_ranges <- tech_cost_table %>%
    group_by(Resolution) %>%
    summarise(lo = min(Tech_Cost_EUR), hi = max(Tech_Cost_EUR), .groups = "drop")
  
  fig6 <- ggplot(be_plot_dat,
                 aes(x = Scenario_Label, y = Tech_Breakeven, fill = Scenario_Label)) +
    geom_rect(data = tech_ranges, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
              fill = "grey80", alpha = 0.45) +
    geom_boxplot(width = 0.55, alpha = 0.85, color = "black",
                 linewidth = 0.4, outlier.size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~Resolution, nrow = 1) +
    scale_fill_manual(values = pub_colors, guide = "none") +
    scale_y_continuous(breaks = nice_breaks(10)) +   # readable axis
    labs(y = "Break-even technology cost (\u20ac/ha)", x = "Decision strategy",
         caption = paste("Maximum annualised technology and information cost at which",
                         "SSWM matches broadcast application, by field.",
                         "The shaded band spans the five technology cost scenarios of Table 4;",
                         "a box below the band means no scenario is economically viable.")) +
    compact_theme
}


# --- Figure 7: tornado ---

tornado_main <- tornado %>% filter(Tech_Cost_Label == main_text_cost_scenario)

fig7 <- ggplot(tornado_main,
               aes(x = Range, y = reorder(Dimension, Range))) +
  geom_col(fill = "grey55", color = "black", width = 0.65, linewidth = 0.3) +
  facet_grid(Resolution ~ Scenario_Label) +
  scale_x_continuous(breaks = nice_breaks(6)) +   # readable axis
  labs(x = "Range in net return difference (\u20ac/ha)", y = NULL,
       caption = paste0("Sensitivity ranking (", main_text_cost_scenario,
                        "). Bars show the spread in mean net return produced by each ",
                        "sensitivity dimension. Longer bars indicate greater influence.")) +
  compact_theme +
  theme(axis.text.y = element_text(size = sz(17), angle = 0, hjust = 1, face = "plain",
                                   lineheight = 0.9),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = sz(17), face = "plain"))


# --- Figure 8: dual-dose background efficacy ---

if (!is.null(dd_sweep) && nrow(dd_sweep) > 0) {
  dd_plot_dat <- dd_sweep %>%
    left_join(tech_cost_table, by = "Resolution", relationship = "many-to-many") %>%
    filter(Tech_Cost_Label == main_text_cost_scenario) %>%
    mutate(DIFF = Net - Tech_Cost_EUR) %>%
    group_by(Resolution, K_Background) %>%
    summarise(Mean_Net = mean(DIFF, na.rm = TRUE),
              SD_Net   = sd(DIFF,   na.rm = TRUE), .groups = "drop")
  
  fig8 <- ggplot(dd_plot_dat, aes(x = K_Background, y = Mean_Net)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
    # mark the efficacy at which the main results are reported, so the
    # reader can see where the rest of the paper sits on this curve.
    {if (!is.na(DD_SIMULATED_K))
      geom_vline(xintercept = DD_SIMULATED_K, linetype = "dotted",
                 linewidth = 0.7, colour = "grey40")} +
    geom_errorbar(aes(ymin = Mean_Net - SD_Net, ymax = Mean_Net + SD_Net),
                  width = 0.012, linewidth = 0.6) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3.5, shape = 21, fill = "white", stroke = 1.1) +
    facet_wrap(~Resolution, nrow = 1) +
    scale_x_reverse(breaks = c(0.95, 0.90, 0.85, 0.80, 0.70)) +  # label all
    scale_y_continuous(breaks = nice_breaks(10)) +               # readable
    
    labs(x = "Efficacy of the 50% background dose (kill fraction)",
         y = "Difference in Net Return (\u20ac/ha)",
         caption = paste("Dual-dose strategy under reduced background-dose efficacy.",
                         "At a kill fraction of 0.95 the background dose controls weeds as",
                         "well as the full dose. The main results are reported at 0.85,",
                         "the mid-range of published reduced-rate efficacies.")) +
    compact_theme
}


# --- Figure 3 (new): break-even technology cost as a function of recall ---
#
# Replaces the boxplot (kept below as fig6 for the supplement). The curve is
# the mean across fields of Net(r) = herbicide saved + year-0 yield difference
# + FWP difference, i.e. the maximum technology cost at which SSWM matches
# broadcast if the sensor had recall r. The shaded band is the five cost
# scenarios (Table 3): where a curve enters the band from below is the
# break-even recall for the cheapest scenario, where it leaves it the dearest.
# The vertical dotted line is the measured recall, so the old Fig. 3 is the
# slice of this figure at that line.

line_colors <- c("ET strategy" = "#8c8c8c", "PB strategy" = "#000000",
                 "DD strategy" = "#4a4a4a", "DynET strategy" = "#b0b0b0")
line_types  <- c("ET strategy" = "dashed", "PB strategy" = "solid",
                 "DD strategy" = "dotdash", "DynET strategy" = "dotted")

if (use_curve) {
  # cost lines: S2 (main-text scenario) dashed, the others thin grey;
  # scenario pairs with (almost) equal costs share one label.
  cost_lines <- tech_cost_table %>%
    mutate(S = paste0("S", str_extract(Tech_Cost_Label, "[0-9]+")),
           Main = Tech_Cost_Label == main_text_cost_scenario)
  cost_labels <- cost_lines %>%
    mutate(Group = case_when(S %in% c("S1", "S4") ~ "S1/S4",
                             S %in% c("S3", "S5") ~ "S3/S5",
                             TRUE ~ S)) %>%
    group_by(Resolution, Group) %>%
    summarise(y = mean(Tech_Cost_EUR), Main = any(Main), .groups = "drop")
  s2_cost <- cost_lines %>% filter(Main) %>% select(Resolution, S2_Cost = Tech_Cost_EUR)

  x_lo <- CURVE_X_MIN - 0.01
  x_hi <- CURVE_X_MAX + 0.055          # room for the cost labels

  curve_zero <- curve_mean %>%
    filter(Tech_Cost_Label == "Zero technology cost",
           Scenario_Label %in% CURVE_STRATEGIES,
           Recall >= CURVE_X_MIN - 1e-9, Recall <= CURVE_X_MAX + 1e-9) %>%
    mutate(Band_Lo = case_when(CURVE_BAND == "iqr" ~ Q25_DIFF,
                               CURVE_BAND == "sd"  ~ Mean_DIFF - SD_DIFF,
                               TRUE ~ NA_real_),
           Band_Hi = case_when(CURVE_BAND == "iqr" ~ Q75_DIFF,
                               CURVE_BAND == "sd"  ~ Mean_DIFF + SD_DIFF,
                               TRUE ~ NA_real_))

  # genuine crossings inside the plotted range only
  be_zero <- be_mean %>%
    filter(Tech_Cost_Label == "Zero technology cost", !is.na(BE_Recall_MeanCurve),
           Scenario_Label %in% CURVE_STRATEGIES,
           BE_Recall_MeanCurve > CURVE_X_MIN + 1e-9,
           BE_Recall_MeanCurve <= CURVE_X_MAX + 1e-9) %>%
    mutate(y = 0)
  be_s2 <- be_mean %>%
    filter(Tech_Cost_Label == main_text_cost_scenario, !is.na(BE_Recall_MeanCurve),
           Scenario_Label %in% CURVE_STRATEGIES,
           BE_Recall_MeanCurve > CURVE_X_MIN + 1e-9,
           BE_Recall_MeanCurve <= CURVE_X_MAX + 1e-9) %>%
    mutate(y = Tech_Cost_EUR)

  band_caption <- switch(CURVE_BAND,
                         iqr  = "Shaded ribbons: interquartile range across the 19 fields.",
                         sd   = "Shaded ribbons: \u00b11 SD across the 19 fields.",
                         none = "")

  fig3_curve <- ggplot(curve_zero, aes(x = Recall, y = Mean_DIFF,
                                       colour = Scenario_Label,
                                       linetype = Scenario_Label)) +
    geom_hline(data = cost_lines %>% filter(!Main),
               aes(yintercept = Tech_Cost_EUR), colour = "grey75", linewidth = 0.5) +
    geom_hline(data = cost_lines %>% filter(Main),
               aes(yintercept = Tech_Cost_EUR), colour = "grey35",
               linetype = "longdash", linewidth = 0.7) +
    geom_text(data = cost_labels, inherit.aes = FALSE,
              aes(x = CURVE_X_MAX + 0.008, y = y, label = Group),
              hjust = 0, size = sz(4.2), colour = "grey35") +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.6) +
    geom_vline(xintercept = r_base, linetype = "dotted", linewidth = 0.7,
               colour = "grey40") +
    {if (CURVE_BAND != "none")
      geom_ribbon(aes(ymin = Band_Lo, ymax = Band_Hi, fill = Scenario_Label),
                  alpha = 0.12, colour = NA, show.legend = FALSE)} +
    geom_line(linewidth = 1.2) +
    geom_point(data = be_zero, inherit.aes = FALSE,
               aes(x = BE_Recall_MeanCurve, y = y), shape = 23, size = 4,
               fill = "white", stroke = 1.1) +
    geom_point(data = be_s2, inherit.aes = FALSE,
               aes(x = BE_Recall_MeanCurve, y = y), shape = 23, size = 4,
               fill = "black", stroke = 1.1) +
    facet_wrap(~Resolution, nrow = 1) +
    scale_colour_manual(values = line_colors, name = "Decision strategy") +
    scale_linetype_manual(values = line_types, name = "Decision strategy") +
    scale_fill_manual(values = line_colors) +
    scale_x_continuous(breaks = CURVE_X_BREAKS,
                       labels = format(CURVE_X_BREAKS, drop0trailing = TRUE, trim = TRUE)) +
    coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
    scale_y_continuous(breaks = nice_breaks(10)) +
    labs(x = "Sensor recall (proportion of weeds detected)",
         y = "Break-even technology cost (\u20ac/ha)",
         caption = str_wrap(paste(
           "Maximum annualised technology and information cost at which SSWM",
           "matches broadcast application, as a function of the recall of the",
           "weed-detection system applied to the same weed population.",
           "Lines: mean of 19 fields.", band_caption,
           "Horizontal lines: annualised cost of the five technology scenarios",
           "(Table 3; dashed: S2, least expensive). Dotted line: measured recall (0.41).",
           "Open diamonds: break-even recall with free technology; filled diamonds:",
           "break-even recall under S2. ET and PB coincide at 25 x 25 cm.",
           sprintf("Curves end at recall %.2f; at recall 1.0 no weed escapes and the curve jumps.",
                   CURVE_X_MAX)), width = 135)) +
    compact_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          panel.spacing = unit(2.2, "lines"))

  # --- Supplementary: field-level curves, strategy x resolution ---
  curve_fields <- recall_curve %>%
    filter(Recall >= CURVE_X_MIN - 1e-9, Recall <= CURVE_X_MAX + 1e-9)
  curve_zero_mean <- curve_zero %>% select(Scenario_Label, Resolution, Recall, Mean_DIFF)
  be_s2_all <- be_mean %>%
    filter(Tech_Cost_Label == main_text_cost_scenario, !is.na(BE_Recall_MeanCurve),
           BE_Recall_MeanCurve > CURVE_X_MIN + 1e-9,
           BE_Recall_MeanCurve <= CURVE_X_MAX + 1e-9) %>%
    mutate(y = Tech_Cost_EUR)

  figS_curve_all <- ggplot(curve_fields, aes(x = Recall, y = Net)) +
    geom_hline(data = cost_lines %>% filter(!Main),
               aes(yintercept = Tech_Cost_EUR), colour = "grey80", linewidth = 0.4) +
    geom_hline(data = cost_lines %>% filter(Main),
               aes(yintercept = Tech_Cost_EUR), colour = "grey35",
               linetype = "longdash", linewidth = 0.6) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    geom_vline(xintercept = r_base, linetype = "dotted", linewidth = 0.6,
               colour = "grey40") +
    geom_line(aes(group = Field_Name), colour = "grey55", linewidth = 0.45) +
    geom_line(data = curve_zero_mean, aes(y = Mean_DIFF), colour = "black",
              linewidth = 1.2) +
    geom_point(data = be_s2_all, inherit.aes = FALSE,
               aes(x = BE_Recall_MeanCurve, y = y), shape = 23, size = 3.5,
               fill = "black", stroke = 1.0) +
    facet_grid(Scenario_Label ~ Resolution, scales = "free_y") +
    geom_text(data = cost_labels, inherit.aes = FALSE,
              aes(x = CURVE_X_MAX + 0.008, y = y, label = Group),
              hjust = 0, size = sz(3.8), colour = "grey35") +
    scale_x_continuous(breaks = CURVE_X_BREAKS,
                       labels = format(CURVE_X_BREAKS, drop0trailing = TRUE, trim = TRUE)) +
    coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
    scale_y_continuous(breaks = nice_breaks(6)) +
    labs(x = "Sensor recall (proportion of weeds detected)",
         y = "Break-even technology cost (\u20ac/ha)",
         caption = str_wrap(paste(
           "Break-even technology cost as a function of sensor recall for each",
           "field (grey) and the mean of 19 fields (black). Horizontal lines:",
           "technology cost scenarios (dashed: S2). Filled diamonds: recall at",
           "which the mean covers S2. Dotted line: measured recall (0.41).",
           "Break-even recalls for all five scenarios are given in Supplementary Table [x]."),
           width = 135)) +
    compact_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = sz(15)),
          strip.text = element_text(size = sz(15)),
          panel.spacing = unit(1.2, "lines"))
}


# --- Supplementary sensitivity figures (same layout as Figure 4) ---

# DD error bars are drawn as for the other strategies; they are collapsed onto
# the mean only if the DD rows were substituted from the sweep.
supp_fig_dat <- dat
if (dd_substituted) {
  supp_fig_dat <- supp_fig_dat %>%
    mutate(across(c(Econ_Weed_Lo, Econ_Weed_Hi,
                    Econ_Split_Lo, Econ_Split_Hi,
                    Econ_Kill_Lo,  Econ_Kill_Hi),
                  ~ ifelse(Scenario_Label == "DD strategy", Econ_Mean, .x)))
  message("Figures S4-S6: dual-dose bars shown without error bars because ",
          "DD_REPORT_K differs from the simulated background efficacy.")
}

# report any non-monotone cases so they are explained rather than discovered
nm <- dat %>%
  filter(Tech_Cost_Label == main_text_cost_scenario,
         NonMono_Weed | NonMono_Split) %>%
  select(Scenario_Label, Resolution, NonMono_Weed, NonMono_Split)
if (nrow(nm) > 0) {
  message("\nNON-MONOTONE sensitivity responses (baseline outside the low-high ",
          "interval):")
  print(as.data.frame(nm), row.names = FALSE)
  message("  The error bars in Figures S4-S6 span the full range across all ",
          "settings\n  INCLUDING the baseline, so the bar always sits inside ",
          "its interval.")
}

# note appended to every supplementary sensitivity caption when DD is excluded
dd_note <- if (dd_substituted)
  paste("The dual-dose strategy is shown without an interval because its",
        "baseline and its sensitivity columns come from different background",
        "efficacies.") else
          paste("For the dual-dose strategy the herbicide-efficacy interval varies the",
                "full-rate kill fraction while the background dose is held at its own",
                "value, so it measures sensitivity to the difference between the two",
                "doses; sensitivity to the background dose itself is shown in Fig. 8.")

make_supp_fig <- function(dat, ymin_col, ymax_col, caption_text) {
  caption_text <- paste(caption_text, dd_note)
  ggplot(dat, aes(x = Resolution, y = Econ_Mean,
                  fill = Scenario_Label, group = Scenario_Label)) +
    geom_bar(stat = "identity", position = position_dodge(0.9),
             color = "black", alpha = 0.78) +
    geom_errorbar(aes(ymin = .data[[ymin_col]], ymax = .data[[ymax_col]],
                      group = Scenario_Label),
                  position = position_dodge(0.9),
                  width = 0.5, linewidth = 1.2, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
    facet_wrap(~Tech_Cost_Label, ncol = 5) +
    scale_fill_manual(values = pub_colors) +
    scale_y_continuous(breaks = nice_breaks(10)) +   # readable axis
    labs(y = "Difference in Net Return (\u20ac/ha)", x = "Resolution",
         fill = "Decision Strategy", caption = caption_text) +
    supp_caption_theme
}

figS4 <- make_supp_fig(supp_fig_dat, "Econ_Weed_Lo", "Econ_Weed_Hi",
                       paste("Sensitivity to weed population dynamics.",
                             "Error bars span the full range of outcomes across a",
                             "\u00b120% variation in the population parameters (Table 3),",
                             "including the baseline. Note that REDUCING the parameters",
                             "widens the gap to broadcast application, because broadcast",
                             "then suppresses the population more effectively while",
                             "undetected weeds continue to escape."))

figS5 <- make_supp_fig(supp_fig_dat, "Econ_Split_Lo", "Econ_Split_Hi",
                       paste("Sensitivity to false-negative spatial allocation.",
                             "Error bars show the net return range between 50% and 90% of",
                             "undetected weeds allocated to detected-weed patches",
                             "(baseline: 70/30 allocation)."))

# caption updated for the RECALIBRATED seed-set values. The model now
# uses gamma = 0.08 / 0.13 / 0.16 (rebased around the calibrated baseline of
# 0.13), not the former 0.20 / 0.33 / 0.40.
figS6 <- make_supp_fig(supp_fig_dat, "Econ_Kill_Lo", "Econ_Kill_Hi",
                       paste("Sensitivity to herbicide efficacy. Error bars span",
                             "k = 0.85 with a seed-set factor of 0.08 and k = 0.99",
                             "with 0.16; the baseline is k = 0.95 with 0.13.",
                             "Note that this dimension varies the kill rate and the",
                             "seed set of survivors together."))

# All four supplementary panels now show the same three strategies.
# The dual-dose recall interval is genuinely zero-width: the strategy treats the
# whole field regardless of what the sensor reports.
figS7 <- make_supp_fig(
  supp_fig_dat, "Econ_Recall_Lo", "Econ_Recall_Hi",
  if (use_curve)
    sprintf(paste("Sensitivity to sensor recall (\u00b120%%). Error bars span a sensor of",
                  "recall %.3f (false-negative rate %.0f%%) and %.3f (false-negative rate %.0f%%)",
                  "applied to the same weed population; baseline recall 0.41 (59%%).",
                  "The year-0 spray map, the year-0 yield loss and the future weed penalty",
                  "all vary with recall."),
            r_lo, 100 * (1 - r_lo), r_hi, 100 * (1 - r_hi))
  else if (use_new_recall)
    paste("Sensitivity to sensor recall (\u00b120%). Validation design: spray map held fixed,",
          "true population varied. Error bars span recall 0.328 and 0.492.")
  else
    paste("Sensitivity to sensor recall (\u00b120%). WARNING: computed from the",
          "liability-only columns, which hold Year-0 revenue fixed and therefore",
          "understate the influence of recall. Re-run the simulation model."))


# --- Print all figures ---
print(fig2)
print(fig4)
print(fig4b)
if (has_decomp) { print(fig5); print(fig5_full); print(fig6) }
if (use_curve) { print(fig3_curve); if (MAKE_FIELD_CURVE_FIGURE) print(figS_curve_all) }
print(fig7)
if (!is.null(dd_sweep) && nrow(dd_sweep) > 0) print(fig8)
print(figS4)
print(figS5)
print(figS6)
print(figS7)

dev.off()

# The 15-panel decomposition figure (all cost scenarios) is exported
# separately, larger and in landscape. The 3-panel version (fig5) is the
# main-text figure.
if (has_decomp) {
  big_path <- file.path(csv_folder, "FigS_decomposition_all_scenarios_LARGE.pdf")
  tryCatch({
    ggsave(big_path, fig5_full, width = 20, height = 12, units = "in",
           device = cairo_pdf, limitsize = FALSE)
    message("Large supplementary decomposition figure: ", big_path)
  }, error = function(e) {
    ggsave(big_path, fig5_full, width = 20, height = 12, units = "in",
           limitsize = FALSE)
    message("Large supplementary decomposition figure (default device): ", big_path)
  })
  
  # recall-curve figures, exported on their own as well
  if (use_curve) {
    tryCatch({
      ggsave(file.path(csv_folder, "Fig3_breakeven_vs_recall.pdf"), fig3_curve,
             width = 16, height = 8, units = "in", device = cairo_pdf, limitsize = FALSE)
      ggsave(file.path(csv_folder, "Fig3_breakeven_vs_recall.png"), fig3_curve,
             width = 16, height = 8, units = "in", dpi = 300, limitsize = FALSE)
      if (MAKE_FIELD_CURVE_FIGURE) {
        ggsave(file.path(csv_folder, "FigS_breakeven_vs_recall_by_field.pdf"),
               figS_curve_all, width = 16, height = 14, units = "in",
               device = cairo_pdf, limitsize = FALSE)
        ggsave(file.path(csv_folder, "FigS_breakeven_vs_recall_by_field.png"),
               figS_curve_all, width = 16, height = 14, units = "in", dpi = 300,
               limitsize = FALSE)
      }
      message("Recall-curve figures exported (Fig3_breakeven_vs_recall.*, ",
              "FigS_breakeven_vs_recall_by_field.*)")
    }, error = function(e) message("Recall-curve figure export skipped: ",
                                    conditionMessage(e)))
  }

  # the tornado also has 9 panels with long labels; export it larger too
  tor_path <- file.path(csv_folder, "Fig7_sensitivity_tornado_LARGE.pdf")
  tryCatch({
    ggsave(tor_path, fig7, width = 18, height = 11, units = "in",
           device = cairo_pdf, limitsize = FALSE)
    message("Large tornado figure: ", tor_path)
  }, error = function(e) message("Tornado export skipped: ", conditionMessage(e)))
  
  # 300 dpi raster version, in case the journal wants an image rather than vector
  png_path <- file.path(csv_folder, "FigS_decomposition_all_scenarios_LARGE.png")
  tryCatch({
    ggsave(png_path, fig5_full, width = 20, height = 12, units = "in",
           dpi = 300, limitsize = FALSE)
    message("                                    (png): ", png_path)
  }, error = function(e) message("PNG export skipped: ", conditionMessage(e)))
}

message("\n", strrep("=", 70))
message("Done. Figures saved to: ", pdf_output)
message("Tables saved to: ", csv_folder)
message(strrep("=", 70))
message("Figures in the PDF (manuscript numbering):")
message("  Fig. 1  herbicide savings with mixed-model letters (fig2)")
message("  additional: net return with weed dynamics sensitivity (fig4)")
message("  additional: between-field distribution (fig4b)")
if (has_decomp) {
  message("  Fig. 2  partial budget decomposition, cost scenario 2 (fig5)")
  message("  Supplementary: decomposition, all cost scenarios (fig5_full)")
  message("  Supplementary: break-even technology cost by field (fig6)")
}
if (use_curve) {
  message("  Fig. 3  break-even technology cost vs recall (fig3_curve)")
  if (MAKE_FIELD_CURVE_FIGURE)
    message("  additional: break-even technology cost vs recall, by field")
}
message("  Fig. 4  sensitivity comparison (fig7)")
if (!is.null(dd_sweep) && nrow(dd_sweep) > 0)
  message("  Fig. 5  efficacy of the background dose (fig8)")
message("  Supplementary: sensitivity to population dynamics, false-negative ",
        "allocation, herbicide efficacy and recall (figS4-figS7)")
