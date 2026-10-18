# =============================================================================
# SSWM bio-economic simulation model                      (Supplementary Code 1)
# Version 2.0
#
# Economic evaluation of site-specific weed management under imperfect sensing:
# a multi-year bio-economic simulation study.
#
# WHAT THE SCRIPT DOES
#   For one field at one grid resolution, the script reads the UAV weed
#   detections per grid cell (dicots and monocots), estimates the true weed
#   population (Eq. 7), applies one decision strategy in year 0, projects the
#   weed population and its costs over ten years with probabilistic detection
#   (Eqs. 8-12), and compares the result with broadcast application in a
#   partial budget (Eq. 6). It also computes the sensitivity analyses, the
#   break-even recall curve and the marginal shadow cost of one additional weed
#   in every grid cell. Per-cell results are written to the grid layer and,
#   optionally, to CSV; field-level results are appended to CSV.
#
# STRATEGIES (SCENARIO_ID)
#   1  ET                  economic threshold on the detected weed counts
#   2  PB                  presence based: every cell with a detection is sprayed
#   3  DD                  dual dose: reduced background dose over the whole
#                          field, full dose in cells with a detection
#   4  forward-looking ET  threshold on the recall-corrected density, including
#                          the discounted future weed penalty
#
# REQUIREMENTS
#   QGIS 3.x (tested with 3.36) and its Python console. No extra packages.
#
# INPUT
#   One QGIS project per field and grid resolution. The project holds a polygon
#   grid layer (GRID_LAYER_NAME) whose cells contain the number of detected
#   dicots and monocots (INPUT_DICOT_FIELD, INPUT_MONOCOT_FIELD), as produced
#   by the grid script (Supplementary Code 3). Edge cells are clipped to the
#   field boundary; the model uses each cell's actual area. The project file
#   name identifies field and resolution in all outputs
#   (e.g. 'R101_econ_model', 'R101_econ_model-50', 'R101_econ_model-25').
#
# HOW TO RUN
#   1. Set RESULTS_CSV_PATH below. All runs append to the same files.
#   2. Open the QGIS project of one field and resolution.
#   3. Set SCENARIO_ID and run this script in the QGIS Python console.
#   4. Repeat for every strategy, resolution and field
#      (19 fields x 3 resolutions x 4 strategies = 228 runs in the study).
#   Running a field twice appends a second row; start from empty files.
#
# OUTPUT
#   Grid layer      per-cell results: true counts and densities, spray decision,
#                   dose, revenue, herbicide cost, future liability, net return
#                   and the marginal shadow costs Marg_SC_D / Marg_SC_M.
#   RESULTS_CSV_PATH               one row per run: partial budget and
#                                  sensitivity analyses.
#   <results>_recall_curve.csv     break-even recall analysis (ET, PB, DD).
#   <results>_dd_kill_sweep.csv    efficacy of the background dose (DD runs).
#   <results>_yearly.csv           undiscounted costs in years 1-10.
#   cells/<project>_<strategy>.csv per-cell results for the marginal
#                                  shadow-cost figure (by default ET at 1 x 1 m).
#   The CSV files are the input of the R scripts (Supplementary Codes 2 and 4).
#
# All reported results use the default settings below. Settings described as
# "validation" reproduce simplified variants of the model and were used only to
# check the implementation; no reported result uses them.
#
# The calculation core does not depend on QGIS; only the layer input and output
# under MAIN EXECUTION does.
# =============================================================================
import os
import csv
import math
from qgis.core import QgsProject, QgsField
from qgis.PyQt.QtCore import QVariant
from qgis.utils import iface

MODEL_VERSION = '2.0'

# =============================================================================
# RUN SETTINGS
# =============================================================================
SCENARIO_ID = 1
scenario_lookup = {
    1: 'SCENARIO_1_ECON',    # ET  - economic threshold on detected counts
    2: 'SCENARIO_2_SPOT',    # PB  - presence based
    3: 'SCENARIO_3_SECURE',  # DD  - dual dose
    4: 'SCENARIO_4_DYNET',   # forward-looking ET
}
SCENARIO = scenario_lookup[SCENARIO_ID]

# Technology and information costs (Table 3) are added in the R analysis
# (Supplementary Code 2); the model reports all results before them.
TECH_SCENARIO = 'BASED'

# Main results file. The other output files are written next to it.
RESULTS_CSV_PATH = r'C:\SSWM_simulation\results\results_thesis.csv'

# ----- OUTPUTS AND OPTIONAL ANALYSES -----
USE_CELL_CACHE        = True    # cells with identical inputs are computed once
ENABLE_RECALL_CURVE   = True    # break-even recall analysis
ENABLE_DD_KILL_SWEEP  = True    # efficacy of the background dose (DD runs only)
ENABLE_YEARLY_EXPORT  = True    # per-year costs (other discount rates / horizons)
ENABLE_CELL_EXPORT    = True    # per-cell CSV for the marginal shadow-cost figure
CELL_EXPORT_SCENARIOS = ('SCENARIO_1_ECON',)   # strategies exported per cell
CELL_EXPORT_CELL_SIZE_M = 1.0   # export only grids of this cell size (m); None = all
ENABLE_RECALL_SWEEP   = False   # validation: recall varied with the year-0 map fixed
ANNUALISE_SHADOW_COST = False   # False: marginal shadow cost as ten-year NPV

# =============================================================================
# MODEL SETTINGS
# =============================================================================

# Weeds present in a cell in years 1-10, given the detection outcome.
#   'thinned'  conditional on what the sensor detected (Poisson thinning): a
#              cell that is not treated holds its undetected weeds,
#              lambda (1 - r) + E[X | X < m]; a treated cell E[N | X >= m].
#              DEFAULT.
#   'v16'      validation: all lambda weeds survive when a cell is missed.
SURVIVOR_CONDITIONING = 'thinned'

# Allocation of the isolated share (1 - PATCH_FRACTION) of undetected weeds
# (Eq. 7).
#   'area'         uniform background density: each cell receives a share
#                  proportional to its area, whether or not it holds a
#                  detection. The estimated true weed map is then independent
#                  of the grid resolution. DEFAULT.
#   'empty_cells'  validation: divided equally among zero-detection cells.
ALLOCATION_RULE = 'area'

# Marginal shadow cost of one additional weed (Marg_SC_D / Marg_SC_M).
#   'fixed_management'  the weed is added to the year-0 population before the
#                       herbicide application; the base cell and the cell with
#                       the added weed follow the same future treatment events
#                       (those of the base cell), and the dormant seedbank of
#                       Eq. 8 is left unchanged. The value is the weed's own
#                       damage and is never negative. DEFAULT.
#   'adaptive'          validation: two separate simulations, so future
#                       treatment may respond to the added weed.
SHADOW_COST_MODE = 'fixed_management'

# Herbicide cost of the dual-dose strategy in years 1-10.
#   'dose_scaled'  background-dose years are charged DD_BASE_DOSE x the full
#                  herbicide cost, full-dose years the full cost. DEFAULT.
#   'v13'          validation: the full cost in every year. Together with
#                  DD_BASE_KILL = 0.95 the DD strategy is then identical to
#                  broadcast application, which checks the implementation.
DD_FUTURE_COST_MODE = 'dose_scaled'

# Decision density of the threshold strategies in years 1-10.
#   'consistent'  the same information as in year 0: detected counts for the ET
#                 (see ET_INFORMATION) and the recall-corrected density for the
#                 forward-looking ET. DEFAULT.
#   'v13'         validation: the detected density for both.
# PB and DD do not use a threshold and are not affected.
ET_FUTURE_RULE = 'consistent'

# Information the ET acts on, in year 0 and in years 1-10.
#   'detected'  the weed counts the sensor reports. The recall of a commercial
#               detection system is unknown to the operator and cannot be
#               measured in the field, so it cannot enter the decision rule.
#               DEFAULT.
#   'true'      validation: detected counts corrected for the known recall.
# The forward-looking ET always uses the recall-corrected density and the
# discounted future weed penalty. At 50 x 50 cm and 25 x 25 cm a single
# detected weed already exceeds the threshold under either setting.
ET_INFORMATION = 'detected'

# Detection in years 1-10.
#   'probabilistic'  the number of detected weeds in a cell is Poisson(lambda r),
#                    so sparsely infested cells can be missed in successive
#                    years and their seedbank compounds. DEFAULT.
#   'expected'       validation: detected weeds = emerged weeds x recall, so
#                    every cell holding weeds is found from year 1 onward.
# Year 0 always uses the observed UAV detections.
FUTURE_DETECTION = 'probabilistic'

# Branching control (FUTURE_DETECTION = 'probabilistic'). Each year a cell is
# either treated or not, so ten years give up to 2^10 = 1024 trajectories.
# Paths with less probability than PRUNE_PROB are dropped and paths whose
# seedbanks have converged are merged. The expected liability was identical to
# four decimal places for MAX_STATES between 128 and 4096 and MERGE_DIGITS
# between 3 and 5.
PRUNE_PROB   = 1e-6
MERGE_DIGITS = 3      # significant digits used when merging similar seedbanks
MAX_STATES   = 1024   # ceiling on simultaneous trajectories

# Share of undetected weeds allocated to cells with detections (Eq. 7); the
# remainder is the isolated share (ALLOCATION_RULE).
PATCH_FRACTION = 0.70

tech_cost_lookup = {
    'BASED':      0,
    'SCENARIO 1': 0,
    'SCENARIO 2': 0,
    'SCENARIO 3': 0,
    'SCENARIO 4': 0,
}
sswm_annuity_per_ha = tech_cost_lookup[TECH_SCENARIO]

GRID_LAYER_NAME     = 'grid_rotated'
INPUT_DICOT_FIELD   = 'DicotCount'
INPUT_MONOCOT_FIELD = 'MonoCount'

# ----- COMPETITION, YIELD LOSS AND PRICES (Table 1, Eqs. 1-5) -----
dicot_competition_index   = 8.0
monocot_competition_index = 2.0
weed_height_modifier      = 0.75    # competitive index modifier (CIM)

# Rectangular hyperbola on the combined competitive load (Eq. 3). The effective
# initial slope per plant is 0.5 x 8 x 0.75 = 3.00 for dicots and
# 0.5 x 2 x 0.75 = 0.75 for monocots; one asymptote applies to the combined load.
yield_loss_initial_slope = 0.50
yield_loss_asymptote_pct = 60.0

weed_free_yield_kg_per_m2          = 1.015
crop_price_eur_per_kg              = 0.229
herbicide_cost_fullrate_eur_per_m2 = 0.0130
full_kill_rate                     = 0.95    # herbicide efficacy k (Eq. 5)

# ----- DUAL-DOSE STRATEGY -----
DD_BASE_DOSE = 0.5     # background dose, share of the full rate
# Efficacy (kill fraction) of the background dose, in year 0 and in years 1-10.
# Full-dose cells always use full_kill_rate. Varied in DD_KILL_SWEEP.
DD_BASE_KILL = 0.85
security_base_dose = DD_BASE_DOSE
security_base_kill = DD_BASE_KILL

# ----- DETECTOR -----
recall_dicot   = 0.41    # measured recall (proportion of weeds detected)
recall_monocot = 0.41

# ================= BIOLOGICAL ENGINE =================
DISCOUNT_RATE        = 0.04
SIM_YEARS            = 10
MAX_HERB_COST_YR     = 130.0     # full-rate herbicide cost, EUR/ha

# Seed set of herbicide survivors (gamma_s, Eq. 9). Calibrated: projected under
# continuous broadcast application at k = 0.95, the management the fields have
# received, the population settles at 17.3 plants m-2 (observed median
# 16.4 plants m-2). 0.16 gives 33 and 0.11 gives 7 plants m-2; below about 0.10
# the population dies out. The check is printed at the start of every run.
STUNTED_SEED_PENALTY = 0.13

# Competitive ability of herbicide survivors relative to untreated weeds
# (gamma_c, Eq. 5).
SURVIVOR_COMPETITIVE_FACTOR = 1.0 / 3.0

# Annuity factor for r = 0.04 and T = 10 years (Eq. 12). The exact value
# (1 - 1.04^-10) / 0.04 = 8.110896 is used rounded, as in all reported results.
ANNUITY_FACTOR = 8.11

# Allocation sensitivity: share of undetected weeds allocated to cells with
# detections.
SPLIT_SCENARIOS = {
    'split_90_10': 0.90,
    'split_70_30': 0.70,
    'split_50_50': 0.50,
}

# Herbicide efficacy sensitivity: the full-rate kill fraction and the seed set
# of survivors are varied together.
KILL_SCENARIOS = {
    'k_low':  {'kill': 0.85, 'stunt': 0.08},
    'k_base': {'kill': 0.95, 'stunt': 0.13},
    'k_high': {'kill': 0.99, 'stunt': 0.16},
}

# Recall +/-20%. In the main results file these vary the future liability only;
# the reported recall sensitivity is taken from the recall curve at the same
# values (0.328 and 0.492), which also re-derives the year-0 spray map.
RECALL_SCENARIOS = {
    'recall_low':  round(recall_dicot * 0.8, 4),
    'recall_base': recall_dicot,
    'recall_high': round(recall_dicot * 1.2, 4),
}

# Validation (ENABLE_RECALL_SWEEP): recall varied with the year-0 spray map held
# at the observed detections.
RECALL_SWEEP = [0.328, 0.41, 0.492]

# ----- RECALL CURVE -----
# A hypothetical sensor of recall r' is applied to the same estimated true
# population, in year 0 and in years 1-10; broadcast application does not
# depend on the sensor. The year-0 spray map is re-derived for each r':
#   r' > r : the sensor finds everything the measured sensor found plus each
#            missed weed with probability q = (r' - r) / (1 - r) (nested sensor)
#   r' < r : each observed detection is kept with probability r' / r.
# The list must contain the measured recall (self-check against the baseline);
# 0.328 and 0.492 are the +/-20% values. The curve is steep just below r' = 1,
# hence the denser spacing there.
RECALL_CURVE = [0.328, 0.41, 0.45, 0.492, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,
                0.85, 0.90, 0.93, 0.95, 0.97, 0.98, 0.99, 1.00]
# Probability mass below which year-0 enumeration branches are dropped.
CURVE_PRUNE_PROB = 1e-7

# ----- EFFICACY OF THE BACKGROUND DOSE (DD runs) -----
# The background kill fraction acts in year 0 (undetected weeds in background-
# dose cells survive at 1 - k_background, causing yield loss and seed return)
# and in years 1-10 (background-dose branches). With DD_FUTURE_COST_MODE =
# 'dose_scaled', the 0.95 row keeps the future chemical saving of the
# background dose.
DD_KILL_SWEEP = [0.95, 0.90, 0.85, 0.80, 0.70]

# ----- WEED POPULATION DYNAMICS (Table 2) -----
#   max_seeds            maximum seed production per plant (S_max, Eq. 9)
#   emergence            annual emergence from the seedbank (e, Eqs. 8 and 10)
#   bank_survival        seedbank persistence net of germination (sigma, Eq. 10)
#   intra_comp           intraspecific competition coefficient (a, Eq. 9)
#   overwinter_survival  overwinter survival of new seeds (phi, Eq. 9)
p_dicot = {
    'name':                'dicot_base',
    'max_seeds':           50000,
    'emergence':           0.05,
    'bank_survival':       0.40,
    'intra_comp':          0.41,
    'overwinter_survival': 0.05,
}

p_mono = {
    'name':                'mono_base',
    'max_seeds':           2500,
    'emergence':           0.20,
    'bank_survival':       0.50,
    'intra_comp':          0.25,
    'overwinter_survival': 0.10,
}


def get_sensitivity_params(p_base, multiplier, tag):
    """+/-20% sensitivity of the five population parameters. The intraspecific
    competition coefficient is divided rather than multiplied, so that all five
    shift weed pressure in the same direction."""
    return {
        'name':                p_base['name'] + '_' + tag,
        'max_seeds':           p_base['max_seeds'] * multiplier,
        'emergence':           min(1.0, p_base['emergence'] * multiplier),
        'bank_survival':       min(1.0, p_base['bank_survival'] * multiplier),
        'intra_comp':          p_base['intra_comp'] / multiplier,
        'overwinter_survival': min(1.0, p_base['overwinter_survival'] * multiplier),
    }


p_dicot_low  = get_sensitivity_params(p_dicot, 0.8, 'low')
p_dicot_high = get_sensitivity_params(p_dicot, 1.2, 'high')
p_mono_low   = get_sensitivity_params(p_mono,  0.8, 'low')
p_mono_high  = get_sensitivity_params(p_mono,  1.2, 'high')


# =====================================================================
# Decision scale-up used by the ET / DynET threshold.
#
# Returns the factor by which a raw detected count is multiplied to obtain the
# estimated true count, using the same patch-allocation logic as Year 0:
#     T / P = 1 + ((1/R) - 1) * PATCH_FRACTION
# At R = 0.41 and PATCH_FRACTION = 0.70 this equals 2.007.
# Returns 1.0 for the ET under ET_INFORMATION = 'detected' and for both
# threshold strategies under ET_FUTURE_RULE = 'v13'.
# =====================================================================
def _advance_one_year(bank_d, bank_m, emerged_d, emerged_m, treated,
                      kill_y, seed_penalty_if_treated, chem_cost, p_d, p_m):
    """One year of the population model for a single trajectory.

    Returns (new_bank_d, new_bank_m, cost) where cost is the undiscounted sum of
    herbicide and yield-loss cost for that year. Shared by the treated and
    untreated continuations of a trajectory.
    """
    if treated:
        survivors_d = emerged_d * (1.0 - kill_y)
        survivors_m = emerged_m * (1.0 - kill_y)
        seed_pen    = seed_penalty_if_treated
    else:
        survivors_d = emerged_d
        survivors_m = emerged_m
        seed_pen    = 1.0

    TCL_surv = weed_height_modifier * (
        survivors_d * dicot_competition_index +
        survivors_m * monocot_competition_index)
    if treated:
        TCL_surv *= SURVIVOR_COMPETITIVE_FACTOR
    YL_pct  = multispecies_yield_loss_pct(
        TCL_surv, yield_loss_initial_slope, yield_loss_asymptote_pct)
    yl_cost = (YL_pct / 100.0) * (weed_free_yield_kg_per_m2 * 10000) * crop_price_eur_per_kg

    # Seed production (Eq. 9) and seedbank update (Eq. 10): surviving bank
    # plus the new seeds
    seeds_d = (((p_d['max_seeds'] * survivors_d) /
                (1.0 + p_d['intra_comp'] * survivors_d)) * p_d['overwinter_survival']
               if survivors_d > 0 else 0.0) * seed_pen
    seeds_m = (((p_m['max_seeds'] * survivors_m) /
                (1.0 + p_m['intra_comp'] * survivors_m)) * p_m['overwinter_survival']
               if survivors_m > 0 else 0.0) * seed_pen

    return (bank_d * p_d['bank_survival'] + seeds_d,
            bank_m * p_m['bank_survival'] + seeds_m,
            chem_cost + yl_cost)


_BRANCH_STATS = {'max_states': 0, 'truncated': 0}


def _poisson_sf(k, mu):
    """P(X > k) for X ~ Poisson(mu). Used for the detection probability.

    Written out rather than imported so the script has no dependency beyond the
    Python standard library and whatever ships inside QGIS.
    """
    if mu <= 0.0:
        return 0.0
    if k < 0:
        return 1.0
    if mu > 700.0:            # exp(-mu) underflows; detection is certain
        return 1.0
    term = math.exp(-mu)
    cdf = term
    for i in range(1, k + 1):
        term *= mu / i
        cdf += term
        if cdf >= 1.0:
            return 0.0
    return max(0.0, 1.0 - cdf)


def _et_threshold_count(lam_d, lam_m, area, su, scenario='SCENARIO_1_ECON'):
    """Smallest number of DETECTED weeds in a cell that clears the economic
    threshold, given the current class mix.

    Inverts the rectangular hyperbola analytically rather than searching:

        YL* = 100 * herbicide_cost / (yield_goal * price)          target loss
        TCL* = YL* / (slope * (1 - YL*/A))                         from Eq. 3
        x*  = TCL* * area / (CIM * su * mean_competitive_index)     from Eq. 1-2

    Returns None when the threshold cannot be reached at any density, i.e. when
    the herbicide cost exceeds the value of the asymptotic yield loss.
    """
    lam_tot = lam_d + lam_m
    if lam_tot <= 0:
        return None

    value_per_ha = weed_free_yield_kg_per_m2 * 10000.0 * crop_price_eur_per_kg

    # The forward-looking ET (years 1-10) also counts the benefit of preventing
    # this cohort from setting seed, approximated as the same loss in the
    # following season discounted by one year. It therefore clears its
    # threshold at a lower weed density than the ET.
    effective_cost = MAX_HERB_COST_YR
    if scenario == 'SCENARIO_4_DYNET':
        effective_cost = MAX_HERB_COST_YR / (1.0 + 1.0 / (1.0 + DISCOUNT_RATE))

    yl_target = 100.0 * effective_cost / value_per_ha
    if yl_target >= yield_loss_asymptote_pct:
        return None                      # unreachable: never economic to spray

    denom = yield_loss_initial_slope * (1.0 - yl_target / yield_loss_asymptote_pct)
    if denom <= 0:
        return None
    tcl_target = yl_target / denom

    # competitive index of the average weed in this cell, given the class mix
    frac_d = lam_d / lam_tot
    ci_mean = (frac_d * dicot_competition_index +
               (1.0 - frac_d) * monocot_competition_index)
    if ci_mean <= 0:
        return None

    x_star = tcl_target * area / (weed_height_modifier * su * ci_mean)
    return max(1, int(math.ceil(x_star)))


def _p_treated(lam_d, lam_m, area, recall, scenario, su):
    """Probability that a cell is treated in a given future year.

    THE DETECTION MODEL
    -------------------
    The number of weeds standing in a cell is not known exactly; the simulation
    carries an expected density. Weeds are taken to be distributed at random
    within the cell, so the count N follows a Poisson distribution with mean
    lambda = density x cell area. The detector finds each weed independently
    with probability R (the empirical recall, 0.41). By the thinning property of
    the Poisson distribution the number DETECTED is then exactly

        X ~ Poisson(lambda * R)

    which handles fractional expected counts naturally, unlike a binomial.

    A presence-based sprayer treats the cell if it sees at least one weed:

        P(treated) = P(X >= 1) = 1 - exp(-lambda * R)

    An economic-threshold sprayer treats it only if enough weeds are seen for
    the expected yield loss to exceed the herbicide cost:

        P(treated) = P(X >= m),  m from _et_threshold_count()

    CONSEQUENCES, which are the point of the change:
      * a cell holding few weeds is usually missed, and is missed again the next
        year, so its seedbank compounds instead of being reset;
      * the probability rises steeply with density, so an escaping patch is
        eventually found and treated - escapes grow but are self-limiting;
      * because lambda scales with CELL AREA, the same weed density is harder to
        detect in a 25 x 25 cm cell (0.0625 m2) than in a 1 x 1 m cell. Finer
        spraying resolutions therefore leak more weeds, independently of the
        Year-0 effect already in the model.
    """
    return _detection_stats(lam_d, lam_m, area, recall, scenario, su)[0]


def _detection_stats(lam_d, lam_m, area, recall, scenario, su):
    """Treatment probability AND the conditional weed counts.

    Returns (p_treat, f_lo, f_hi) where

        p_treat = P(X >= m)                     X ~ Poisson(lambda * r)
        f_lo    = E[N | X <  m] / lambda        weeds present if NOT treated
        f_hi    = E[N | X >= m] / lambda        weeds present if treated

    N is the (Poisson) number of weeds in the cell and X the number detected.
    By the thinning property N - X ~ Poisson(lambda (1 - r)) is independent of
    X, so E[N | branch] = lambda (1 - r) + E[X | branch]. The two factors are
    applied to each class separately (class thinning is also independent, and
    E[X_d | X] = X * lambda_d / lambda), and they satisfy
    p_treat * f_hi + (1 - p_treat) * f_lo = 1 exactly.

    Under SURVIVOR_CONDITIONING = 'v16' (validation) both factors are 1.
    """
    lam_tot = lam_d + lam_m
    if lam_tot <= 0.0:
        return 0.0, 1.0, 1.0
    mu = lam_tot * recall              # mean number DETECTED in this cell

    if scenario in ('SCENARIO_2_SPOT', 'SCENARIO_3_SECURE'):
        m = 1
    else:
        m = _et_threshold_count(lam_d, lam_m, area, su, scenario)
        if m is None:
            return 0.0, 1.0, 1.0

    if mu <= 0.0:
        return 0.0, 1.0, 1.0

    if m == 1:
        p = (1.0 - math.exp(-mu)) if mu < 700.0 else 1.0
        e_lo = 0.0                                  # E[X ; X < 1] = 0
    else:
        p = _poisson_sf(m - 1, mu)
        # partial mean  E[X ; X < m] = sum_{k<m} k * pmf(k)
        e_lo = 0.0
        if mu < 700.0:
            term = math.exp(-mu)
            for k in range(1, m):
                term *= mu / k
                e_lo += k * term

    if SURVIVOR_CONDITIONING != 'thinned':
        return p, 1.0, 1.0

    # E[N | X < m]  = lam (1-r) + E[X | X < m]  with E[X | X<m] = e_lo / (1-p)
    # E[N | X >= m] = lam (1-r) + (mu - e_lo) / p
    if p >= 1.0 - 1e-15:
        return 1.0, 1.0, 1.0
    if p <= 1e-300:
        return 0.0, 1.0, 1.0
    f_lo = (1.0 - recall) + (e_lo / (1.0 - p)) / lam_tot
    f_hi = (1.0 - recall) + ((mu - e_lo) / p) / lam_tot
    return p, f_lo, f_hi


def _merge_states(states):
    """Collapse trajectories whose seedbanks have converged.

    Without this the state count doubles every year. Keys are the seedbank
    values rounded to MERGE_DIGITS significant figures; merged states carry the
    summed probability and the probability-weighted mean seedbank, so no
    probability mass is lost.
    """
    def sig(x):
        if x <= 0:
            return 0.0
        d = MERGE_DIGITS - int(math.floor(math.log10(x))) - 1
        return round(x, d)

    merged = {}
    for bd, bm, pr in states:
        key = (sig(bd), sig(bm))
        if key in merged:
            obd, obm, opr = merged[key]
            tot = opr + pr
            merged[key] = ((obd * opr + bd * pr) / tot,
                           (obm * opr + bm * pr) / tot, tot)
        else:
            merged[key] = (bd, bm, pr)

    out = list(merged.values())
    if len(out) > MAX_STATES:
        out.sort(key=lambda t: -t[2])
        keep, drop = out[:MAX_STATES], out[MAX_STATES:]
        lost = sum(t[2] for t in drop)
        _BRANCH_STATS['truncated'] += 1
        if lost > 0:                     # redistribute so probabilities sum to 1
            scale = 1.0 + lost / sum(t[2] for t in keep)
            keep = [(bd, bm, pr * scale) for bd, bm, pr in keep]
        out = keep
    return out


def et_decision_scaleup(recall, scenario='SCENARIO_1_ECON'):
    """Factor converting a detected count into an estimated true count.

    Returns 1.0 - i.e. no correction, the threshold acts on what the sensor
    reports - for the economic-threshold strategy when ET_INFORMATION is
    'detected'. The forward-looking strategy (DynET) always applies the
    correction, since it is defined as the informed rule.
    """
    if recall <= 0 or ET_FUTURE_RULE != 'consistent':
        return 1.0
    if scenario == 'SCENARIO_1_ECON' and ET_INFORMATION == 'detected':
        return 1.0
    return 1.0 + ((1.0 / recall) - 1.0) * PATCH_FRACTION


# =====================================================================
# CORE FUNCTION: two-class 10-year forward simulation
# =====================================================================
def calculate_10yr_liability_v2(
    surv_dens_d,
    surv_dens_m,
    p_d, p_m,
    yr0_penalty,
    pre_d=None,
    pre_m=None,
    scenario='SCENARIO_1_ECON',
    is_broadcast=False,
    recall=0.41,
    stunt_future=None,
    future_kill=0.95,
    dd_base_kill=None,
    cell_area_m2=1.0,         # needed to convert density -> count
    year_costs_out=None       # optional list[SIM_YEARS]; receives the
                              #         undiscounted expected cost of each year
):
    """Combined two-class 10-year discounted future weed liability (EUR ha-1).

    Returns an NPV. Divide by ANNUITY_FACTOR for the equivalent annual value
    used in the annualised partial budget (Eq. 6).
    """

    if stunt_future is None:
        stunt_future = STUNTED_SEED_PENALTY
    if dd_base_kill is None:
        dd_base_kill = DD_BASE_KILL

    # -- Year 0 seedbank initialisation (Eq. 8, uses PRE-herbicide density) --
    ref_d = pre_d if (pre_d is not None and pre_d > 0) else surv_dens_d
    if ref_d > 0:
        native_bank_d = ref_d / p_d['emergence']
        yr0_seeds_d   = ((p_d['max_seeds'] * surv_dens_d) /
                         (1.0 + p_d['intra_comp'] * surv_dens_d)
                         if surv_dens_d > 0 else 0.0)
        bank_d = native_bank_d * p_d['bank_survival'] + yr0_seeds_d * yr0_penalty * p_d['overwinter_survival']
    else:
        bank_d = 0.0

    ref_m = pre_m if (pre_m is not None and pre_m > 0) else surv_dens_m
    if ref_m > 0:
        native_bank_m = ref_m / p_m['emergence']
        yr0_seeds_m   = ((p_m['max_seeds'] * surv_dens_m) /
                         (1.0 + p_m['intra_comp'] * surv_dens_m)
                         if surv_dens_m > 0 else 0.0)
        bank_m = native_bank_m * p_m['bank_survival'] + yr0_seeds_m * yr0_penalty * p_m['overwinter_survival']
    else:
        bank_m = 0.0

    # ---------------------------------------------------------------------
    # YEARS 1-10
    #
    # Under FUTURE_DETECTION = 'expected' (validation) this is a single
    # trajectory. Under 'probabilistic' the cell branches each year into a
    # treated and an untreated continuation, weighted by the probability that
    # the sprayer actually sees the weeds present (see _p_treated). The
    # reported liability is the probability-weighted expectation over all
    # surviving trajectories.
    # ---------------------------------------------------------------------
    su = et_decision_scaleup(recall, scenario)
    probabilistic = (FUTURE_DETECTION == 'probabilistic'
                     and not is_broadcast)

    # each state: (seedbank_dicot, seedbank_monocot, probability)
    states = [(bank_d, bank_m, 1.0)]
    total_npv_cost = 0.0

    for year in range(1, SIM_YEARS + 1):
        discount = 1.0 / ((1.0 + DISCOUNT_RATE) ** year)
        year_cost = 0.0
        next_states = []

        for bd, bm, prob in states:
            emerged_d = bd * p_d['emergence']
            emerged_m = bm * p_m['emergence']

            # -- probability that the cell is treated this year --------------
            f_lo = f_hi = 1.0            # no conditioning unless set below
            if is_broadcast or scenario == 'SCENARIO_3_SECURE':
                # both always pass over the whole field; for DD the BRANCH is
                # full dose vs background dose, not treated vs untreated
                p_treat = 1.0
            elif not probabilistic:
                # 'expected': detected = emerged * recall, above zero whenever
                # any weed emerged, so the decision is deterministic
                det_d = emerged_d * recall * su
                det_m = emerged_m * recall * su
                TCL_det = weed_height_modifier * (
                    det_d * dicot_competition_index +
                    det_m * monocot_competition_index)
                YL_det = multispecies_yield_loss_pct(
                    TCL_det, yield_loss_initial_slope, yield_loss_asymptote_pct)
                loss_ha = (YL_det / 100.0) * weed_free_yield_kg_per_m2 * 10000 * crop_price_eur_per_kg
                if scenario == 'SCENARIO_2_SPOT':
                    fires = (emerged_d + emerged_m) * recall > 0
                elif scenario == 'SCENARIO_4_DYNET':
                    fires = (loss_ha + loss_ha / (1.0 + DISCOUNT_RATE)) > MAX_HERB_COST_YR
                else:
                    fires = loss_ha > MAX_HERB_COST_YR
                p_treat = 1.0 if fires else 0.0
            else:
                p_treat, f_lo, f_hi = _detection_stats(
                    emerged_d * cell_area_m2, emerged_m * cell_area_m2,
                    cell_area_m2, recall, scenario, su)

            # f_lo / f_hi are the weeds present in each branch relative
            # to the expected count. Under 'v16', 'expected' and broadcast they
            # are 1.

            # -- evaluate both continuations ---------------------------------
            for treated, branch_p in ((True, p_treat), (False, 1.0 - p_treat)):
                if branch_p <= 0.0:
                    continue
                w = prob * branch_p
                if w < PRUNE_PROB:
                    continue

                if scenario == 'SCENARIO_3_SECURE' and not is_broadcast:
                    # DD: 'treated' means the cell was SEEN and got the full
                    # dose; otherwise it still receives the background dose.
                    if probabilistic:
                        p_full, g_lo, g_hi = _detection_stats(
                            emerged_d * cell_area_m2, emerged_m * cell_area_m2,
                            cell_area_m2, recall, 'SCENARIO_2_SPOT', su)
                    else:
                        p_full = 1.0 if (emerged_d + emerged_m) * recall > 0 else 0.0
                        g_lo = g_hi = 1.0
                    # branch over dose instead of over treatment
                    for full_dose, dose_p in ((True, p_full), (False, 1.0 - p_full)):
                        wd = prob * dose_p
                        if dose_p <= 0.0 or wd < PRUNE_PROB:
                            continue
                        if full_dose:
                            kill_y = future_kill
                            chem   = MAX_HERB_COST_YR
                            g      = g_hi
                        else:
                            kill_y = dd_base_kill
                            chem   = (MAX_HERB_COST_YR if DD_FUTURE_COST_MODE == 'v13'
                                      else MAX_HERB_COST_YR * DD_BASE_DOSE)
                            g      = g_lo
                        nbd, nbm, cost = _advance_one_year(
                            bd, bm, emerged_d * g, emerged_m * g, True, kill_y,
                            stunt_future, chem, p_d, p_m)
                        year_cost += wd * cost
                        next_states.append((nbd, nbm, wd))
                    break      # dose branching replaces treatment branching

                kill_y = future_kill
                chem   = MAX_HERB_COST_YR if treated else 0.0
                f      = f_hi if treated else f_lo
                nbd, nbm, cost = _advance_one_year(
                    bd, bm, emerged_d * f, emerged_m * f, treated, kill_y,
                    stunt_future, chem, p_d, p_m)
                year_cost += w * cost
                next_states.append((nbd, nbm, w))

        total_npv_cost += year_cost * discount
        if year_costs_out is not None:
            year_costs_out[year - 1] += year_cost
        states = _merge_states(next_states) if len(next_states) > 1 else next_states
        _BRANCH_STATS['max_states'] = max(_BRANCH_STATS['max_states'], len(states))

    return total_npv_cost


def _init_bank(surv_dens, pre_dens, yr0_penalty, p):
    """Year-1 seedbank of one class, exactly as in calculate_10yr_liability_v2."""
    ref = pre_dens if (pre_dens is not None and pre_dens > 0) else surv_dens
    if ref <= 0:
        return 0.0
    seeds = ((p['max_seeds'] * surv_dens) / (1.0 + p['intra_comp'] * surv_dens)
             if surv_dens > 0 else 0.0)
    return ref / p['emergence'] * p['bank_survival'] + seeds * yr0_penalty * p['overwinter_survival']


def _merge_paired(states):
    """_merge_states for paired trajectories (bd, bm, bd1, bm1, prob): keyed on
    the BASE seedbanks, so the base part merges exactly as in the unpaired
    model; the +1-weed seedbanks are carried as probability-weighted means."""
    def sig(x):
        if x <= 0:
            return 0.0
        d = MERGE_DIGITS - int(math.floor(math.log10(x))) - 1
        return round(x, d)
    merged = {}
    for bd, bm, bd1, bm1, pr in states:
        key = (sig(bd), sig(bm))
        if key in merged:
            o = merged[key]; tot = o[4] + pr
            merged[key] = ((o[0] * o[4] + bd * pr) / tot, (o[1] * o[4] + bm * pr) / tot,
                           (o[2] * o[4] + bd1 * pr) / tot, (o[3] * o[4] + bm1 * pr) / tot, tot)
        else:
            merged[key] = (bd, bm, bd1, bm1, pr)
    out = list(merged.values())
    if len(out) > MAX_STATES:
        out.sort(key=lambda t: -t[4])
        keep, drop = out[:MAX_STATES], out[MAX_STATES:]
        lost = sum(t[4] for t in drop)
        _BRANCH_STATS['truncated'] += 1
        if lost > 0:
            scale = 1.0 / (1.0 - lost)
            keep = [(a, b, c, d, pr * scale) for a, b, c, d, pr in keep]
        out = keep
    return out


def liability_difference_fixed_mgmt(
    base_surv_d, base_surv_m, base_pre_d, base_pre_m,
    plus_surv_d, plus_surv_m, plus_pre_d, plus_pre_m,
    p_d, p_m, yr0_penalty, scenario='SCENARIO_1_ECON', recall=0.41,
    cell_area_m2=1.0, return_base=False
):
    """NPV (EUR/ha) of liability(+1 weed) - liability(base) with the
    future management of the BASE cell applied to both.

    Both cells are simulated along the same tree of treated / untreated years.
    Branch probabilities, the detection-conditioned fractions (f, g) and the
    ET / DD decisions are computed from the base cell only. With return_base the
    base liability is returned as well; it equals calculate_10yr_liability_v2
    for the same inputs (self-check)."""
    stunt_future = STUNTED_SEED_PENALTY
    future_kill = full_kill_rate
    dd_base_kill = DD_BASE_KILL
    states = [(_init_bank(base_surv_d, base_pre_d, yr0_penalty, p_d),
               _init_bank(base_surv_m, base_pre_m, yr0_penalty, p_m),
               _init_bank(plus_surv_d, plus_pre_d, yr0_penalty, p_d),
               _init_bank(plus_surv_m, plus_pre_m, yr0_penalty, p_m), 1.0)]
    su = et_decision_scaleup(recall, scenario)
    probabilistic = FUTURE_DETECTION == 'probabilistic'
    npv_base = 0.0; npv_diff = 0.0
    for year in range(1, SIM_YEARS + 1):
        discount = 1.0 / ((1.0 + DISCOUNT_RATE) ** year)
        nxt = []
        for bd, bm, bd1, bm1, prob in states:
            ed, em = bd * p_d['emergence'], bm * p_m['emergence']
            ed1, em1 = bd1 * p_d['emergence'], bm1 * p_m['emergence']
            f_lo = f_hi = 1.0
            if scenario == 'SCENARIO_3_SECURE':
                p_treat = 1.0
            elif not probabilistic:
                det_d = ed * recall * su; det_m = em * recall * su
                TCL_det = weed_height_modifier * (det_d * dicot_competition_index +
                                                  det_m * monocot_competition_index)
                YL_det = multispecies_yield_loss_pct(TCL_det, yield_loss_initial_slope,
                                                     yield_loss_asymptote_pct)
                loss_ha = (YL_det / 100.0) * weed_free_yield_kg_per_m2 * 10000 * crop_price_eur_per_kg
                if scenario == 'SCENARIO_2_SPOT':
                    fires = (ed + em) * recall > 0
                elif scenario == 'SCENARIO_4_DYNET':
                    fires = (loss_ha + loss_ha / (1.0 + DISCOUNT_RATE)) > MAX_HERB_COST_YR
                else:
                    fires = loss_ha > MAX_HERB_COST_YR
                p_treat = 1.0 if fires else 0.0
            else:
                p_treat, f_lo, f_hi = _detection_stats(
                    ed * cell_area_m2, em * cell_area_m2, cell_area_m2, recall, scenario, su)
            for treated, branch_p in ((True, p_treat), (False, 1.0 - p_treat)):
                if branch_p <= 0.0:
                    continue
                w = prob * branch_p
                if w < PRUNE_PROB:
                    continue
                if scenario == 'SCENARIO_3_SECURE':
                    if probabilistic:
                        p_full, g_lo, g_hi = _detection_stats(
                            ed * cell_area_m2, em * cell_area_m2, cell_area_m2,
                            recall, 'SCENARIO_2_SPOT', su)
                    else:
                        p_full = 1.0 if (ed + em) * recall > 0 else 0.0
                        g_lo = g_hi = 1.0
                    for full_dose, dose_p in ((True, p_full), (False, 1.0 - p_full)):
                        wd = prob * dose_p
                        if dose_p <= 0.0 or wd < PRUNE_PROB:
                            continue
                        if full_dose:
                            kill_y, chem, g = future_kill, MAX_HERB_COST_YR, g_hi
                        else:
                            kill_y = dd_base_kill
                            chem = (MAX_HERB_COST_YR if DD_FUTURE_COST_MODE == 'v13'
                                    else MAX_HERB_COST_YR * DD_BASE_DOSE)
                            g = g_lo
                        a = _advance_one_year(bd, bm, ed * g, em * g, True, kill_y,
                                              stunt_future, chem, p_d, p_m)
                        b = _advance_one_year(bd1, bm1, ed1 * g, em1 * g, True, kill_y,
                                              stunt_future, chem, p_d, p_m)
                        npv_base += wd * a[2] * discount
                        npv_diff += wd * (b[2] - a[2]) * discount
                        nxt.append((a[0], a[1], b[0], b[1], wd))
                    break
                chem = MAX_HERB_COST_YR if treated else 0.0
                f = f_hi if treated else f_lo
                a = _advance_one_year(bd, bm, ed * f, em * f, treated, future_kill,
                                      stunt_future, chem, p_d, p_m)
                b = _advance_one_year(bd1, bm1, ed1 * f, em1 * f, treated, future_kill,
                                      stunt_future, chem, p_d, p_m)
                npv_base += w * a[2] * discount
                npv_diff += w * (b[2] - a[2]) * discount
                nxt.append((a[0], a[1], b[0], b[1], w))
        states = _merge_paired(nxt) if len(nxt) > 1 else nxt
    return (npv_diff, npv_base) if return_base else npv_diff


def multispecies_yield_loss_pct(TCL, slope, A):
    if TCL <= 0:
        return 0.0
    linear_impact = slope * TCL
    return min(A, linear_impact / (1.0 + (linear_impact / A)))


def year0_revenue(S_d, S_m, area, sprayed):
    """Year-0 revenue of a cell (Eqs. 3-5) from its survivors S_d, S_m."""
    TCL_surv = weed_height_modifier * ((S_d * dicot_competition_index +
                                        S_m * monocot_competition_index) / area)
    if sprayed:
        TCL_surv *= SURVIVOR_COMPETITIVE_FACTOR
    YL = multispecies_yield_loss_pct(TCL_surv, yield_loss_initial_slope,
                                     yield_loss_asymptote_pct)
    return (weed_free_yield_kg_per_m2 * (1.0 - YL / 100.0)) * area * crop_price_eur_per_kg


def ensure_field_exists(layer, name, qtype):
    if layer.fields().indexFromName(name) == -1:
        layer.dataProvider().addAttributes([QgsField(name, qtype)])
        layer.updateFields()


# =====================================================================
# YEAR-0 DECISION RULE AND RECALL-CURVE BRANCHES
# =====================================================================
def _year0_decision(det_d, det_m, area, scenario):
    """Year-0 spray decision of the ET, PB and DD strategies from the DETECTED
    counts (det_d, det_m) in a cell. Returns (spray, dose). Mirrors the rule in
    _compute_cell_uncached exactly, so the baseline and the recall curve use the
    same code path. DynET is not covered (it needs the liability engine and is
    excluded from the curve).
    """
    det_tot = det_d + det_m
    if scenario == 'SCENARIO_2_SPOT':
        return (det_tot > 0), (1.0 if det_tot > 0 else 0.0)
    if scenario == 'SCENARIO_3_SECURE':
        return True, (1.0 if det_tot > 0 else DD_BASE_DOSE)
    if scenario == 'SCENARIO_1_ECON':
        if det_tot <= 0:
            return False, 0.0
        if ET_INFORMATION == 'detected':
            dd, dm = det_d / area, det_m / area
        else:
            sc = 1.0 + ((1.0 / recall_dicot) - 1.0) * PATCH_FRACTION
            dd, dm = det_d * sc / area, det_m * sc / area
        tcl = weed_height_modifier * (dd * dicot_competition_index +
                                      dm * monocot_competition_index)
        yl  = multispecies_yield_loss_pct(tcl, yield_loss_initial_slope,
                                          yield_loss_asymptote_pct)
        loss = (yl / 100.0) * weed_free_yield_kg_per_m2 * crop_price_eur_per_kg
        return (loss > herbicide_cost_fullrate_eur_per_m2), \
               (1.0 if loss > herbicide_cost_fullrate_eur_per_m2 else 0.0)
    raise ValueError('no year-0 rule for ' + scenario)


def _poisson_support(mu, cutoff):
    """[(k, p)] for k in a window of a Poisson(mu) holding 1 - cutoff of the
    mass. Log-space so large mu does not underflow."""
    if mu <= 0.0:
        return [(0, 1.0)]
    sd = math.sqrt(mu)
    lo = max(0, int(math.floor(mu - 8.0 * sd)))
    hi = int(math.ceil(mu + 8.0 * sd)) + 2
    out = []
    lmu = math.log(mu)
    for k in range(lo, hi + 1):
        lp = -mu + k * lmu - math.lgamma(k + 1)
        if lp > -745.0:
            p = math.exp(lp)
            if p >= cutoff:
                out.append((k, p))
    tot = sum(p for _, p in out)
    return [(k, p / tot) for k, p in out] if tot > 0 else [(0, 1.0)]


def _binomial_support(n, s, cutoff):
    """[(k, p)] for k ~ Binomial(n, s), n integer, in an 8-sd window."""
    if n <= 0:
        return [(0, 1.0)]
    if s >= 1.0:
        return [(n, 1.0)]
    if s <= 0.0:
        return [(0, 1.0)]
    mean = n * s
    sd = math.sqrt(n * s * (1.0 - s))
    lo = max(0, int(math.floor(mean - 8.0 * sd)))
    hi = min(n, int(math.ceil(mean + 8.0 * sd)) + 1)
    out = []
    ls, l1s = math.log(s), math.log(1.0 - s)
    lgn = math.lgamma(n + 1)
    for k in range(lo, hi + 1):
        lp = lgn - math.lgamma(k + 1) - math.lgamma(n - k + 1) + k * ls + (n - k) * l1s
        if lp > -745.0:
            p = math.exp(lp)
            if p >= cutoff:
                out.append((k, p))
    tot = sum(p for _, p in out)
    return [(k, p / tot) for k, p in out] if tot > 0 else [(n, 1.0)]


def _year0_branches(P_d, P_m, T_d, T_m, area, r_new, scenario):
    """Year-0 outcome of a cell if a sensor of recall r_new had mapped
    the SAME true population (T_d, T_m) that the measured sensor (recall r0,
    detections P_d, P_m) implies.

    Returns a list of branches (prob, spray, dose, Td, Tm) where Td/Tm are the
    expected weeds present in the cell given that branch. Branches with equal
    (spray, dose) are merged. At r_new == r0 the single branch reproduces the
    baseline exactly.

    r_new > r0  (nested sensor): the new sensor finds each of the M = T - P
        missed weeds independently with probability q = (r' - r0)/(1 - r0).
        The missed count is taken as Poisson(M) so that fractional expected
        counts (the isolated-escape 'drip' in empty cells) are handled, and the
        weeds present in a branch are P + M(1-q) + E[Y | branch].
    r_new < r0  (thinned sensor): each observed detection is kept with
        probability s = r'/r0. The true population is unchanged, so both
        branches hold T weeds; only the spray decision changes.
    """
    r0 = recall_dicot
    base_spray, base_dose = _year0_decision(P_d, P_m, area, scenario)

    if abs(r_new - r0) < 1e-12:
        return [(1.0, base_spray, base_dose, T_d, T_m)]

    M_d = max(0.0, T_d - P_d)
    M_m = max(0.0, T_m - P_m)

    if r_new > r0:
        q = (r_new - r0) / (1.0 - r0)
        # More detections can only raise a monotone rule; a cell already at
        # full dose cannot change.
        if base_spray and base_dose == 1.0:
            return [(1.0, True, 1.0, T_d, T_m)]
        sup_d = _poisson_support(M_d * q, CURVE_PRUNE_PROB)
        sup_m = _poisson_support(M_m * q, CURVE_PRUNE_PROB)
        acc = {}
        for yd, pd_ in sup_d:
            for ym, pm_ in sup_m:
                pr = pd_ * pm_
                if pr < CURVE_PRUNE_PROB:
                    continue
                spray, dose = _year0_decision(P_d + yd, P_m + ym, area, scenario)
                key = (spray, dose)
                a = acc.get(key)
                if a is None:
                    acc[key] = [pr, pr * yd, pr * ym]
                else:
                    a[0] += pr; a[1] += pr * yd; a[2] += pr * ym
        tot = sum(a[0] for a in acc.values())
        out = []
        for (spray, dose), (pr, syd, sym) in acc.items():
            Td = P_d + M_d * (1.0 - q) + syd / pr
            Tm = P_m + M_m * (1.0 - q) + sym / pr
            out.append((pr / tot, spray, dose, Td, Tm))
        return out

    # r_new < r0: thin the observed detections
    s = r_new / r0
    if not base_spray:
        # ET/PB cell not sprayed at r0: fewer detections cannot change that
        return [(1.0, False, 0.0, T_d, T_m)]
    if scenario == 'SCENARIO_3_SECURE' and base_dose < 1.0:
        # DD background-dose cell (no detection): stays background dose
        return [(1.0, True, DD_BASE_DOSE, T_d, T_m)]
    n_d, n_m = int(round(P_d)), int(round(P_m))
    if scenario in ('SCENARIO_2_SPOT', 'SCENARIO_3_SECURE'):
        p_det = 1.0 - (1.0 - s) ** (n_d + n_m)
        if scenario == 'SCENARIO_2_SPOT':
            return [(p_det, True, 1.0, T_d, T_m), (1.0 - p_det, False, 0.0, T_d, T_m)]
        return [(p_det, True, 1.0, T_d, T_m), (1.0 - p_det, True, DD_BASE_DOSE, T_d, T_m)]
    sup_d = _binomial_support(n_d, s, CURVE_PRUNE_PROB)
    sup_m = _binomial_support(n_m, s, CURVE_PRUNE_PROB)
    acc = {}
    for xd, pd_ in sup_d:
        for xm, pm_ in sup_m:
            pr = pd_ * pm_
            if pr < CURVE_PRUNE_PROB:
                continue
            key = _year0_decision(xd, xm, area, scenario)
            acc[key] = acc.get(key, 0.0) + pr
    tot = sum(acc.values())
    return [(pr / tot, spray, dose, T_d, T_m) for (spray, dose), pr in acc.items()]


# =====================================================================
# CELL-LEVEL MEMOISATION
#
# Everything below is a deterministic function of (P_d, P_m, area) plus the
# module-level configuration. Distinct states are few, so each is computed once.
#
# Field-level constants (drip rates etc.) are set once before the loop and are
# constant within a run, so they do not need to be part of the key.
# =====================================================================
_CELL_CACHE = {}
_CACHE_STATS = {'hit': 0, 'miss': 0}


def compute_cell(P_d, P_m, area):
    if USE_CELL_CACHE:
        key = (P_d, P_m, round(area, 6))
        cached = _CELL_CACHE.get(key)
        if cached is not None:
            _CACHE_STATS['hit'] += 1
            return cached
    _CACHE_STATS['miss'] += 1
    res = _compute_cell_uncached(P_d, P_m, area)
    if USE_CELL_CACHE:
        _CELL_CACHE[key] = res
    return res


# ALLOCATION OF UNDETECTED WEEDS (Eq. 7)
def _iso_param(total_detected, r, patch_frac, n_empty):
    """Allocation parameter for the isolated share of one weed class.

    'area'        : background density of isolated escapes, weeds per m2 of field.
    'empty_cells' : weeds per zero-detection cell (validation).
    """
    missed = (total_detected / r - total_detected) if r > 0 else 0.0
    iso = missed * (1.0 - patch_frac)
    if ALLOCATION_RULE == 'area':
        return iso / total_area_m2 if total_area_m2 > 0 else 0.0
    return iso / n_empty if n_empty > 0 else 0.0


def _true_count(P, r, patch_frac, iso, area):
    """Estimated true count of one weed class in one cell (Eq. 7)."""
    patch = (P + P * ((1.0 / r) - 1.0) * patch_frac) if P > 0 else 0.0
    if ALLOCATION_RULE == 'area':
        return patch + iso * area
    return patch if P > 0 else iso


def _compute_cell_uncached(P_d, P_m, area):
    """All per-cell economics. Returns a dict of contributions to accumulate."""
    area_ha = area / 10000.0
    out = {}

    # -- True count estimation, baseline 70/30 split (Eq. 7) ----------------
    T_d = _true_count(P_d, recall_dicot,   0.70, isolated_drip_d, area)
    T_m = _true_count(P_m, recall_monocot, 0.70, isolated_drip_m, area)

    bio_dens_d = T_d / area
    bio_dens_m = T_m / area

    # -- Year 0 spray decision ----------------------------------------------
    # The DECISION density depends on what information the rule is
    # assumed to act on. Under ET_INFORMATION = 'detected' the economic
    # threshold sees only the raw detected counts, as a commercial system would;
    # under 'true' it sees the recall-corrected estimate. The BIOLOGY always
    # uses the estimated true counts (bio_dens_*): the weeds that are present
    # cause yield loss and set seed whether or not the sensor found them.
    if SCENARIO == 'SCENARIO_1_ECON' and ET_INFORMATION == 'detected':
        dec_dens_d, dec_dens_m = P_d / area, P_m / area
    else:
        dec_dens_d, dec_dens_m = bio_dens_d, bio_dens_m

    map_TCL = weed_height_modifier * (dec_dens_d * dicot_competition_index +
                                      dec_dens_m * monocot_competition_index)
    YL_decision = multispecies_yield_loss_pct(map_TCL, yield_loss_initial_slope,
                                              yield_loss_asymptote_pct)
    loss_eur_per_m2 = (YL_decision / 100.0) * weed_free_yield_kg_per_m2 * crop_price_eur_per_kg

    spray = False; dose = 0.0; kill = 0.0

    if SCENARIO == 'SCENARIO_1_ECON':
        if (loss_eur_per_m2 > herbicide_cost_fullrate_eur_per_m2) and ((P_d + P_m) > 0):
            spray = True; dose = 1.0

    elif SCENARIO == 'SCENARIO_2_SPOT':
        if (P_d + P_m) > 0:
            spray = True; dose = 1.0

    elif SCENARIO == 'SCENARIO_3_SECURE':
        spray = True
        dose  = 1.0 if (P_d + P_m) > 0 else DD_BASE_DOSE

    elif SCENARIO == 'SCENARIO_4_DYNET':
        # FORWARD-LOOKING ECONOMIC THRESHOLD
        #
        # Compare full economics of both branches instead of Year-0 yield loss
        # alone. This internalises the shadow cost of letting the cell set seed
        # without needing it as a separate input:
        #
        #   spray if (rev_spray - herb - liab_spray) > (rev_nospray - liab_nospray)
        #
        # As with ET, a cell can only be sprayed if at least one weed was
        # actually DETECTED; a cell without detections cannot be targeted.
        if (P_d + P_m) > 0:
            S_d_s = T_d * (1.0 - full_kill_rate)
            S_m_s = T_m * (1.0 - full_kill_rate)
            rev_s = year0_revenue(S_d_s, S_m_s, area, True)
            liab_s = calculate_10yr_liability_v2(
                S_d_s / area, S_m_s / area, p_dicot, p_mono, STUNTED_SEED_PENALTY,
                pre_d=bio_dens_d, pre_m=bio_dens_m,
                scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha / ANNUITY_FACTOR

            rev_n = year0_revenue(T_d, T_m, area, False)
            liab_n = calculate_10yr_liability_v2(
                T_d / area, T_m / area, p_dicot, p_mono, 1.0,
                pre_d=bio_dens_d, pre_m=bio_dens_m,
                scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha / ANNUITY_FACTOR

            herb_if_spray = herbicide_cost_fullrate_eur_per_m2 * area
            if (rev_s - herb_if_spray - liab_s) > (rev_n - liab_n):
                spray = True; dose = 1.0

    if spray:
        kill      = full_kill_rate if dose == 1.0 else DD_BASE_KILL
        real_herb = (herbicide_cost_fullrate_eur_per_m2 * area) * dose
        yr0_penalty = STUNTED_SEED_PENALTY
    else:
        kill = 0.0; dose = 0.0; real_herb = 0.0
        yr0_penalty = 1.0

    S_d    = T_d * (1.0 - kill)
    S_m    = T_m * (1.0 - kill)
    dens_d = S_d / area
    dens_m = S_m / area

    # -- FUTURE LIABILITY, SSWM ---------------------------------------------
    _yc_s = [0.0] * SIM_YEARS
    shad_sswm_base = calculate_10yr_liability_v2(
        dens_d, dens_m, p_dicot, p_mono, yr0_penalty,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area,
        year_costs_out=_yc_s) * area_ha

    shad_sswm_low = calculate_10yr_liability_v2(
        dens_d, dens_m, p_dicot_low, p_mono_low, yr0_penalty,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha

    shad_sswm_high = calculate_10yr_liability_v2(
        dens_d, dens_m, p_dicot_high, p_mono_high, yr0_penalty,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha

    # -- MARGINAL SHADOW COST ------------------------------------------------
    # Ten-year NPV (EUR per weed) of one additional dicot (marg_sc_d) or monocot
    # (marg_sc_m) added to the year-0 population before the herbicide
    # application. The added weed receives the same year-0 treatment as the
    # cell (survives at 1 - kill, seed set yr0_penalty).
    if SHADOW_COST_MODE == 'fixed_management':
        # Same future treatment events as the base cell; the dormant seedbank
        # back-calculated from density (Eq. 8) is left at its base value.
        marg_sc_d = liability_difference_fixed_mgmt(
            dens_d, dens_m, bio_dens_d, bio_dens_m,
            (T_d + 1.0) * (1.0 - kill) / area, dens_m, bio_dens_d, bio_dens_m,
            p_dicot, p_mono, yr0_penalty, scenario=SCENARIO, recall=recall_dicot,
            cell_area_m2=area) * area_ha
        marg_sc_m = liability_difference_fixed_mgmt(
            dens_d, dens_m, bio_dens_d, bio_dens_m,
            dens_d, (T_m + 1.0) * (1.0 - kill) / area, bio_dens_d, bio_dens_m,
            p_dicot, p_mono, yr0_penalty, scenario=SCENARIO, recall=recall_dicot,
            cell_area_m2=area) * area_ha
    else:   # 'adaptive' (validation): two separate simulations
        npv_plus_d = calculate_10yr_liability_v2(
            (T_d + 1.0) * (1.0 - kill) / area, dens_m,
            p_dicot, p_mono, yr0_penalty,
            pre_d=bio_dens_d + (1.0 / area), pre_m=bio_dens_m,
            scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha
        marg_sc_d = npv_plus_d - shad_sswm_base
        npv_plus_m = calculate_10yr_liability_v2(
            dens_d, (T_m + 1.0) * (1.0 - kill) / area,
            p_dicot, p_mono, yr0_penalty,
            pre_d=bio_dens_d, pre_m=bio_dens_m + (1.0 / area),
            scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha
        marg_sc_m = npv_plus_m - shad_sswm_base

    if ANNUALISE_SHADOW_COST:
        marg_sc_d /= ANNUITY_FACTOR
        marg_sc_m /= ANNUITY_FACTOR

    # -- SPLIT SENSITIVITY, SSWM --------------------------------------------
    split_sswm = {}; split_rev_sswm = {}
    for split_name, (patch_frac, drip_d_var, drip_m_var) in drip_rates.items():
        T_d_v = _true_count(P_d, recall_dicot,   patch_frac, drip_d_var, area)
        T_m_v = _true_count(P_m, recall_monocot, patch_frac, drip_m_var, area)
        split_rev_sswm[split_name] = year0_revenue(
            T_d_v * (1.0 - kill), T_m_v * (1.0 - kill), area, spray)
        split_sswm[split_name] = calculate_10yr_liability_v2(
            T_d_v * (1.0 - kill) / area, T_m_v * (1.0 - kill) / area,
            p_dicot, p_mono, yr0_penalty,
            pre_d=T_d_v / area, pre_m=T_m_v / area,
            scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha

    # -- KILL SENSITIVITY, SSWM ---------------------------------------------
    kill_sswm = {}
    for k_name, k_params in KILL_SCENARIOS.items():
        # Year-0 kill must respect the applied DOSE. Under DD, a
        # background-dose cell kills at DD_BASE_KILL, not at the full-rate
        # efficacy being varied here; the background efficacy has its own
        # dedicated sensitivity (DD_KILL_SWEEP). Without this guard
        # the k_base column would no longer reproduce the DD baseline once
        # DD_BASE_KILL != full_kill_rate. ET/PB/DynET are unaffected: they
        # always spray at dose == 1.0. SELF-CHECK: for every strategy,
        # Kill_Base_Liab must equal the baseline liability of the same run.
        kill_yr0  = ((k_params['kill'] if dose == 1.0 else DD_BASE_KILL)
                     if spray else 0.0)
        stunt_yr0 = k_params['stunt'] if spray else 1.0
        kill_sswm[k_name] = calculate_10yr_liability_v2(
            T_d * (1.0 - kill_yr0) / area, T_m * (1.0 - kill_yr0) / area,
            p_dicot, p_mono, stunt_yr0,
            pre_d=bio_dens_d, pre_m=bio_dens_m,
            scenario=SCENARIO, recall=recall_dicot,
            stunt_future=k_params['stunt'],
            future_kill=k_params['kill'], cell_area_m2=area) * area_ha

    # -- RECALL SENSITIVITY (+/-20%, future liability only), SSWM ----------
    recall_sswm = {}
    for rname, (rv, rdrip_d, rdrip_m) in recall_precomp.items():
        T_d_r = _true_count(P_d, rv, 0.70, rdrip_d, area)
        T_m_r = _true_count(P_m, rv, 0.70, rdrip_m, area)
        recall_sswm[rname] = calculate_10yr_liability_v2(
            T_d_r * (1.0 - kill) / area, T_m_r * (1.0 - kill) / area,
            p_dicot, p_mono, yr0_penalty,
            pre_d=T_d_r / area, pre_m=T_m_r / area,
            scenario=SCENARIO, recall=rv, cell_area_m2=area) * area_ha

    # -- RECALL SWEEP: liability AND revenue ------------------------
    sweep_l_sswm = {}; sweep_r_sswm = {}
    sweep_l_bc   = {}; sweep_r_bc   = {}
    if ENABLE_RECALL_SWEEP:
        for rv in RECALL_SWEEP:
            rdrip_d, rdrip_m = sweep_precomp[rv]
            T_d_r = _true_count(P_d, rv, 0.70, rdrip_d, area)
            T_m_r = _true_count(P_m, rv, 0.70, rdrip_m, area)

            S_d_r = T_d_r * (1.0 - kill)
            S_m_r = T_m_r * (1.0 - kill)
            sweep_r_sswm[rv] = year0_revenue(S_d_r, S_m_r, area, spray)
            sweep_l_sswm[rv] = calculate_10yr_liability_v2(
                S_d_r / area, S_m_r / area, p_dicot, p_mono, yr0_penalty,
                pre_d=T_d_r / area, pre_m=T_m_r / area,
                scenario=SCENARIO, recall=rv, cell_area_m2=area) * area_ha

            S_d_b = T_d_r * (1.0 - full_kill_rate)
            S_m_b = T_m_r * (1.0 - full_kill_rate)
            sweep_r_bc[rv] = year0_revenue(S_d_b, S_m_b, area, True)
            sweep_l_bc[rv] = calculate_10yr_liability_v2(
                S_d_b / area, S_m_b / area, p_dicot, p_mono, STUNTED_SEED_PENALTY,
                pre_d=T_d_r / area, pre_m=T_m_r / area,
                is_broadcast=True, recall=rv, cell_area_m2=area) * area_ha

    # -- RECALL CURVE: sensor of recall r' on the SAME population ----
    # For each r' the year-0 map is re-derived (see _year0_branches) and the
    # years 1-10 are simulated at r'. Broadcast is independent of the sensor
    # and is NOT recomputed. Each entry: dict(herb, rev, liab, spr_area, vol).
    curve = {}
    if ENABLE_RECALL_CURVE and SCENARIO in ('SCENARIO_1_ECON', 'SCENARIO_2_SPOT',
                                           'SCENARIO_3_SECURE'):
        for rv in RECALL_CURVE:
            c_herb = c_rev = c_liab = c_spr = c_vol = 0.0
            for pr, b_spray, b_dose, Td_b, Tm_b in _year0_branches(
                    P_d, P_m, T_d, T_m, area, rv, SCENARIO):
                if pr <= 0.0:
                    continue
                if b_spray:
                    b_kill = full_kill_rate if b_dose == 1.0 else DD_BASE_KILL
                    b_herb = herbicide_cost_fullrate_eur_per_m2 * area * b_dose
                    b_pen  = STUNTED_SEED_PENALTY
                else:
                    b_kill = 0.0; b_herb = 0.0; b_pen = 1.0
                S_d_b = Td_b * (1.0 - b_kill)
                S_m_b = Tm_b * (1.0 - b_kill)
                b_rev  = year0_revenue(S_d_b, S_m_b, area, b_spray)
                b_liab = calculate_10yr_liability_v2(
                    S_d_b / area, S_m_b / area, p_dicot, p_mono, b_pen,
                    pre_d=Td_b / area, pre_m=Tm_b / area,
                    scenario=SCENARIO, recall=rv, cell_area_m2=area) * area_ha
                c_herb += pr * b_herb
                c_rev  += pr * b_rev
                c_liab += pr * b_liab
                if b_spray:
                    c_spr += pr * area
                c_vol  += pr * b_dose * area
            curve[rv] = {'herb': c_herb, 'rev': c_rev, 'liab': c_liab,
                         'spr_area': c_spr, 'vol': c_vol}

    # -- DUAL-DOSE BACKGROUND KILL SWEEP ----------------------------
    ddk_l = {}; ddk_r = {}; ddk_h = {}
    if ENABLE_DD_KILL_SWEEP and SCENARIO == 'SCENARIO_3_SECURE':
        for ddk in DD_KILL_SWEEP:
            kill_dd = full_kill_rate if dose == 1.0 else ddk
            S_d_dd  = T_d * (1.0 - kill_dd)
            S_m_dd  = T_m * (1.0 - kill_dd)
            ddk_r[ddk] = year0_revenue(S_d_dd, S_m_dd, area, True)
            ddk_h[ddk] = (herbicide_cost_fullrate_eur_per_m2 * area) * dose
            ddk_l[ddk] = calculate_10yr_liability_v2(
                S_d_dd / area, S_m_dd / area, p_dicot, p_mono, STUNTED_SEED_PENALTY,
                pre_d=bio_dens_d, pre_m=bio_dens_m,
                scenario=SCENARIO, recall=recall_dicot,
                dd_base_kill=ddk, cell_area_m2=area) * area_ha

    # -- CURRENT-YEAR REVENUE (SSWM) ----------------------------------------
    rev_sswm = year0_revenue(S_d, S_m, area, spray)

    # -- PATH B: BROADCAST ---------------------------------------------------
    bc_cost   = herbicide_cost_fullrate_eur_per_m2 * area
    dens_d_bc = T_d * (1.0 - full_kill_rate) / area
    dens_m_bc = T_m * (1.0 - full_kill_rate) / area

    _yc_b = [0.0] * SIM_YEARS
    shad_bc_base = calculate_10yr_liability_v2(
        dens_d_bc, dens_m_bc, p_dicot, p_mono, STUNTED_SEED_PENALTY,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        is_broadcast=True, recall=recall_dicot, cell_area_m2=area,
        year_costs_out=_yc_b) * area_ha
    shad_bc_low = calculate_10yr_liability_v2(
        dens_d_bc, dens_m_bc, p_dicot_low, p_mono_low, STUNTED_SEED_PENALTY,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        is_broadcast=True, recall=recall_dicot, cell_area_m2=area) * area_ha
    shad_bc_high = calculate_10yr_liability_v2(
        dens_d_bc, dens_m_bc, p_dicot_high, p_mono_high, STUNTED_SEED_PENALTY,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        is_broadcast=True, recall=recall_dicot, cell_area_m2=area) * area_ha

    split_bc = {}; split_rev_bc = {}
    for split_name, (patch_frac, drip_d_var, drip_m_var) in drip_rates.items():
        T_d_v = _true_count(P_d, recall_dicot,   patch_frac, drip_d_var, area)
        T_m_v = _true_count(P_m, recall_monocot, patch_frac, drip_m_var, area)
        split_rev_bc[split_name] = year0_revenue(
            T_d_v * (1.0 - full_kill_rate), T_m_v * (1.0 - full_kill_rate), area, True)
        split_bc[split_name] = calculate_10yr_liability_v2(
            T_d_v * (1.0 - full_kill_rate) / area,
            T_m_v * (1.0 - full_kill_rate) / area,
            p_dicot, p_mono, STUNTED_SEED_PENALTY,
            pre_d=T_d_v / area, pre_m=T_m_v / area,
            is_broadcast=True, recall=recall_dicot, cell_area_m2=area) * area_ha

    kill_bc = {}
    for k_name, k_params in KILL_SCENARIOS.items():
        kill_bc[k_name] = calculate_10yr_liability_v2(
            T_d * (1.0 - k_params['kill']) / area,
            T_m * (1.0 - k_params['kill']) / area,
            p_dicot, p_mono, k_params['stunt'],
            pre_d=bio_dens_d, pre_m=bio_dens_m,
            is_broadcast=True, recall=recall_dicot,
            stunt_future=k_params['stunt'],
            future_kill=k_params['kill'], cell_area_m2=area) * area_ha

    recall_bc = {}
    for rname, (rv, rdrip_d, rdrip_m) in recall_precomp.items():
        T_d_r = _true_count(P_d, rv, 0.70, rdrip_d, area)
        T_m_r = _true_count(P_m, rv, 0.70, rdrip_m, area)
        recall_bc[rname] = calculate_10yr_liability_v2(
            T_d_r * (1.0 - full_kill_rate) / area,
            T_m_r * (1.0 - full_kill_rate) / area,
            p_dicot, p_mono, STUNTED_SEED_PENALTY,
            pre_d=T_d_r / area, pre_m=T_m_r / area,
            is_broadcast=True, recall=rv, cell_area_m2=area) * area_ha

    rev_bc = year0_revenue(T_d * (1.0 - full_kill_rate),
                           T_m * (1.0 - full_kill_rate), area, True)

    if ENABLE_DD_KILL_SWEEP and SCENARIO == 'SCENARIO_3_SECURE':
        ddk_bc = {ddk: shad_bc_base for ddk in DD_KILL_SWEEP}
    else:
        ddk_bc = {}

    out.update(
        T_d=T_d, T_m=T_m, bio_dens_d=bio_dens_d, bio_dens_m=bio_dens_m,
        S_d=S_d, S_m=S_m, dens_d=dens_d, dens_m=dens_m,
        spray=spray, dose=dose, real_herb=real_herb,
        rev_sswm=rev_sswm, rev_bc=rev_bc, bc_cost=bc_cost,
        shad_sswm_base=shad_sswm_base, shad_sswm_low=shad_sswm_low,
        shad_sswm_high=shad_sswm_high,
        shad_bc_base=shad_bc_base, shad_bc_low=shad_bc_low,
        shad_bc_high=shad_bc_high,
        marg_sc_d=marg_sc_d, marg_sc_m=marg_sc_m,
        split_sswm=split_sswm, split_bc=split_bc,
        split_rev_sswm=split_rev_sswm, split_rev_bc=split_rev_bc,
        yc_sswm=[v * area_ha for v in _yc_s],
        yc_bc=[v * area_ha for v in _yc_b],
        kill_sswm=kill_sswm, kill_bc=kill_bc,
        recall_sswm=recall_sswm, recall_bc=recall_bc,
        sweep_l_sswm=sweep_l_sswm, sweep_r_sswm=sweep_r_sswm,
        sweep_l_bc=sweep_l_bc, sweep_r_bc=sweep_r_bc,
        ddk_l=ddk_l, ddk_r=ddk_r, ddk_h=ddk_h, ddk_bc=ddk_bc,
        curve=curve,
    )
    return out


# ================= MAIN EXECUTION =================
lyr = QgsProject.instance().mapLayersByName(GRID_LAYER_NAME)[0]

out_fields = [
    ('Area_m2',   QVariant.Double), ('P_Dicot',   QVariant.Double),
    ('P_Monocot', QVariant.Double), ('TDicot',    QVariant.Double),
    ('TMonocot',  QVariant.Double), ('TrueDenD',  QVariant.Double),
    ('TrueDenM',  QVariant.Double), ('SDicot',    QVariant.Double),
    ('SMonocot',  QVariant.Double), ('SurvDenD',  QVariant.Double),
    ('SurvDenM',  QVariant.Double), ('Spray_Dec', QVariant.String),
    ('Dose_Rate', QVariant.Double), ('Rev_Eur',   QVariant.Double),
    ('Herb_Cost', QVariant.Double), ('Fut_Liab',  QVariant.Double),
    ('Net_Ret',   QVariant.Double), ('Marg_SC_D', QVariant.Double),
    ('Marg_SC_M', QVariant.Double)
]
for nm, t in out_fields:
    ensure_field_exists(lyr, nm, t)

fi    = {f.name(): i for i, f in enumerate(lyr.fields())}
idx_d = lyr.fields().indexFromName(INPUT_DICOT_FIELD)
idx_m = lyr.fields().indexFromName(INPUT_MONOCOT_FIELD)

# -- single pass over the features for all field-level totals --------------
print("Loading features and computing field totals (single pass)...")
features_dict = {}
total_detected_d = 0.0
total_detected_m = 0.0
empty_cells_d = 0
empty_cells_m = 0
empty_cells_both = 0          # cells with NO detection of either class
total_area_m2 = 0.0           # field area for the background density

for f in lyr.getFeatures():
    features_dict[f.id()] = f
    total_area_m2 += f.geometry().area()
    a = f.attributes()
    pd_ = float(a[idx_d] or 0.0)
    pm_ = float(a[idx_m] or 0.0)
    total_detected_d += pd_
    total_detected_m += pm_
    if pd_ == 0:
        empty_cells_d += 1
    if pm_ == 0:
        empty_cells_m += 1
    if pd_ == 0 and pm_ == 0:
        empty_cells_both += 1

total_missed_d = (total_detected_d / recall_dicot   - total_detected_d) if recall_dicot   > 0 else 0.0
total_missed_m = (total_detected_m / recall_monocot - total_detected_m) if recall_monocot > 0 else 0.0

# drip_rates holds (patch_frac, iso_d, iso_m); under 'area' iso_* is a
# background density in weeds per m2, under 'empty_cells' weeds per empty cell.
drip_rates = {}
for split_name, patch_frac in SPLIT_SCENARIOS.items():
    drip_rates[split_name] = (
        patch_frac,
        _iso_param(total_detected_d, recall_dicot,   patch_frac, empty_cells_d),
        _iso_param(total_detected_m, recall_monocot, patch_frac, empty_cells_m),
    )

isolated_drip_d = drip_rates['split_70_30'][1]
isolated_drip_m = drip_rates['split_70_30'][2]

recall_precomp = {}
for rname, rv in RECALL_SCENARIOS.items():
    recall_precomp[rname] = (
        rv,
        _iso_param(total_detected_d, rv, 0.70, empty_cells_d),
        _iso_param(total_detected_m, rv, 0.70, empty_cells_m),
    )

sweep_precomp = {}
for rv in RECALL_SWEEP:
    sweep_precomp[rv] = (
        _iso_param(total_detected_d, rv, 0.70, empty_cells_d),
        _iso_param(total_detected_m, rv, 0.70, empty_cells_m),
    )

if ALLOCATION_RULE == 'area':
    print(f"Allocation: area | isolated background {isolated_drip_d:.4f} dicots/m2, "
          f"{isolated_drip_m:.4f} monocots/m2 over {total_area_m2:.0f} m2")
else:
    print(f"Allocation: empty_cells (validation) | {isolated_drip_d:.4f} dicots/cell | "
          f"{isolated_drip_m:.4f} monocots/cell")
print(f"Annuity factor: {ANNUITY_FACTOR:.4f}")
print(f"Cells: {len(features_dict)} | scenario={SCENARIO} | cache={'ON' if USE_CELL_CACHE else 'OFF'}")
print(f"Future detection: {FUTURE_DETECTION} | ET future rule: {ET_FUTURE_RULE} "
      f"| DD future cost: {DD_FUTURE_COST_MODE} | DD_BASE_KILL: {DD_BASE_KILL}")
print(f"Survivor conditioning: {SURVIVOR_CONDITIONING} | recall curve: "
      f"{'ON (' + str(len(RECALL_CURVE)) + ' values)' if ENABLE_RECALL_CURVE else 'off'} "
      f"| validation recall sweep: {'ON' if ENABLE_RECALL_SWEEP else 'off'}")
if SURVIVOR_CONDITIONING == 'v16':
    print("WARNING: SURVIVOR_CONDITIONING = 'v16' is a validation setting. The recall "
          "curve does not converge at r = 1 under it.")


def equilibrium_check(start_density=13.4, years=100, verbose=True):
    """CALIBRATION DIAGNOSTIC - printed at the start of every run.

    Projects the weed population forward under CONTINUOUS BROADCAST application
    at the baseline kill rate, which is the management these fields have
    actually received, and reports where it settles. If the model is calibrated,
    that equilibrium should sit close to the observed field densities. If it
    runs away, the ten-year projection is not representing the system and the
    seed-set parameter needs revisiting.
    """
    pdx = p_dicot
    bank = start_density / pdx['emergence'] * pdx['bank_survival']
    for _ in range(years):
        em   = bank * pdx['emergence']
        surv = em * (1.0 - full_kill_rate)
        seeds = (((pdx['max_seeds'] * surv) / (1.0 + pdx['intra_comp'] * surv))
                 * pdx['overwinter_survival'] * STUNTED_SEED_PENALTY)
        bank = bank * pdx['bank_survival'] + seeds
    eq = bank * pdx['emergence']
    if verbose:
        print(f"[CALIBRATION] seed set of survivors (gamma) = {STUNTED_SEED_PENALTY}")
        print(f"[CALIBRATION] equilibrium under continuous broadcast: "
              f"{eq:.1f} plants m-2")
        if eq > 60:
            print("              WARNING: the population runs away even under "
                  "broadcast control. The projection is not representing the "
                  "system; lower STUNTED_SEED_PENALTY.")
        elif eq < 1:
            print("              WARNING: the population is driven to extinction "
                  "under broadcast control. Raise STUNTED_SEED_PENALTY.")
    return eq


equilibrium_check()
print(f"Starting simulation (model v{MODEL_VERSION}) ...")

# Accumulators
field_area_m2 = 0.0;  sprayed_area_m2 = 0.0
vol_herb_sswm = 0.0;  vol_herb_bc     = 0.0
sum_rev_sswm  = 0.0;  sum_herb_sswm   = 0.0
sum_rev_bc    = 0.0;  sum_herb_bc     = 0.0

sum_shad_sswm_base = 0.0; sum_shad_bc_base = 0.0
sum_yc_sswm = [0.0] * SIM_YEARS; sum_yc_bc = [0.0] * SIM_YEARS
sum_shad_sswm_low  = 0.0; sum_shad_bc_low  = 0.0
sum_shad_sswm_high = 0.0; sum_shad_bc_high = 0.0

sum_split_sswm  = {n: 0.0 for n in SPLIT_SCENARIOS}
sum_split_bc    = {n: 0.0 for n in SPLIT_SCENARIOS}
sum_split_rev_sswm = {n: 0.0 for n in SPLIT_SCENARIOS}
sum_split_rev_bc   = {n: 0.0 for n in SPLIT_SCENARIOS}
sum_kill_sswm   = {n: 0.0 for n in KILL_SCENARIOS}
sum_kill_bc     = {n: 0.0 for n in KILL_SCENARIOS}
sum_recall_sswm = {n: 0.0 for n in RECALL_SCENARIOS}
sum_recall_bc   = {n: 0.0 for n in RECALL_SCENARIOS}

sweep_liab_sswm = {r: 0.0 for r in RECALL_SWEEP}
sweep_liab_bc   = {r: 0.0 for r in RECALL_SWEEP}
sweep_rev_sswm  = {r: 0.0 for r in RECALL_SWEEP}
sweep_rev_bc    = {r: 0.0 for r in RECALL_SWEEP}

ddk_liab_sswm = {k: 0.0 for k in DD_KILL_SWEEP}
ddk_liab_bc   = {k: 0.0 for k in DD_KILL_SWEEP}
ddk_rev_sswm  = {k: 0.0 for k in DD_KILL_SWEEP}
ddk_herb_sswm = {k: 0.0 for k in DD_KILL_SWEEP}

# recall curve accumulators
curve_sum = {rv: {'herb': 0.0, 'rev': 0.0, 'liab': 0.0, 'spr_area': 0.0, 'vol': 0.0}
             for rv in RECALL_CURVE}

bulk_updates = {}

# =====================================================================
# MAIN CELL LOOP: compute (cached) then accumulate
# =====================================================================
for feat_id, feat in features_dict.items():
    area = feat.geometry().area()
    if area <= 0:
        continue

    P_d = float(feat.attributes()[idx_d] or 0.0)
    P_m = float(feat.attributes()[idx_m] or 0.0)

    c = compute_cell(P_d, P_m, area)

    field_area_m2 += area
    if c['spray']:
        sprayed_area_m2 += area
    vol_herb_sswm += c['dose'] * area
    vol_herb_bc   += area

    sum_rev_sswm  += c['rev_sswm']
    sum_herb_sswm += c['real_herb']
    sum_rev_bc    += c['rev_bc']
    sum_herb_bc   += c['bc_cost']

    sum_shad_sswm_base += c['shad_sswm_base']
    sum_shad_sswm_low  += c['shad_sswm_low']
    sum_shad_sswm_high += c['shad_sswm_high']
    sum_shad_bc_base   += c['shad_bc_base']
    for _t in range(SIM_YEARS):
        sum_yc_sswm[_t] += c['yc_sswm'][_t]
        sum_yc_bc[_t]   += c['yc_bc'][_t]
    sum_shad_bc_low    += c['shad_bc_low']
    sum_shad_bc_high   += c['shad_bc_high']

    for n in SPLIT_SCENARIOS:
        sum_split_sswm[n] += c['split_sswm'][n]
        sum_split_bc[n]   += c['split_bc'][n]
        sum_split_rev_sswm[n] += c['split_rev_sswm'][n]
        sum_split_rev_bc[n]   += c['split_rev_bc'][n]
    for n in KILL_SCENARIOS:
        sum_kill_sswm[n] += c['kill_sswm'][n]
        sum_kill_bc[n]   += c['kill_bc'][n]
    for n in RECALL_SCENARIOS:
        sum_recall_sswm[n] += c['recall_sswm'][n]
        sum_recall_bc[n]   += c['recall_bc'][n]

    if ENABLE_RECALL_SWEEP:
        for rv in RECALL_SWEEP:
            sweep_liab_sswm[rv] += c['sweep_l_sswm'][rv]
            sweep_rev_sswm[rv]  += c['sweep_r_sswm'][rv]
            sweep_liab_bc[rv]   += c['sweep_l_bc'][rv]
            sweep_rev_bc[rv]    += c['sweep_r_bc'][rv]

    if ENABLE_DD_KILL_SWEEP and SCENARIO == 'SCENARIO_3_SECURE':
        for k in DD_KILL_SWEEP:
            ddk_liab_sswm[k] += c['ddk_l'][k]
            ddk_rev_sswm[k]  += c['ddk_r'][k]
            ddk_herb_sswm[k] += c['ddk_h'][k]
            ddk_liab_bc[k]   += c['ddk_bc'][k]

    if c['curve']:
        for rv, v in c['curve'].items():
            cs = curve_sum[rv]
            for kk in ('herb', 'rev', 'liab', 'spr_area', 'vol'):
                cs[kk] += v[kk]

    bulk_updates[feat_id] = {
        fi['Area_m2']:   area,             fi['P_Dicot']:   P_d,
        fi['P_Monocot']: P_m,              fi['TDicot']:    c['T_d'],
        fi['TMonocot']:  c['T_m'],         fi['TrueDenD']:  c['bio_dens_d'],
        fi['TrueDenM']:  c['bio_dens_m'],  fi['SDicot']:    c['S_d'],
        fi['SMonocot']:  c['S_m'],         fi['SurvDenD']:  c['dens_d'],
        fi['SurvDenM']:  c['dens_m'],
        fi['Spray_Dec']: 'Spray' if c['spray'] else 'No',
        fi['Dose_Rate']: c['dose'],        fi['Rev_Eur']:   c['rev_sswm'],
        fi['Herb_Cost']: c['real_herb'],   fi['Fut_Liab']:  c['shad_sswm_base'],
        fi['Net_Ret']:   (c['rev_sswm'] - c['real_herb'] - c['shad_sswm_base']),
        fi['Marg_SC_D']: c['marg_sc_d'],   fi['Marg_SC_M']: c['marg_sc_m']
    }

print("Writing results to map...")
lyr.dataProvider().changeAttributeValues(bulk_updates)
lyr.triggerRepaint()
iface.mapCanvas().refresh()

hits = _CACHE_STATS['hit']; misses = _CACHE_STATS['miss']; tot = hits + misses
print(f"cell cache: {hits} hits / {misses} full computations "
      f"({100.0*hits/tot if tot else 0:.2f}% hit rate)")
if FUTURE_DETECTION == 'probabilistic':
    print(f"max simultaneous trajectories in any cell-year: "
          f"{_BRANCH_STATS['max_states']} (ceiling {MAX_STATES})")
    if _BRANCH_STATS['truncated'] > 0:
        print("       WARNING: hit the ceiling. Raise MAX_STATES or lower "
              "MERGE_DIGITS if you want a finer expectation.")

# ================= SUMMARY =================
annuity_factor = ANNUITY_FACTOR
field_area_ha  = field_area_m2 / 10000.0


def calc_extra_sc(shad_sswm, shad_bc):
    """Annual-equivalent EXTRA liability of SSWM relative to broadcast.
    The NPV difference is converted with the capital recovery factor
    (1/ANNUITY_FACTOR), Eq. 12."""
    return ((shad_sswm - shad_bc) / field_area_ha) / annuity_factor


sc_extra_base = calc_extra_sc(sum_shad_sswm_base, sum_shad_bc_base)

# per-year expected costs (undiscounted, EUR/ha) and a self-check that
# discounting them at 4% over 10 years reproduces the baseline FWP exactly.
yearly_rows = []
for _t in range(SIM_YEARS):
    _s = sum_yc_sswm[_t] / field_area_ha
    _b = sum_yc_bc[_t] / field_area_ha
    yearly_rows.append([_t + 1, _s, _b, _s - _b])
_fwp_chk = sum(r[3] / (1.0 + DISCOUNT_RATE) ** r[0] for r in yearly_rows) / annuity_factor
print(f"yearly export self-check: FWP from yearly file {_fwp_chk:.4f} vs "
      f"baseline {sc_extra_base:.4f} EUR/ha "
      f"({'OK' if abs(_fwp_chk - sc_extra_base) < 1e-6 else 'MISMATCH'})")
sc_extra_low  = calc_extra_sc(sum_shad_sswm_low,  sum_shad_bc_low)
sc_extra_high = calc_extra_sc(sum_shad_sswm_high, sum_shad_bc_high)

net_ann_sswm  = ((sum_rev_sswm / field_area_ha) - (sum_herb_sswm / field_area_ha)
                 - sc_extra_base - sswm_annuity_per_ha)
net_ann_bc    = (sum_rev_bc / field_area_ha) - (sum_herb_bc / field_area_ha)
net_sswm_low  = ((sum_rev_sswm / field_area_ha) - (sum_herb_sswm / field_area_ha)
                 - sc_extra_low  - sswm_annuity_per_ha)
net_sswm_high = ((sum_rev_sswm / field_area_ha) - (sum_herb_sswm / field_area_ha)
                 - sc_extra_high - sswm_annuity_per_ha)

sprayed_pct    = (sprayed_area_m2 / field_area_m2) * 100 if field_area_m2 > 0 else 0
vol_reduct_pct = (1.0 - vol_herb_sswm / vol_herb_bc) * 100 if vol_herb_bc > 0 else 0

# -- PARTIAL BUDGET DECOMPOSITION -----------------------------------
# The four components of the difference vs broadcast (Fig. 2).
d_herb_saved = (sum_herb_bc - sum_herb_sswm) / field_area_ha   # (a) positive
d_yield_yr0  = (sum_rev_sswm - sum_rev_bc)  / field_area_ha    # (b) negative
d_liab       = -sc_extra_base                                  # (c) negative
d_tech       = -sswm_annuity_per_ha                            # (d) negative
d_net_check  = d_herb_saved + d_yield_yr0 + d_liab + d_tech    # (e)

# Break-even technology cost — needs NO extra simulation.
tech_breakeven = d_herb_saved + d_yield_yr0 + d_liab

print(f"\n{'='*70}")
print(f"RESULTS (model v{MODEL_VERSION}): {SCENARIO}  |  Tech: {TECH_SCENARIO}")
print(f"Field: {field_area_ha:.2f} ha  |  Sprayed: {sprayed_pct:.1f}%")
print(f"{'='*70}")
print(f"{'METRIC':<24} | {'SSWM':>10} | {'BROADCAST':>10} | {'DIFF':>9}")
print(f"{'-'*70}")
print(f"{'Revenue (EUR/ha)':<24} | {sum_rev_sswm/field_area_ha:>10.2f} | {sum_rev_bc/field_area_ha:>10.2f} | {d_yield_yr0:>+9.2f}")
print(f"{'Herbicide (EUR/ha)':<24} | {sum_herb_sswm/field_area_ha:>10.2f} | {sum_herb_bc/field_area_ha:>10.2f} | {d_herb_saved:>+9.2f}")
print(f"{'Future Liability':<24} | {sc_extra_base:>10.2f} | {'0.00':>10} | {d_liab:>+9.2f}")
print(f"{'Tech Annuity':<24} | {sswm_annuity_per_ha:>10.2f} | {'0.00':>10} | {d_tech:>+9.2f}")
print(f"{'-'*70}")
print(f"{'NET RETURN (EUR/ha)':<24} | {net_ann_sswm:>10.2f} | {net_ann_bc:>10.2f} | {net_ann_sswm-net_ann_bc:>+9.2f}")
print(f"{'='*70}")
print(f"Herbicide volume reduced by {vol_reduct_pct:.1f}%")

print(f"\n===== PARTIAL BUDGET DECOMPOSITION (EUR/ha vs broadcast) =====")
print(f"  (a) Herbicide cost saved      : {d_herb_saved:>+8.2f}")
print(f"  (b) Year-0 yield loss         : {d_yield_yr0:>+8.2f}")
print(f"  (c) Future weed penalty       : {d_liab:>+8.2f}")
print(f"  (d) Technology cost           : {d_tech:>+8.2f}")
print(f"  {'-'*44}")
print(f"  (e) NET vs broadcast          : {d_net_check:>+8.2f}")
print(f"\n  Break-even technology cost    : {tech_breakeven:>8.2f} EUR/ha")
print(f"  (max annualised tech+info cost at which SSWM matches broadcast)")
print(f"{'='*70}\n")

print("===== SENSITIVITY 1: Weed Population Dynamics (+/-20%) =====")
print(f"  Optimistic: Liability={sc_extra_low:>7.2f}  Net={net_sswm_low:>8.2f}")
print(f"  Baseline:   Liability={sc_extra_base:>7.2f}  Net={net_ann_sswm:>8.2f}")
print(f"  Pessimist:  Liability={sc_extra_high:>7.2f}  Net={net_sswm_high:>8.2f}")
print(f"{'='*70}\n")

sc_split = {n: calc_extra_sc(sum_split_sswm[n], sum_split_bc[n]) for n in SPLIT_SCENARIOS}
sc_kill  = {n: calc_extra_sc(sum_kill_sswm[n],  sum_kill_bc[n])  for n in KILL_SCENARIOS}


def net_from_sc(sc):
    return (sum_rev_sswm/field_area_ha) - (sum_herb_sswm/field_area_ha) - sc - sswm_annuity_per_ha


def net_split(n):
    """Full-channel allocation sensitivity. Year-0 revenue of SSWM and
    of broadcast are both recomputed under allocation n; the value is expressed
    so that (net_split - net_ann_bc) is the difference vs broadcast on the same
    allocation, which is how the R script uses Split_*_Net."""
    d = ((sum_split_rev_sswm[n] - sum_split_rev_bc[n])
         - (sum_herb_sswm - sum_herb_bc)) / field_area_ha
    return net_ann_bc + d - sc_split[n] - sswm_annuity_per_ha


print("===== SENSITIVITY 2: False Negative Spatial Allocation (full channel) =====")
print(f"  allocation rule: {ALLOCATION_RULE}")
for sn in SPLIT_SCENARIOS:
    print(f"  {sn:<14}: Liability={sc_split[sn]:>7.2f}  Net={net_split(sn):>8.2f}  "
          f"(vs broadcast {net_split(sn) - net_ann_bc:>+8.2f})")
_chk = abs(net_split('split_70_30') - net_ann_sswm)
print(f"  self-check 70/30 vs baseline: {_chk:.6f} EUR/ha ({'OK' if _chk < 1e-6 else 'MISMATCH'})")
print(f"{'='*70}\n")

print("===== SENSITIVITY 3: Herbicide Efficacy =====")
for kn, kp in KILL_SCENARIOS.items():
    print(f"  {kn:<10} (k={kp['kill']}, stunt={kp['stunt']}): "
          f"Liability={sc_kill[kn]:>7.2f}  Net={net_from_sc(sc_kill[kn]):>8.2f}")
print(f"{'='*70}\n")

sc_recall = {n: calc_extra_sc(sum_recall_sswm[n], sum_recall_bc[n]) for n in RECALL_SCENARIOS}

print("===== SENSITIVITY 4: Sensor Recall (+/-20%, LIABILITY CHANNEL ONLY) =====")
print("  Varies the future liability only and holds year-0 revenue at the")
print("  measured recall. The reported recall sensitivity comes from the recall")
print("  curve below, which also re-derives the year-0 spray map.")
for rn, rv in RECALL_SCENARIOS.items():
    print(f"  {rn} (recall={rv}, FN={round((1-rv)*100,1)}%): "
          f"Liability={sc_recall[rn]:>7.2f}  Net={net_from_sc(sc_recall[rn]):>8.2f}")
print(f"{'='*70}\n")

sweep_rows = []
if ENABLE_RECALL_SWEEP:
    print("===== RECALL SWEEP: full effect (liability + Year-0 yield) =====")
    print("  Spray map held FIXED at observed detections => LOWER bound on the")
    print("  value of better detection, so the break-even recall shown is an")
    print("  UPPER bound on what is actually required.")
    print(f"  {'recall':>7} | {'FN%':>5} | {'d_herb':>8} | {'d_yield':>8} | {'d_liab':>8} | {'NET':>8}")
    print(f"  {'-'*58}")
    for rv in RECALL_SWEEP:
        d_h = d_herb_saved
        d_y = (sweep_rev_sswm[rv] - sweep_rev_bc[rv]) / field_area_ha
        d_l = -calc_extra_sc(sweep_liab_sswm[rv], sweep_liab_bc[rv])
        net = d_h + d_y + d_l - sswm_annuity_per_ha
        sweep_rows.append((rv, d_h, d_y, d_l, net))
        print(f"  {rv:>7.2f} | {100*(1-rv):>5.1f} | {d_h:>+8.2f} | {d_y:>+8.2f} | "
              f"{d_l:>+8.2f} | {net:>+8.2f}")
    print(f"{'='*70}\n")

ddk_rows = []
if ENABLE_DD_KILL_SWEEP and SCENARIO == 'SCENARIO_3_SECURE':
    print("===== DUAL-DOSE BACKGROUND KILL SWEEP =====")
    print("  Tests the assumption that a 50% background dose controls weeds as")
    print("  well as a full dose. Under 'v13' the k_bg = 0.95 differential is")
    print("  zero by construction; under 'dose_scaled' it keeps the future")
    print(f"  chemical saving. Mode: {DD_FUTURE_COST_MODE}. SELF-CHECK: the")
    print(f"  k_bg = {DD_BASE_KILL:.2f} row must equal the baseline of this run.")
    print(f"  {'k_bg':>6} | {'d_herb':>8} | {'d_yield':>8} | {'d_liab':>8} | {'NET':>8}")
    print(f"  {'-'*50}")
    for ddk in DD_KILL_SWEEP:
        d_h = (sum_herb_bc - ddk_herb_sswm[ddk]) / field_area_ha
        d_y = (ddk_rev_sswm[ddk] - sum_rev_bc) / field_area_ha
        d_l = -calc_extra_sc(ddk_liab_sswm[ddk], ddk_liab_bc[ddk])
        net = d_h + d_y + d_l - sswm_annuity_per_ha
        ddk_rows.append((ddk, d_h, d_y, d_l, net))
        print(f"  {ddk:>6.2f} | {d_h:>+8.2f} | {d_y:>+8.2f} | {d_l:>+8.2f} | {net:>+8.2f}")
    print(f"{'='*70}\n")

# ================= RECALL CURVE SUMMARY =================
curve_rows = []
curve_breakeven_recall = None
if ENABLE_RECALL_CURVE and SCENARIO in ('SCENARIO_1_ECON', 'SCENARIO_2_SPOT',
                                       'SCENARIO_3_SECURE'):
    print("===== RECALL CURVE: sensor of recall r' on the same population =====")
    print("  Year-0 map re-derived at each r'; years 1-10 at r'; broadcast unchanged.")
    print("  SELF-CHECK: the row at the measured recall must equal the baseline")
    print("  decomposition above (D_Herb, D_Yield, D_Liab).")
    print(f"  {'recall':>7} | {'spray%':>6} | {'volRed%':>7} | {'d_herb':>8} | "
          f"{'d_yield':>8} | {'d_liab':>8} | {'BE_tech':>8}")
    print(f"  {'-'*70}")
    for rv in RECALL_CURVE:
        cs  = curve_sum[rv]
        d_h = (sum_herb_bc - cs['herb']) / field_area_ha
        d_y = (cs['rev'] - sum_rev_bc) / field_area_ha
        d_l = -calc_extra_sc(cs['liab'], sum_shad_bc_base)
        net = d_h + d_y + d_l                    # = break-even tech cost at r'
        spr = 100.0 * cs['spr_area'] / field_area_m2 if field_area_m2 > 0 else 0.0
        vol = 100.0 * (1.0 - cs['vol'] / vol_herb_bc) if vol_herb_bc > 0 else 0.0
        curve_rows.append((rv, d_h, d_y, d_l, net, spr, vol))
        flag = '  <- baseline' if abs(rv - recall_dicot) < 1e-12 else ''
        print(f"  {rv:>7.3f} | {spr:>6.1f} | {vol:>7.1f} | {d_h:>+8.2f} | "
              f"{d_y:>+8.2f} | {d_l:>+8.2f} | {net:>+8.2f}{flag}")
    base_row = [r for r in curve_rows if abs(r[0] - recall_dicot) < 1e-12]
    if base_row:
        dev = max(abs(base_row[0][1] - d_herb_saved), abs(base_row[0][2] - d_yield_yr0),
                  abs(base_row[0][3] - d_liab))
        print(f"  self-check deviation at measured recall: {dev:.6f} EUR/ha "
              f"({'OK' if dev < 1e-6 else 'MISMATCH - investigate'})")
    # break-even recall at ZERO technology cost, linear interpolation on the grid
    # (R adds the five cost scenarios).
    rs  = [r[0] for r in curve_rows]; ns = [r[4] for r in curve_rows]
    for i in range(1, len(rs)):
        if ns[i - 1] < 0.0 <= ns[i]:
            curve_breakeven_recall = rs[i - 1] + (0.0 - ns[i - 1]) * (rs[i] - rs[i - 1]) / (ns[i] - ns[i - 1])
            break
    if ns and ns[0] >= 0.0:
        curve_breakeven_recall = rs[0]
    print(f"  Break-even recall at zero technology cost: "
          f"{'none below %.2f' % rs[-1] if curve_breakeven_recall is None else '%.3f' % curve_breakeven_recall}")
    print(f"{'='*70}\n")

# ================= CSV EXPORT =================
project_name = QgsProject.instance().baseName() or "Untitled_Project"

csv_row = [
    project_name, SCENARIO, round(field_area_ha, 2), round(sprayed_pct, 2),
    round(sum_rev_sswm  / field_area_ha, 2), round(sum_herb_sswm / field_area_ha, 2),
    round(sc_extra_base, 2), round(sswm_annuity_per_ha, 2), round(net_ann_sswm, 2),
    round(sum_rev_bc  / field_area_ha, 2), round(sum_herb_bc / field_area_ha, 2),
    0.00, 0.00, round(net_ann_bc, 2),
    round(vol_reduct_pct, 2),
    round(sc_extra_low, 2),  round(net_sswm_low, 2),
    round(sc_extra_high, 2), round(net_sswm_high, 2),
    round(sc_split['split_90_10'], 2), round(net_split('split_90_10'), 2),
    round(sc_split['split_70_30'], 2), round(net_split('split_70_30'), 2),
    round(sc_split['split_50_50'], 2), round(net_split('split_50_50'), 2),
    round(sc_kill['k_low'],  2), round(net_from_sc(sc_kill['k_low']),  2),
    round(sc_kill['k_base'], 2), round(net_from_sc(sc_kill['k_base']), 2),
    round(sc_kill['k_high'], 2), round(net_from_sc(sc_kill['k_high']), 2),
    round(sc_recall['recall_low'],  2), round(net_from_sc(sc_recall['recall_low']),  2),
    round(sc_recall['recall_base'], 2), round(net_from_sc(sc_recall['recall_base']), 2),
    round(sc_recall['recall_high'], 2), round(net_from_sc(sc_recall['recall_high']), 2),
    round(d_herb_saved, 2), round(d_yield_yr0, 2), round(d_liab, 2),
    round(d_tech, 2), round(d_net_check, 2), round(tech_breakeven, 2),
    # run configuration, so a results file always records the model
    # settings that produced it. Prevents mixing runs from different modes.
    FUTURE_DETECTION, ET_FUTURE_RULE, DD_FUTURE_COST_MODE,
    DD_BASE_KILL, recall_dicot, round(ANNUITY_FACTOR, 4),
    SURVIVOR_CONDITIONING,
    ALLOCATION_RULE,
    # descriptive statistics for the supplementary field table.
    len(features_dict),
    round(total_detected_d / field_area_ha, 1),
    round(total_detected_m / field_area_ha, 1),
    round(100.0 * (1.0 - empty_cells_both / len(features_dict)), 2),
]

CSV_HEADER = [
    'Project', 'Scenario', 'Area_Ha', 'Sprayed_Pct',
    'SSWM_Rev', 'SSWM_Herb', 'SSWM_Liab_Base', 'SSWM_Tech', 'SSWM_Net_Base',
    'BC_Rev', 'BC_Herb', 'BC_Liab', 'BC_Tech', 'BC_Net', 'Vol_Red_Pct',
    'Weed_Low_Liab', 'Weed_Low_Net', 'Weed_High_Liab', 'Weed_High_Net',
    'Split_90_10_Liab', 'Split_90_10_Net',
    'Split_70_30_Liab', 'Split_70_30_Net',
    'Split_50_50_Liab', 'Split_50_50_Net',
    'Kill_Low_Liab', 'Kill_Low_Net',
    'Kill_Base_Liab', 'Kill_Base_Net',
    'Kill_High_Liab', 'Kill_High_Net',
    'Recall_Low_Liab',  'Recall_Low_Net',
    'Recall_Base_Liab', 'Recall_Base_Net',
    'Recall_High_Liab', 'Recall_High_Net',
    'D_Herb_Saved', 'D_Yield_Yr0', 'D_Liab', 'D_Tech', 'D_Net', 'Tech_Breakeven',
    'Future_Detection', 'ET_Future_Rule', 'DD_Future_Cost', 'DD_Base_Kill',
    'Recall', 'Annuity_Factor',
    'Survivor_Cond',
    'Alloc_Rule',
    'N_Cells', 'Det_Dicot_Per_Ha', 'Det_Mono_Per_Ha', 'Pct_Cells_Infested',
]

file_exists = os.path.isfile(RESULTS_CSV_PATH)
try:
    os.makedirs(os.path.dirname(RESULTS_CSV_PATH), exist_ok=True)
    with open(RESULTS_CSV_PATH, 'a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(CSV_HEADER)
        writer.writerow(csv_row)
    print(f"[SUCCESS] Results appended to:\n{RESULTS_CSV_PATH}")
except Exception as e:
    print(f"[ERROR] Could not write to CSV: {e}")


def _append_rows(path, header, rows, prefix):
    exists = os.path.isfile(path)
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'a', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            if not exists:
                w.writerow(header)
            for r in rows:
                w.writerow(list(prefix) + [round(v, 4) for v in r])
        print(f"[SUCCESS] appended {len(rows)} rows to {path}")
    except Exception as exc:
        print(f"[ERROR] could not write {path}: {exc}")


if ENABLE_RECALL_SWEEP and sweep_rows:
    _append_rows(
        RESULTS_CSV_PATH.replace('.csv', '_recall_sweep.csv'),
        ['Project', 'Scenario', 'Area_Ha', 'Recall', 'D_Herb', 'D_Yield', 'D_Liab', 'Net'],
        sweep_rows, (project_name, SCENARIO, round(field_area_ha, 2)))

if ENABLE_DD_KILL_SWEEP and ddk_rows:
    _append_rows(
        RESULTS_CSV_PATH.replace('.csv', '_dd_kill_sweep.csv'),
        ['Project', 'Scenario', 'Area_Ha', 'K_Background', 'D_Herb', 'D_Yield', 'D_Liab', 'Net'],
        ddk_rows, (project_name, SCENARIO, round(field_area_ha, 2)))

if ENABLE_YEARLY_EXPORT and yearly_rows:
    _append_rows(
        RESULTS_CSV_PATH.replace('.csv', '_yearly.csv'),
        ['Project', 'Scenario', 'Area_Ha', 'Year', 'SSWM_Cost', 'BC_Cost', 'Extra_Cost'],
        yearly_rows, (project_name, SCENARIO, round(field_area_ha, 2)))

if ENABLE_RECALL_CURVE and curve_rows:
    # 'Net' is the difference vs broadcast BEFORE technology cost, i.e.
    # the break-even technology cost at that recall. R subtracts the cost
    # scenarios and locates the break-even recall.
    _append_rows(
        RESULTS_CSV_PATH.replace('.csv', '_recall_curve.csv'),
        ['Project', 'Scenario', 'Area_Ha', 'Recall', 'D_Herb', 'D_Yield', 'D_Liab',
         'Net', 'Sprayed_Pct', 'Vol_Red_Pct'],
        curve_rows, (project_name, SCENARIO, round(field_area_ha, 2)))


# ================= PER-CELL EXPORT (marginal shadow costs) =================
# The per-cell results written to the grid layer, as one CSV per run. By default
# only the ET strategy at 1 x 1 m is exported; these files are the input of the
# marginal shadow-cost figure (Supplementary Code 4).
CELL_EXPORT_COLUMNS = ['Area_m2', 'P_Dicot', 'P_Monocot', 'TDicot', 'TMonocot',
                       'TrueDenD', 'TrueDenM', 'SurvDenD', 'SurvDenM', 'Spray_Dec',
                       'Dose_Rate', 'Rev_Eur', 'Herb_Cost', 'Fut_Liab', 'Net_Ret',
                       'Marg_SC_D', 'Marg_SC_M']


def _nominal_cell_size(feats):
    """Side length (m) of a full grid cell, from the largest cell area in the
    layer (edge cells are clipped and therefore smaller)."""
    max_area = max((f.geometry().area() for f in feats.values()), default=0.0)
    return math.sqrt(max_area)


if ENABLE_CELL_EXPORT and SCENARIO in CELL_EXPORT_SCENARIOS:
    cell_size = _nominal_cell_size(features_dict)
    if (CELL_EXPORT_CELL_SIZE_M is None
            or abs(cell_size - CELL_EXPORT_CELL_SIZE_M) <= 0.05 * CELL_EXPORT_CELL_SIZE_M):
        cell_dir  = os.path.join(os.path.dirname(RESULTS_CSV_PATH), 'cells')
        cell_path = os.path.join(cell_dir, f"{project_name}_{SCENARIO}.csv")
        try:
            os.makedirs(cell_dir, exist_ok=True)
            with open(cell_path, 'w', newline='', encoding='utf-8') as fh:
                w = csv.writer(fh)
                w.writerow(['Project', 'Scenario', 'Cell_ID'] + CELL_EXPORT_COLUMNS)
                for feat_id, vals in bulk_updates.items():
                    w.writerow([project_name, SCENARIO, feat_id] +
                               [vals[fi[c]] for c in CELL_EXPORT_COLUMNS])
            print(f"[SUCCESS] per-cell results ({len(bulk_updates)} cells) written to {cell_path}")
        except Exception as exc:
            print(f"[ERROR] could not write {cell_path}: {exc}")
    else:
        print(f"Per-cell export skipped: cell size {cell_size:.2f} m "
              f"(CELL_EXPORT_CELL_SIZE_M = {CELL_EXPORT_CELL_SIZE_M}).")
