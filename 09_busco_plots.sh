#!/bin/bash
# =============================================================================
# SLURM JOB: BUSCO COMPLETENESS, KARYOTYPE AND SYNTENY PLOTS
#            (PTARMIGAN-SCAFFOLDED ASSEMBLIES)
# Step 09 — requires
# 07_busco_array.sh, and 08_ptarmigan_reference_busco.sh if you want
# the ptarmigan reference itself drawn alongside the grouse.
#
# NOT an array job. This is an aggregate step: one task reads every BUSCO run
# at once, because the completeness barplots compare assemblies against each
# other and the synteny plots compare them pairwise.
#
# THIS REUSES 09_busco_plots.py UNCHANGED. The Python script takes its input and
# output directories as arguments, so there is exactly one copy of the plotting
# logic for both trees — the undocumented BUSCO-Plot-Py workarounds in it
# (karyotype column requirements, the chromosome-ordering fix, the lineage-name
# correction, the SVG hook) are fixed in one place and cannot drift between the
# chicken and ptarmigan figures. Keep 09_busco_plots.py in the submission
# directory; this script will not run without it.
#
# ONE EDIT WAS MADE TO 09_busco_plots.py for this step: the species order now
# appends any species it discovers that is not in its built-in list, instead of
# silently dropping it. Without that, a LAGMUT_*_ref reference would have been
# left out of the per-species barplot and out of the default synteny pairs — the
# very comparison this step exists for. With only the original four species
# present the output is byte-identical to before, so figures made without a
# reference genome included are reproduced unchanged.
#
# WHAT TO LOOK AT FIRST
#   karyotype/ — the ptarmigan-ordered grouse assemblies. If the chr6/chr8
#     fusion is ancestral to Tetraoninae, ptarmigan should place these as ONE
#     chromosome, so the grouse chr_6 should now be ~65 Mb WITHOUT a stranded
#     ~0.1-1.2 Mb chr_8 beside it. That stranded fragment under chicken was the
#     clearest sign the reference, not the assembly, was imposing the split.
#   synteny/ — grouse vs the ptarmigan reference should be a far cleaner
#     diagonal than grouse vs chicken was. Where it is NOT clean is where the
#     real rearrangements are.
#
# Plots are produced with BUSCO-Plot-Py:
#   https://github.com/lorenzo-arcioni/BUSCO-Plot-Py
#
# OUTPUTS, under ${QC_DIR}/plots:
#   completeness/busco_barplot_<SPECIES>_completeness.png   one per species
#   completeness/busco_barplot_ALL_completeness.png         all assemblies
#   karyotype/<SPECIES>_<SAMPLE>_<HAP>.png                  one per assembly
#   synteny/<A>__vs__<B>.png                                one per pair
#
# Every PNG is written with a matching .svg (vector, for figures that will be
# resized or edited in Illustrator/Inkscape). Set SAVE_SVG=false below for
# PNG only.
#
# USAGE:
#   sbatch 09_busco_plots.sh
#
#   Options are set in the USER-DEFINED VARIABLES block below: which synteny
#   pairs to draw, plot orientation, how many chromosomes to show, and which
#   plot types to skip on a rerun.
#
# INSTALL NOTE (why this does not pip-install the package):
#   `pip install buscoplotpy` gives PyPI 0.0.2, which has NO synteny module —
#   only the barplot and karyoplot. Synteny exists only in the GitHub tree.
#   `pip install git+https://github.com/lorenzo-arcioni/BUSCO-Plot-Py` also
#   fails: the repo has no package-discovery configuration and a top-level
#   images/ directory, so setuptools refuses with "Multiple top-level packages
#   discovered in a flat-layout". The package is pure Python with no build
#   step, so this script clones it and puts the clone on PYTHONPATH, and
#   installs only its three real dependencies (matplotlib, seaborn, pandas).
# =============================================================================
#SBATCH --job-name=grouse_buscoplots_ptarmigan
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
# 32G: plotting is I/O- and matplotlib-bound, not compute-bound. The largest
# object in memory is the concatenation of 46 BUSCO full tables (~8,300 rows
# each), which is trivial; the synteny figures dominate, and those are a few
# hundred MB at most at 300 dpi.
#SBATCH --mem=32G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

module unload anaconda 2>/dev/null || true

# unset LD_PRELOAD: RCAC's XALT usage-tracking library is injected via
# LD_PRELOAD and fails on some nodes (GLIBC_2.33/2.34 mismatch), which can
# kill subshells under `set -e`. It's accounting only — safe to drop.
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"

# Must match REF_TAG in 05_ragtag_liftoff_array.sh. The conda env and
# the BUSCO-Plot-Py checkout below are deliberately NOT suffixed: they are tools,
# shared across reference trees, and rebuilding them per reference would be waste.
REF_TAG="ptarmigan"

FINAL_DIR="${PROJECT_DIR}/final_${REF_TAG}"
QC_DIR="${PROJECT_DIR}/qc_${REF_TAG}"
BUSCO_DIR="${QC_DIR}/busco"
PLOT_DIR="${QC_DIR}/plots"

CONDA_ENVS_DIR="${PROJECT_DIR}/conda_envs"
PLOT_ENV_DIR="${CONDA_ENVS_DIR}/buscoplot"
PLOT_PYTHON="${PLOT_ENV_DIR}/bin/python3"

BUSCOPLOTPY_DIR="${PROJECT_DIR}/tools/BUSCO-Plot-Py"
BUSCOPLOTPY_REPO="https://github.com/lorenzo-arcioni/BUSCO-Plot-Py.git"

PLOT_SCRIPT="${SLURM_SUBMIT_DIR}/09_busco_plots.py"

# ---- Plot options -----------------------------------------------------------
# Maximum chromosomes drawn per assembly. Only chr_* sequences are ever
# plotted: ~39 autosomes plus Z, W and MT, so 42 at most. The default is set
# above that deliberately. If an assembly exceeds this limit, the karyotype
# plot keeps only the chromosomes with the most BUSCO hits, which would drop
# real chromosomes from some assemblies and not others — and the comparison
# between plots is the whole point. It still caps the SYNTENY panels, where
# the largest chromosomes are kept for readability.
CHRS_LIMIT=50

# horizontal: the two assemblies stacked top and bottom (wide figure).
# vertical:   side by side (tall figure).
SYNTENY_ORIENTATION="horizontal"

# Which synteny comparisons to draw.
#   empty  -> one representative per species (lowest sample id, hap1),
#             compared pairwise: LEPC-GRPC, LEPC-STGR, GRPC-STGR.
#   custom -> comma-separated A:B pairs of assembly directory names, e.g.
#             "LEPC_F5540_hap1:GRPC_F5595_hap1,LEPC_F5540_hap1:STGR_F5457_hap1"
SYNTENY_PAIRS=""

# true also draws hap1 vs hap2 for every individual (23 extra figures). This is
# a useful phasing check — a well-phased pair should be a clean diagonal — but
# it roughly quadruples runtime, so it is off by default.
SYNTENY_HAPS=false

# Comma-separated list of plot types to skip on a rerun: barplot,karyotype,synteny
SKIP=""

# Write an .svg beside every .png. Vector output is what you want for anything
# going into a manuscript or the USFWS report; PNG stays for quick viewing.
SAVE_SVG=true

THREADS=$SLURM_CPUS_PER_TASK
export MPLBACKEND=Agg          # belt and braces; the .py also forces it
export OMP_NUM_THREADS="$THREADS"

mkdir -p logs "$PLOT_DIR" "$CONDA_ENVS_DIR" "$(dirname "$BUSCOPLOTPY_DIR")"

# =============================================================================
# PRE-FLIGHT
# =============================================================================
if [[ ! -f "$PLOT_SCRIPT" ]]; then
    echo "ERROR: companion script not found: ${PLOT_SCRIPT}"
    echo "       09_busco_plots.py must sit next to this file in the submission directory."
    exit 1
fi
if [[ ! -d "$BUSCO_DIR" ]]; then
    echo "ERROR: BUSCO output directory not found: ${BUSCO_DIR}"
    echo "       Run 07_busco_array.sh first."
    exit 1
fi

N_RUNS=$(find "$BUSCO_DIR" -maxdepth 2 -name 'short_summary.specific.*.json' | wc -l)
if [[ "$N_RUNS" -eq 0 ]]; then
    echo "ERROR: no completed BUSCO runs under ${BUSCO_DIR}"
    echo "       Expected <SPECIES>_<SAMPLE>_<HAP>/short_summary.specific.*.json"
    exit 1
fi
echo ">>> ${N_RUNS} completed BUSCO run(s) found"

# =============================================================================
# ONE-TIME SETUP — plotting env + BUSCO-Plot-Py checkout
# =============================================================================
if [[ ! -x "$PLOT_PYTHON" ]]; then
    echo ">>> Creating plotting environment: ${PLOT_ENV_DIR}"
    ml anaconda/2025.12-py313
    conda create --yes --override-channels --prefix "$PLOT_ENV_DIR" \
        -c conda-forge python=3.11 matplotlib seaborn pandas
    module unload anaconda 2>/dev/null || true
fi
if [[ ! -x "$PLOT_PYTHON" ]]; then
    echo "ERROR: plotting environment was not created at ${PLOT_ENV_DIR}"
    exit 1
fi

if [[ ! -d "${BUSCOPLOTPY_DIR}/buscoplotpy" ]]; then
    echo ">>> Cloning BUSCO-Plot-Py"
    rm -rf "$BUSCOPLOTPY_DIR"
    git clone --depth 1 "$BUSCOPLOTPY_REPO" "$BUSCOPLOTPY_DIR" || true
fi
if [[ ! -f "${BUSCOPLOTPY_DIR}/buscoplotpy/graphics/synteny.py" ]]; then
    echo "ERROR: BUSCO-Plot-Py checkout is missing or incomplete:"
    echo "       ${BUSCOPLOTPY_DIR}"
    echo "       Expected buscoplotpy/graphics/synteny.py (GitHub only — the"
    echo "       PyPI release does not ship the synteny module)."
    echo "       If the compute node cannot reach GitHub, clone it from a login"
    echo "       node and resubmit:"
    echo "         git clone --depth 1 ${BUSCOPLOTPY_REPO} ${BUSCOPLOTPY_DIR}"
    exit 1
fi

# The package is pure Python and its pyproject cannot be built (see INSTALL
# NOTE above), so it is used straight from the checkout.
export PYTHONPATH="${BUSCOPLOTPY_DIR}:${PYTHONPATH:-}"

echo ">>> python       : ${PLOT_PYTHON}"
echo ">>> buscoplotpy  : ${BUSCOPLOTPY_DIR}"
"$PLOT_PYTHON" - <<'PYCHECK'
import sys
try:
    import matplotlib, seaborn, pandas
    from buscoplotpy.graphics.synteny import horizontal_synteny_plot
    from buscoplotpy.graphics.karyoplot import karyoplot
    from buscoplotpy.graphics.organism_busco_barplot import organism_busco_barplot
except Exception as exc:                      # noqa: BLE001
    sys.exit(f"ERROR: plotting imports failed: {exc}")
print(f">>> matplotlib {matplotlib.__version__} | seaborn {seaborn.__version__} "
      f"| pandas {pandas.__version__} | buscoplotpy OK (synteny present)")
PYCHECK

# =============================================================================
# PLOTS
# =============================================================================
ARGS=(
    --busco-dir "$BUSCO_DIR"
    --final-dir "$FINAL_DIR"
    --out-dir   "$PLOT_DIR"
    --chrs-limit "$CHRS_LIMIT"
    --synteny-orientation "$SYNTENY_ORIENTATION"
)
# Plain `if` rather than `[[ ... ]] && ...`: the short-circuit form leaves the
# list's exit status at 1 when the test is false, which is harmless mid-script
# but silently marks a job FAILED if it ever ends up on the last line.
if [[ -n "$SYNTENY_PAIRS" ]]; then ARGS+=(--synteny-pairs "$SYNTENY_PAIRS"); fi
if [[ -n "$SKIP" ]];          then ARGS+=(--skip "$SKIP"); fi
if [[ "$SYNTENY_HAPS" == true ]]; then ARGS+=(--synteny-haps); fi
if [[ "$SAVE_SVG" != true ]]; then ARGS+=(--no-svg); fi

echo ">>> Started: $(date)"
"$PLOT_PYTHON" "$PLOT_SCRIPT" "${ARGS[@]}"

# =============================================================================
# DONE
# =============================================================================
echo ""
echo "============================================================"
echo ">>> Complete: $(date)"
echo "  completeness : ${PLOT_DIR}/completeness/"
echo "  karyotype    : ${PLOT_DIR}/karyotype/"
echo "  synteny      : ${PLOT_DIR}/synteny/"
find "$PLOT_DIR" -name '*.png' | wc -l | xargs printf "  %s PNG files written\n"
find "$PLOT_DIR" -name '*.svg' | wc -l | xargs printf "  %s SVG files written\n"
echo "============================================================"
