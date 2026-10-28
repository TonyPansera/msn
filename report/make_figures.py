"""Simulations, figures and tables of the report and the slides.

Run from the project root:
    .venv/bin/python -B report/make_figures.py           # simulate once (~6 min), then plot
    .venv/bin/python -B report/make_figures.py --rerun   # force a new simulation

The raw data are cached in report/figures/data.npz: without --rerun only the
checks (a few seconds) and the plots are redone. Everything is written in
report/figures/.
"""
import math
import platform
import random
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import model

OUT = ROOT / "report" / "figures"
OUT.mkdir(exist_ok=True)

N = model.NB_SITES
SWEEP = model.NB_SITES  # 1 sweep = NB_SITES exchange attempts

# Simulation parameters (temperatures in K, N_EQ in steps, N_MEAS in sweeps)
T_SCAN = [200, 300, 400, 500, 550, 600, 650, 675, 700, 710, 720, 730, 740, 750, 760, 775, 800, 850, 900, 1000, 1200, 1600]
N_EQ = 50000          # equilibration steps before measuring, at every temperature
N_MEAS = 4000         # measured sweeps, one measurement per sweep
N_RELAX = 200000      # steps of the equilibration test runs
T_RELAX_SHOWN = [500, 700, 730, 1000]  # relaxation curves drawn in the figure
SEED = 2026

# Colours (categorical slots in fixed order, ink and chrome for a white page)
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
SERIES = [BLUE, ORANGE, AQUA, YELLOW]
CU_COLOR, ZN_COLOR = ORANGE, BLUE
INK, INK_2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"

CM = 1 / 2.54
plt.rcParams.update({
    "font.size": 9.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9.5,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "mathtext.fontset": "dejavusans",
    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "lines.linewidth": 1.5,
    "lines.markersize": 5,
    "legend.frameon": False,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
})


# ---------------------------------------------------------------- checks

def run_checks() -> dict:
    """Tests of the program (question 4) and cohesive energy (question 5)."""
    random.seed(SEED)
    np.random.seed(SEED)
    c = {}

    # Neighbour tables: 8 and 6 distinct neighbours, symmetric (j in N(i) <=> i in N(j))
    first_ok = all(len(set(nb)) == 8 for nb in model.FIRST_NEIGHBOURS)
    second_ok = all(len(set(nb)) == 6 for nb in model.SECOND_NEIGHBOURS)
    sym_1 = all(i in model.FIRST_NEIGHBOURS[j] for i in range(N) for j in model.FIRST_NEIGHBOURS[i])
    sym_2 = all(i in model.SECOND_NEIGHBOURS[j] for i in range(N) for j in model.SECOND_NEIGHBOURS[i])
    other_sublattice = all((i + j) % 2 == 1 for i in range(N) for j in model.FIRST_NEIGHBOURS[i])
    same_sublattice = all((i + j) % 2 == 0 for i in range(N) for j in model.SECOND_NEIGHBOURS[i])
    c["neighbour_tables_ok"] = first_ok and second_ok and sym_1 and sym_2 and other_sublattice and same_sublattice

    # Perfect B2 at T = 0 K
    model.init_at_t0()
    c["E_B2_total"] = model.get_full_energy()
    c["E_B2"] = c["E_B2_total"] / N
    c["alpha_1_B2"], c["alpha_2_B2"] = model.cowley()
    c["eta_B2"] = model.order_parameter()

    # Pure elements in the same bcc pair model: 4 first-neighbour bonds per atom
    c["E_pure_Cu"] = 4 * model.J[model.CU][model.CU]
    c["E_pure_Zn"] = 4 * model.J[model.ZN][model.ZN]
    c["E_formation"] = c["E_B2"] - (c["E_pure_Cu"] + c["E_pure_Zn"]) / 2

    # Ideal random alloy (c = 1/2): pairs Cu-Cu, Zn-Zn, Cu-Zn with probabilities 1/4, 1/4, 1/2
    J_mean = (model.J[0][0] + model.J[1][1] + 2 * model.J[0][1]) / 4
    I_mean = (model.I[0][0] + model.I[1][1] + 2 * model.I[0][1]) / 4
    c["E_random"] = 4 * J_mean + 3 * I_mean
    c["dE_order"] = c["E_random"] - c["E_B2"]
    # Bragg-Williams mean field: E(eta) = E_random - dE eta^2 and S(eta) ~ kB (ln 2 - eta^2/2), so kB Tc = 2 dE
    c["Tc_BW"] = 2 * c["dE_order"] / model.kB

    # Random lattices: no order expected, |eta| only keeps a finite-size floor
    random_eta, random_a1, random_a2 = [], [], []
    for k in range(2000):
        np.random.shuffle(model.sous_reseau)
        random_eta.append(abs(model.order_parameter()))
        if k < 200:
            a_1, a_2 = model.cowley()
            random_a1.append(a_1)
            random_a2.append(a_2)
    c["eta_random"] = float(np.mean(random_eta))
    c["alpha_1_random"] = float(np.mean(random_a1))
    c["alpha_2_random"] = float(np.mean(random_a2))
    # Number of Cu on alpha is hypergeometric (216 alpha sites drawn among 432, 216 Cu)
    half = N // 2
    sigma_n = math.sqrt(half * half * (N - half) * (N - half) / (N**2 * (N - 1)))
    c["eta_random_theory"] = math.sqrt(2 / math.pi) * 2 * sigma_n / half

    # Energy tracked step by step vs full recomputation, and composition
    energies = model.metropolis_before_eq(100000, 730)
    c["energy_drift"] = abs(energies[-1] - model.get_full_energy())
    c["n_zn_after_run"] = int(np.count_nonzero(model.sous_reseau == model.ZN))
    return c


# ---------------------------------------------------------------- simulations

def run_simulations() -> dict:
    """Equilibration test (question 6) and temperature scan (question 7)."""
    t_start = time.perf_counter()

    # Equilibration: energy from the B2 start, at every temperature of the scan
    random.seed(SEED + 1)
    n_eq = []
    relax_curves = []
    for T in T_SCAN:
        energies = np.array(model.metropolis_before_eq(N_RELAX, T)) / N
        E_eq = energies[N_RELAX // 2:].mean()  # convergence value: mean over the second half
        sigma = energies[N_RELAX // 2:].std()  # size of the equilibrium fluctuations
        n_eq.append(int(np.argmax(energies >= E_eq - sigma)))  # first step inside the fluctuations
        if T in T_RELAX_SHOWN:
            relax_curves.append(energies[:100000:50])
        print(f"relaxation T = {T:5d} K: E_eq = {E_eq:.4f} eV/atom, n_eq = {n_eq[-1]} steps", flush=True)

    # Temperature scan: every temperature starts again from B2
    random.seed(SEED + 2)
    E_series, eta_series, a1_series, a2_series, snapshots = [], [], [], [], []
    for T in T_SCAN:
        energies, eta, alpha_1, alpha_2 = model.metropolis(N_EQ, N_MEAS * SWEEP, T, measure_every=SWEEP)
        E_series.append(energies)
        eta_series.append(eta)
        a1_series.append(alpha_1)
        a2_series.append(alpha_2)
        snapshots.append(model.sous_reseau.copy())  # last configuration of the run
        print(f"scan T = {T:5d} K: E = {energies.mean() / N:.4f} eV/atom, |eta| = {eta.mean():.3f}", flush=True)

    return {
        "T": np.array(T_SCAN, dtype=float),
        "n_eq": np.array(n_eq),
        "relax_T": np.array(T_RELAX_SHOWN),
        "relax_curves": np.array(relax_curves),
        "E": np.array(E_series),
        "eta": np.array(eta_series),
        "alpha_1": np.array(a1_series),
        "alpha_2": np.array(a2_series),
        "snapshots": np.array(snapshots),
        "run_time": time.perf_counter() - t_start,
    }


# ---------------------------------------------------------------- analysis

def analyse(d: dict) -> dict:
    T = d["T"]
    r = {"T": T}
    r["E"] = d["E"].mean(axis=1) / N
    r["C"] = d["E"].var(axis=1) / (N * (model.kB * T) ** 2)  # C / (N kB)
    r["eta"] = d["eta"].mean(axis=1)
    r["eta2"] = (d["eta"] ** 2).mean(axis=1)
    r["chi"] = N * d["eta"].var(axis=1) / (model.kB * T)  # eV^-1
    r["a1"] = d["alpha_1"].mean(axis=1)
    r["a2"] = d["alpha_2"].mean(axis=1)
    r["Tc_C"] = T[np.argmax(r["C"])]
    r["Tc_chi"] = T[np.argmax(r["chi"])]
    r["Tc_eta"] = T[np.argmin(np.gradient(r["eta"], T))]
    return r


# ---------------------------------------------------------------- figures

def style_axis(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(axis="x", visible=False)
    ax.tick_params(length=3, color=AXIS)


def mark_tc(ax, Tc, label):
    ax.axvline(Tc, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.annotate(label, (Tc, 1), xycoords=("data", "axes fraction"), xytext=(4, -10),
                textcoords="offset points", fontsize=8, color=INK_2)


def plot_series(ax, T, y, color, marker="o", label=None):
    ax.plot(T, y, marker + "-", color=color, ms=4.5, mec="white", mew=0.8, label=label)


def fig_equilibration(d):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15 * CM, 5.6 * CM), width_ratios=[1.15, 1])
    sweeps = np.arange(d["relax_curves"].shape[1]) * 50 / SWEEP
    for k, (T, curve) in enumerate(zip(d["relax_T"], d["relax_curves"])):
        ax1.plot(sweeps, curve, color=SERIES[k], lw=1, label=f"{T} K")
    ax1.axhline(model.J[0][1] * 4, color=MUTED, lw=0.8, ls=":")
    ax1.annotate("B2 (T = 0 K)", (120, model.J[0][1] * 4), xytext=(0, -9), textcoords="offset points",
                 ha="right", fontsize=8, color=INK_2)
    ax1.set_xlim(0, 120)
    ax1.set_ylim(bottom=model.J[0][1] * 4 - 0.004)
    ax1.legend(loc="center left", bbox_to_anchor=(1, 0.5), handlelength=1.2)
    style_axis(ax1, "sweeps from the B2 start", "E / N (eV/atom)")
    ax1.set_title("Energy from the B2 start", loc="left")

    ax2.plot(d["T"], d["n_eq"] / 1000, "o-", color=BLUE, ms=4.5, mec="white", mew=0.8)
    ax2.axhline(N_EQ / 1000, color=INK_2, lw=1, ls="--")
    ax2.annotate(f"$N_{{eq}}$ = {N_EQ} steps used", (d["T"][0], N_EQ / 1000), xytext=(0, 3),
                 textcoords="offset points", ha="left", fontsize=8, color=INK_2)
    ax2.set_ylim(0, 1.25 * max(N_EQ, d["n_eq"].max()) / 1000)
    style_axis(ax2, "T (K)", r"$n_{eq}$ (10$^3$ steps)")
    ax2.set_title("Steps needed to reach equilibrium", loc="left")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(OUT / "equilibration.pdf")
    plt.close(fig)


def fig_snapshots(d, r):
    T_shown = [300, int(r["Tc_chi"]), 1200]
    fig, axes = plt.subplots(1, 3, figsize=(15 * CM, 5.6 * CM))
    L = model.L
    for ax, T in zip(axes, T_shown):
        config = d["snapshots"][T_SCAN.index(T)]
        xs, ys, colors = [], [], []
        for y in range(L):
            for x in range(L):
                for s in (0, 1):  # alpha site at (x, y, 0), beta site at (x+1/2, y+1/2, 1/2)
                    xs.append(x + 0.5 * s)
                    ys.append(y + 0.5 * s)
                    colors.append(CU_COLOR if config[model.site_index(x, y, 0, s)] == model.CU else ZN_COLOR)
        ax.scatter(xs, ys, c=colors, s=38, edgecolors="white", linewidths=0.8, zorder=3)
        ax.set_xticks(range(L))
        ax.set_yticks(range(L))
        ax.set_xlim(-0.5, L - 0.25)
        ax.set_ylim(-0.5, L - 0.25)
        ax.set_aspect("equal")
        ax.tick_params(labelbottom=False, labelleft=False, length=0)
        for side in ("left", "bottom"):
            ax.spines[side].set_visible(False)
        n_cu_alpha = np.count_nonzero(config[0::2] == model.CU)
        n_cu_beta = np.count_nonzero(config[1::2] == model.CU)
        eta = abs(n_cu_alpha - n_cu_beta) / (n_cu_alpha + n_cu_beta)
        ax.set_title(f"T = {T} K   (|η| = {eta:.2f})")
    handles = [plt.Line2D([], [], ls="", marker="o", ms=7, mfc=col, mec="white", label=lab)
               for col, lab in ((CU_COLOR, "Cu"), (ZN_COLOR, "Zn"))]
    fig.legend(handles=handles, loc="lower center", ncols=2, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(OUT / "snapshots.pdf")
    plt.close(fig)


def fig_energy_heat(c, r):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15 * CM, 5.6 * CM))
    plot_series(ax1, r["T"], r["E"], BLUE)
    for value, label, x, ha in ((c["E_B2"], "perfect B2", r["T"][-1], "right"),
                                (c["E_random"], "random alloy", r["T"][0], "left")):
        ax1.axhline(value, color=MUTED, lw=0.8, ls=":")
        ax1.annotate(label, (x, value), xytext=(0, 3), textcoords="offset points",
                     ha=ha, fontsize=8, color=INK_2)
    style_axis(ax1, "T (K)", r"$\langle E \rangle / N$ (eV/atom)")
    ax1.set_title("Mean energy", loc="left")

    plot_series(ax2, r["T"], r["C"], BLUE)
    mark_tc(ax2, r["Tc_C"], f"{r['Tc_C']:.0f} K")
    style_axis(ax2, "T (K)", r"$C / (N k_B)$")
    ax2.set_title("Heat capacity (energy fluctuations)", loc="left")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(OUT / "energy_heat.pdf")
    plt.close(fig)


def fig_order_chi(c, r):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15 * CM, 5.6 * CM))
    plot_series(ax1, r["T"], r["eta"], BLUE)
    ax1.axhline(c["eta_random"], color=MUTED, lw=0.8, ls=":")
    ax1.annotate("random lattice", (r["T"][0], c["eta_random"]), xytext=(0, 3),
                 textcoords="offset points", ha="left", fontsize=8, color=INK_2)
    mark_tc(ax1, r["Tc_eta"], f"{r['Tc_eta']:.0f} K")
    ax1.set_ylim(0, 1.05)
    style_axis(ax1, "T (K)", r"$\langle |\eta| \rangle$")
    ax1.set_title("Long-range order parameter", loc="left")

    plot_series(ax2, r["T"], r["chi"], BLUE)
    mark_tc(ax2, r["Tc_chi"], f"{r['Tc_chi']:.0f} K")
    style_axis(ax2, "T (K)", r"$\chi$ (eV$^{-1}$)")
    ax2.set_title(r"Susceptibility ($|\eta|$ fluctuations)", loc="left")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(OUT / "order_chi.pdf")
    plt.close(fig)


def fig_cowley(r):
    fig, ax = plt.subplots(figsize=(15 * CM, 6 * CM))
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.plot(r["T"], -r["eta2"], "--", color=BLUE, lw=1, alpha=0.7, label=r"$-\langle \eta^2 \rangle$")
    ax.plot(r["T"], r["eta2"], "--", color=ORANGE, lw=1, alpha=0.7, label=r"$+\langle \eta^2 \rangle$")
    plot_series(ax, r["T"], r["a1"], BLUE, "o", r"$\alpha_1$")
    plot_series(ax, r["T"], r["a2"], ORANGE, "s", r"$\alpha_2$")
    ax.set_ylim(-1.08, 1.08)
    handles, labels = ax.get_legend_handles_labels()
    order = [2, 3, 0, 1]
    ax.legend([handles[k] for k in order], [labels[k] for k in order], loc="center left", ncols=2,
              bbox_to_anchor=(0.01, 0.5), columnspacing=1.2)
    style_axis(ax, "T (K)", "Cowley parameter")
    fig.tight_layout()
    fig.savefig(OUT / "cowley.pdf")
    plt.close(fig)


def fig_eta_hist(d, r):
    T_shown = [650, int(r["Tc_chi"]), 900]
    # |eta| only takes the values 2k/216: each bin holds 3 of them, so that no bin is favoured
    step = 2 / (N // 2)
    bins = (np.arange(0, N // 4 + 3, 3) - 0.5) * step
    fig, axes = plt.subplots(1, 3, figsize=(15 * CM, 4.8 * CM), sharey=True)
    for ax, T in zip(axes, T_shown):
        ax.hist(d["eta"][T_SCAN.index(T)], bins=bins, color=BLUE, rwidth=0.85)
        ax.set_title(f"T = {T} K", loc="left")
        style_axis(ax, r"$|\eta|$", "")
        ax.grid(axis="y", visible=True)
    axes[0].set_ylabel("number of sweeps")
    fig.tight_layout(w_pad=1)
    fig.savefig(OUT / "eta_hist.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- tables and summary

def number(value: float, decimals: int) -> str:
    return f"{value:.{decimals}f}".replace("-", "$-$")


def write_table(r):
    rows = []
    for k, T in enumerate(r["T"]):
        rows.append(" & ".join([
            f"{T:.0f}",
            number(r["E"][k], 4),
            number(r["C"][k], 2),
            number(r["eta"][k], 3),
            number(r["chi"][k], 0),
            number(r["a1"][k], 3),
            number(r["a2"][k], 3),
        ]) + r" \\")
    # Whole tabular in the file: \input between two table rows breaks the alignment
    header = [r"\begin{tabular}{rllllll}", r"\toprule",
              r"$T$ (K) & $\langle E\rangle/N$ & $C/(Nk_B)$ & $\langle|\eta|\rangle$ & $\chi$ "
              r"& $\langle\alpha_1\rangle$ & $\langle\alpha_2\rangle$ \\", r"\midrule"]
    footer = [r"\bottomrule", r"\end{tabular}"]
    (OUT / "scan_table.tex").write_text("\n".join(header + rows + footer) + "\n")


def write_summary(c, d, r):
    lines = [
        f"Python {platform.python_version()}, NumPy {np.__version__}, matplotlib {matplotlib.__version__}",
        f"Simulation time (equilibration test + scan): {float(d['run_time']):.0f} s",
        f"Parameters: L = {model.L}, N = {N}, N_EQ = {N_EQ} steps, N_MEAS = {N_MEAS} sweeps, "
        f"N_RELAX = {N_RELAX} steps, seed = {SEED}",
        "",
        "Checks",
    ]
    lines += [f"  {key} = {value}" for key, value in c.items()]
    lines += ["", "Equilibration (steps needed from the B2 start)"]
    lines += [f"  T = {T:5.0f} K: n_eq = {n}" for T, n in zip(d["T"], d["n_eq"])]
    lines += [f"  max n_eq = {d['n_eq'].max()} steps at T = {d['T'][np.argmax(d['n_eq'])]:.0f} K", ""]
    lines += ["Transition temperature",
              f"  Tc from the heat capacity peak:     {r['Tc_C']:.0f} K",
              f"  Tc from the susceptibility peak:    {r['Tc_chi']:.0f} K",
              f"  Tc from the steepest drop of |eta|: {r['Tc_eta']:.0f} K",
              f"  Tc Bragg-Williams mean field:       {c['Tc_BW']:.0f} K", ""]
    lines.append(f"{'T':>6} {'E':>8} {'C':>6} {'|eta|':>6} {'chi':>6} {'a1':>7} {'a2':>7}")
    for k, T in enumerate(r["T"]):
        lines.append(f"{T:6.0f} {r['E'][k]:8.4f} {r['C'][k]:6.3f} {r['eta'][k]:6.3f} {r['chi'][k]:6.1f} "
                     f"{r['a1'][k]:+7.3f} {r['a2'][k]:+7.3f}")
    text = "\n".join(lines) + "\n"
    (OUT / "results.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    checks = run_checks()
    data_file = OUT / "data.npz"
    if "--rerun" in sys.argv or not data_file.exists():
        data = run_simulations()
        np.savez_compressed(data_file, **data)
    data = dict(np.load(data_file))
    results = analyse(data)
    fig_equilibration(data)
    fig_snapshots(data, results)
    fig_energy_heat(checks, results)
    fig_order_chi(checks, results)
    fig_cowley(results)
    fig_eta_hist(data, results)
    write_table(results)
    write_summary(checks, data, results)
