# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Student project (Centrale, "Modélisation et simulation numérique"): a Metropolis Monte Carlo code on a rigid bcc lattice for the order-disorder transition of CuZn (β-brass, B2 structure). The assignment and its questions (1 to 9: lattice and neighbour tables, Cowley parameters `cowley`, long-range order `ordre`, cohesive energy at T = 0, equilibrium runs, Tc, discussion) are in `descriptif_du_projet_msn.md`.

The students write the code themselves. When asked to check something, report the problems and verify them numerically (scripts in the scratchpad) rather than editing `src/model.py`, unless they explicitly ask for changes.

## Commands

Everything runs with the project venv `.venv`, set up from `requirements.txt` (setup steps for a new machine are in `README.md`). There is no test suite, linter or build step.

```bash
.venv/bin/python src/model.py                      # small demo run (__main__ block)
.venv/bin/jupyter nbconvert --to notebook --execute --inplace analysis.ipynb --ExecutePreprocessor.timeout=1800   # re-run the whole analysis (~3 min), outputs saved in the notebook
.venv/bin/python -B report/make_figures.py          # report figures, tables and report/figures/results.txt from cached report/figures/data.npz (--rerun: new simulation, ~6 min)
cd report && latexmk -pdf report.tex slides.tex     # build the report and the slides (pdflatex + biber)
```

To use the model from a script or notebook at the project root: `sys.path.insert(0, "src"); import model`. Run ad-hoc checks with `python -B` to avoid writing `src/__pycache__`.

## Architecture of `src/model.py`

Single module, procedural, with **global state**:

- `sous_reseau` (numpy int array, length `NB_SITES` = 432) holds the **species** on each site (`CU = 0`, `ZN = 1`), despite its name. The sublattice of site `i` is `i % 2`: even = α (cube corner, Cu in B2), odd = β (cube centre, Zn in B2). Every function reads or mutates this global; `init_at_t0()` resets it to perfect B2, and both `metropolis*` drivers call it first.
- Site index: `i = 2*(x + L*y + L**2*z) + s`, `L = 6`. `site_index()` applies `% L` on each coordinate, which is what gives periodic boundaries. Never build neighbours by adding offsets to `i` with `% NB_SITES`: that gives helical boundaries and, for first neighbours, a wrong lattice (an earlier bug).
- `FIRST_NEIGHBOURS` (8 per site, other sublattice; an α site looks at cells x-1/x, a β site at x/x+1) and `SECOND_NEIGHBOURS` (6 per site, same sublattice) are built once at import. They must stay **symmetric** (j in N(i) ⇔ i in N(j)): `get_full_energy()` sums local energies and divides by 2, and `energy_after_flip` relies on it.
- Energy: pair model with `J` (1st shell) and `I` (2nd shell) indexed `[species][species]`, in eV; temperatures in K, `kB` in eV/K.
- Moves are **exchanges** (Kawasaki) of a random site with a random first neighbour, even though functions say "flip". `energy_after_flip(E, i, j)` returns the energy after exchanging i and j **without modifying** `sous_reseau` (skips the i-j bond, returns `E` for same species); `flipping_proba` gives the Metropolis probability; `flipping_action` draws and performs the swap.
- Observables: `cowley()` returns `(alpha_1, alpha_2)`; its constant `4` assumes c_Zn = 1/2, valid because exchanges conserve composition. `order_parameter()` is η (the spec's `ordre`); the drivers store |η| because η can flip sign near Tc.
- Drivers (arguments in **steps**; 1 sweep = `NB_SITES` steps):
  - `metropolis_before_eq(n, T)`: energy after every step, starting from B2, used to find the equilibration length.
  - `metropolis(n_until_eq, n_after_eq, T, measure_every=1)`: returns numpy series `(E_total_eV, |eta|, alpha_1, alpha_2)` measured after equilibration. One `cowley()` + `order_parameter()` costs about 60 moves, so pass `measure_every=NB_SITES` for real runs.
- RNG is Python's global `random` module: seed it with `random.seed(...)` for reproducible runs.

Cost: about 7 µs per step, 3 ms per sweep; a temperature with ~1000 sweeps takes a few seconds.

## `analysis.ipynb`

Uses `model.py` to answer the project questions: T = 0 checks and cohesive energy, equilibration (energy from the B2 start: long runs at 8 temperatures, then `N_EQ` = 20000 steps at every scan temperature, with convergence value and equilibrium step marked; `N_EQ` is in steps, `N_MEAS` in sweeps), temperature scan (18 temperatures, `N_MEAS` = 800 sweeps) with E(T), heat capacity from energy fluctuations, ⟨|η|⟩, susceptibility χ = N(⟨η²⟩ − ⟨|η|⟩²)/k_BT, Cowley parameters vs ∓⟨η²⟩, |η| histograms and Tc estimates. Seeds are fixed, so results are reproducible. The heat capacity is the noisiest quantity; increase `N_MEAS` near Tc for final values.

## Reference values for checking changes

- Perfect B2: `get_full_energy()` = −1073.088 eV (−2.484 eV/atom), `cowley()` = (−1, 1), `order_parameter()` = 1.
- Random lattice: α₁, α₂ ≈ 0, ⟨|η|⟩ ≈ 0.04 (finite-size floor of the 432-site box).
- During any run, the incrementally tracked energy must equal `get_full_energy()` to ~1e-9 eV, and the Zn count must stay 216.
- Tc from the notebook: about 700 K (heat capacity peak), 725 K (steepest drop of |η|), 725–750 K (susceptibility peak); each moves by ±25 K between runs; experimental β-brass ≈ 740 K.
