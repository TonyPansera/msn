# CuZn order-disorder transition by Monte Carlo

Metropolis Monte Carlo simulation on a rigid lattice of the order-disorder transition of CuZn (β-brass), project of the course *Modélisation et simulation numérique* (Centrale). The full assignment is in [`descriptif_du_projet_msn.md`](descriptif_du_projet_msn.md).

## Contents

| File | Role |
| :--- | :--- |
| `src/model.py` | The model: bcc lattice and neighbour tables, pair energy, Metropolis exchanges, Cowley parameters, long-range order parameter, simulation drivers |
| `analysis.ipynb` | Analysis: checks at T = 0 K, equilibration time, temperature scan, plots, estimate of the transition temperature |
| `descriptif_du_projet_msn.md` | Project description and questions |
| `requirements.txt` | Python packages to install |

## Installation on a new machine

You need **Python 3.10 or newer** and **git**.

1. Get the project:

   ```bash
   git clone <repository-url>
   cd <project-folder>
   ```

2. Create a virtual environment in the project folder and install the packages.

   Linux / macOS / WSL:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt
   ```

   Windows (PowerShell):

   ```powershell
   py -m venv .venv
   .venv\Scripts\Activate.ps1
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt
   ```

   If PowerShell refuses to run `Activate.ps1`, run `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once, or use `.venv\Scripts\activate.bat` from `cmd`.

3. Check that the model works (it runs a short simulation and prints the lattice and the energies):

   ```bash
   python src/model.py
   ```

The virtual environment has to be activated (step 2, `activate` line) in every new terminal before running the commands below.

## Usage

### Analysis notebook

Open `analysis.ipynb` and run all the cells. A full run takes about 3 minutes.

- **VS Code**: install the *Python* and *Jupyter* extensions, open the project folder, open `analysis.ipynb`, click *Select Kernel* at the top right and choose the `.venv` environment, then *Run All*.
- **JupyterLab**: from the project folder, with the environment activated, run `jupyter lab analysis.ipynb`.
- **Without opening it** (runs every cell and saves the results and plots in the notebook):

  ```bash
  jupyter nbconvert --to notebook --execute --inplace analysis.ipynb --ExecutePreprocessor.timeout=1800
  ```

The notebook must be opened from the project folder: it imports the model with `sys.path.insert(0, "src")`.

### Using the model in your own script

Run from the project folder:

```python
import sys
sys.path.insert(0, "src")
import model

SWEEP = model.NB_SITES  # 1 sweep = 432 exchange attempts

# 150 sweeps to reach equilibrium, then 800 measured sweeps at 700 K, one measurement per sweep
energies, eta, alpha_1, alpha_2 = model.metropolis(150 * SWEEP, 800 * SWEEP, 700, measure_every=SWEEP)
print(energies.mean() / model.NB_SITES, eta.mean(), alpha_1.mean(), alpha_2.mean())
```

The main functions of `src/model.py`:

| Function | Returns |
| :--- | :--- |
| `metropolis(n_until_eq, n_after_eq, T, measure_every=1)` | Series measured after equilibrium: total energy (eV), \|η\|, α₁, α₂ |
| `metropolis_before_eq(n, T)` | Energy after every step, starting from the ordered B2 structure |
| `get_full_energy()` | Total energy of the current configuration (eV) |
| `cowley()` | Cowley short-range order parameters (α₁, α₂) of the current configuration |
| `order_parameter()` | Long-range order parameter η of the current configuration |

The numbers of steps are given in exchange attempts; a sweep is `NB_SITES` = 432 attempts. Both drivers start from the ordered B2 structure. Results are random: call `random.seed(...)` before a run to make it reproducible.

## Model

- **Lattice**: bcc, 6 × 6 × 6 cubic cells, 432 sites, periodic boundary conditions. Site `i = 2*(x + 6*y + 36*z) + s`: even sites (s = 0) are the cube corners (sublattice α), odd sites (s = 1) the cube centres (sublattice β).
- **Energy**: sum over 1st-neighbour pairs (8 per site, distance √3·a/2) of J and over 2nd-neighbour pairs (6 per site, distance a) of I, with the values of the assignment (in eV).
- **Moves**: exchange of an atom with a random first neighbour, accepted with the Metropolis probability min(1, exp(−ΔE / k_B T)). The composition (216 Cu, 216 Zn) is conserved.
- **Observables**: long-range order η = (N_Cu(α) − N_Cu(β)) / N_Cu, Cowley parameters α_n = 1 − p_n(Zn|Cu) / c_Zn for the first two shells.

Reference values (useful to check a modification): at T = 0 K (perfect B2), E = −1073.088 eV = −2.484 eV/atom, (α₁, α₂) = (−1, +1), η = 1.
