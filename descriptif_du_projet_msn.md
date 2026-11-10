# Chapitre 2
# Study of the order-disorder transition of the CuZn compound by the Monte Carlo method on a rigid lattice (IK) **

## Objective
You will here write yourselves a Metropolis Monte Carlo code on a rigid lattice, based on a simplified pair model (a two-species generalized Ising-type model), restricted to first and second neighbors.

## 2.1 Physical model
The $\text{CuZn}$ compound ($\beta$-brass) crystallizes at low temperature in the ordered $\text{B2}$ structure ($\text{CsCl}$ type), made of two interpenetrating simple cubic sublattices.

The elementary cubic cell of parameter $a$ contains $2$ sites:
* $\text{Cu}$: $(0, 0, 0)$
* $\text{Zn}$: $(\frac{1}{2}, \frac{1}{2}, \frac{1}{2})$

The underlying lattice (site positions, independent of the species occupying them) is a body-centered cubic (bcc) lattice. The order-disorder transition consists in randomly permuting $\text{Cu}$ and $\text{Zn}$ atoms on this fixed bcc lattice (frozen positions — hence "rigid lattice"); only the occupation ($\text{Cu}$ or $\text{Zn}$) of each site varies.

Build a cubic box of $6 \times 6 \times 6$ cells (in units of $a$), i.e. $6^3 \times 2 = 432$ sites in total, with periodic boundary conditions (PBC) in the three directions.

At $T = 0\text{ K}$, the reference configuration is the ordered $\text{B2}$ structure described above (one $\text{Cu}$ site and one $\text{Zn}$ site per elementary cell, distributed exactly as in the cell repeated $6 \times 6 \times 6$ times).

For a given bcc site, we consider two neighbor shells:

| Shell | Distance | Number of neighbors $z$ |
| :--- | :---: | :---: |
| $1^{\text{st}}$ neighbors | $\frac{\sqrt{3}a}{2}$ | $8$ |
| $2^{\text{nd}}$ neighbors | $a$ | $6$ |

**Practical advice:** the simplest approach is to generate, once and for all at the beginning of the program, the list of the $8$ indices of the $1^{\text{st}}$ neighbors and of the $6$ indices of the $2^{\text{nd}}$ neighbors for each site $i$, applying PBC (modulo $6$ on each cell coordinate). This neighbor table never changes during the simulation — only the occupation of the sites changes.

The total energy of the system is a sum over pairs of first- and second-neighbor sites, each pair contributing according to the nature of the two atoms occupying it:

$$E = \sum_{\langle i,j \rangle_1} J_{\sigma_i \sigma_j} + \sum_{\langle i,j \rangle_2} I_{\sigma_i \sigma_j} \tag{2.1}$$

where $\sigma_i \in \{\text{Cu}, \text{Zn}\}$ is the species occupying site $i$, the first sum runs over pairs of $1^{\text{st}}$ neighbors, the second over pairs of $2^{\text{nd}}$ neighbors, and

$$J_{\sigma_i \sigma_j} \in \{J_{\text{CuCu}}, J_{\text{ZnZn}}, J_{\text{CuZn}}\}, \qquad I_{\sigma_i \sigma_j} \in \{I_{\text{CuCu}}, I_{\text{ZnZn}}, I_{\text{CuZn}}\} \tag{2.2}$$

according to the nature of the two atoms of the pair.

| Parameter | Value ($\text{eV}$) |
| :--- | :---: |
| $J_{\text{CuCu}}$ | $-0.873$ |
| $J_{\text{ZnZn}}$ | $-0.338$ |
| $J_{\text{CuZn}}$ | $-0.621$ |
| $I_{\text{CuCu}}$ | $0.000$ |
| $I_{\text{ZnZn}}$ | $0.000$ |
| $I_{\text{CuZn}}$ | $+0.005$ |

## 2.2 Simulation of disorder

1. Write the Monte Carlo program on a rigid lattice adapted to the study of the order-disorder transition of $\text{CuZn}$: construction of the crystal lattice and of the neighbor table ($1^{\text{st}}$ and $2^{\text{nd}}$ neighbors with PBC), potential energy of the system from the pair model above, exchange of atoms, Metropolis algorithm,...
2. Write the subprogram computing the Cowley short-range order parameters at first and second neighbors (function `cowley`).
   *Reminder:* the Cowley parameter for shell $n$ is written, for a reference site of species $\text{Cu}$, as
   $$\alpha_n = 1 - \frac{p_n^{\text{Zn}|\text{Cu}}}{c_{\text{Zn}}} \tag{2.3}$$
   where $p_n^{\text{Zn}|\text{Cu}}$ is the probability that a site of shell $n$ around a $\text{Cu}$ site is occupied by $\text{Zn}$, and $c_{\text{Zn}}$ is the overall $\text{Zn}$ concentration ($c_{\text{Zn}} = 1/2$ in $\text{CuZn}$).
3. Write the subprogram computing the long-range order parameter (function `ordre`).
   For the $\text{B2}$ structure of $\text{CuZn}$, the long-range order parameter is defined by
   $$\eta = \frac{N_{\text{Cu}}^{(\alpha)} - N_{\text{Cu}}^{(\beta)}}{N_{\text{Cu}}^{(\alpha)} + N_{\text{Cu}}^{(\beta)}} \tag{2.4}$$
   where $N_{\text{Cu}}^{(\alpha)}$ and $N_{\text{Cu}}^{(\beta)}$ denote, respectively, the number of $\text{Cu}$ atoms present on the $\alpha$ and $\beta$ sublattices.
4. Test the complete program.
5. Determine the associated cohesive energy (at $T = 0\text{ K}$), starting from the ordered $\text{B2}$ structure.
6. Generate equilibrium configurations at different temperatures.
7. Calculate the different physical quantities as a function of temperature.
8. Determine the order-disorder transition temperature.
9. Discuss the results.

---

# Bibliographie
[1] W.L. Bragg et E.J. Williams, *Proc. Roy. Soc. A* **145**, 699 (1934).

[2] C. Kittel, *Introduction to Solid State Physics*, Wiley, chapitres consacrés aux alliages ordonnés.

[3] S.C. Moss, *J. Appl. Phys.* **35**, 3547 (1964).

[4] T.B. Massalski (Ed.), *Binary Alloy Phase Diagrams*, ASM International.

[5] K. Binder et D.W. Heermann, *Monte Carlo Simulation in Statistical Physics*, Springer.