import math

import random
import numpy as np

# Constants of the system for CuZn

L = 6
NB_SITES = 2 * (L ** 3)
kB = 8.617333262e-5

CU = 0
ZN = 1

J = [[ -0.873, -0.621 ],
     [ -0.621, -0.338]]

I = [[ 0, 0.005 ], 
     [ 0.005, 0 ]]

sous_reseau = np.zeros(NB_SITES, dtype=int)

# Init at T=0, even are CU and odd are ZN
def init_at_t0():
    sous_reseau[0::2] = CU
    sous_reseau[1::2] = ZN

# Site index: i = 2*(x + L*y + L**2*z) + s
# (x, y, z) = cell coordinates in 0..L-1, s = 0 for the corner site (alpha), s = 1 for the centre site (beta)
def site_index(x: int, y: int, z: int, s: int) -> int:
    return 2 * ((x % L) + L * (y % L) + L**2 * (z % L)) + s

def site_coords(i: int) -> tuple[int, int, int, int]:
    c = i // 2
    return c % L, (c // L) % L, c // L**2, i % 2

def get_first_neighbours(i: int) -> list[int]:
    x, y, z, s = site_coords(i)
    if s == 0:
        # alpha site (corner): the 8 beta neighbours are the centres of the 8 cells sharing this corner
        return [site_index(x,   y,   z,   1), site_index(x-1, y,   z,   1),
                site_index(x,   y-1, z,   1), site_index(x-1, y-1, z,   1),
                site_index(x,   y,   z-1, 1), site_index(x-1, y,   z-1, 1),
                site_index(x,   y-1, z-1, 1), site_index(x-1, y-1, z-1, 1)]
    else:
        # beta site (centre): the 8 alpha neighbours are the corners of its own cell
        return [site_index(x,   y,   z,   0), site_index(x+1, y,   z,   0),
                site_index(x,   y+1, z,   0), site_index(x+1, y+1, z,   0),
                site_index(x,   y,   z+1, 0), site_index(x+1, y,   z+1, 0),
                site_index(x,   y+1, z+1, 0), site_index(x+1, y+1, z+1, 0)]

def get_second_neighbours(i: int) -> list[int]:
    x, y, z, s = site_coords(i)
    return [site_index(x+1, y, z, s), site_index(x-1, y, z, s),
            site_index(x, y+1, z, s), site_index(x, y-1, z, s),
            site_index(x, y, z+1, s), site_index(x, y, z-1, s)]


def choose_random_atom() -> int:
    return random.randint(0, NB_SITES-1)


FIRST_NEIGHBOURS = [get_first_neighbours(i) for i in range(NB_SITES)]
SECOND_NEIGHBOURS = [get_second_neighbours(i) for i in range(NB_SITES)]


def choose_random_neighbour(i: int) -> int:
    return FIRST_NEIGHBOURS[i][random.randint(0, len(FIRST_NEIGHBOURS[i]) - 1)]

def get_energy_neighbours(i: int) -> float:
    E = 0   
    for j in FIRST_NEIGHBOURS[i]:
        E += J[sous_reseau[i]][sous_reseau[j]]
    for j in SECOND_NEIGHBOURS[i]:
        E += I[sous_reseau[i]][sous_reseau[j]]
    return E



def get_full_energy() -> float:
    E = 0
    for i in range(0, NB_SITES):
        E += get_energy_neighbours(i)
    return E/2

# Energy after exchanging the atoms at sites i and j, computed without touching sous_reseau
def energy_after_flip(E: float, i: int, j: int) -> float:
    a = sous_reseau[i]
    b = sous_reseau[j]
    if a == b:
        return E  # same species: the exchange changes nothing
    delta = 0
    # bonds of i: species a becomes b. The i-j bond is skipped, it stays a Cu-Zn bond
    for k in FIRST_NEIGHBOURS[i]:
        if k != j:
            delta += J[b][sous_reseau[k]] - J[a][sous_reseau[k]]
    for k in SECOND_NEIGHBOURS[i]:
        if k != j:
            delta += I[b][sous_reseau[k]] - I[a][sous_reseau[k]]
    # bonds of j: species b becomes a
    for k in FIRST_NEIGHBOURS[j]:
        if k != i:
            delta += J[a][sous_reseau[k]] - J[b][sous_reseau[k]]
    for k in SECOND_NEIGHBOURS[j]:
        if k != i:
            delta += I[a][sous_reseau[k]] - I[b][sous_reseau[k]]
    return E + delta

def flipping_proba(E_i: float, E_j: float, T: float) -> float:
    if E_j > E_i:
        return math.exp(-(1/(kB*T)) * (E_j - E_i))
    else:
        return 1

def flipping_action(P: float, i: int, j: int) -> bool:
    p = random.uniform(0, 1)
    if P == 1 or p <= P:
        sous_reseau[i], sous_reseau[j] = sous_reseau[j], sous_reseau[i]
        return True
    return False

def cowley() -> tuple[float, float]:
    nb_zn_1 = 0
    nb_zn_2 = 0
    for i in range(NB_SITES):
        if sous_reseau[i] == CU:
            for j in FIRST_NEIGHBOURS[i]:
                if sous_reseau[j] == ZN:
                    nb_zn_1 += 1
            for j in SECOND_NEIGHBOURS[i]:
                if sous_reseau[j] == ZN:
                    nb_zn_2 += 1
    return 1 - (4 * nb_zn_1 / (8 * NB_SITES)), 1 - (4 * nb_zn_2 / (6 * NB_SITES))
    # We multiply by 4 because we divide by 1/2 with cZn and NB_SITES must be divided by 2 (we only check the Cu proba) so we multiply by 2 instead.

def order_parameter() -> float:
    N_alpha_Cu = np.count_nonzero(sous_reseau[::2] == CU)
    N_beta_Cu = np.count_nonzero(sous_reseau[1::2] == CU)
    return (N_alpha_Cu - N_beta_Cu) / (N_alpha_Cu + N_beta_Cu)

def metropolis(n_until_eq: int, n_after_eq: int, T: float) -> tuple[float, float, float, float]:
    init_at_t0()
    alpha_1 = []
    alpha_2 = []
    eta = []
    energies = []
    E = get_full_energy()
    for n in range(n_until_eq + n_after_eq):
        i = choose_random_atom()
        j = choose_random_neighbour(i)
        E_after = energy_after_flip(E, i, j)
        P = flipping_proba(E, E_after, T)
        is_flipped = flipping_action(P, i, j)
        if is_flipped:
            E = E_after
        if n > n_until_eq:
            a_1, a_2 = cowley()
            energies.append(E)
            eta.append(abs(order_parameter()))
            alpha_1.append(a_1)
            alpha_2.append(a_2)
    return np.mean(energies), np.mean(eta), np.mean(alpha_1), np.mean(alpha_2)


def metropolis_before_eq(n: int, T: float) -> list[float]:
    init_at_t0()
    energies = []
    E = get_full_energy()
    for k in range(n):
        i = choose_random_atom()
        j = choose_random_neighbour(i)
        E_after = energy_after_flip(E, i, j)
        P = flipping_proba(E, E_after, T)
        is_flipped = flipping_action(P, i, j)
        if is_flipped:
            E = E_after
        energies.append(E)
    return energies


if __name__ == "__main__":
    energies = metropolis_before_eq(1000, 1000)
    print(sous_reseau)
    print(energies)

