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
    for i in range (0, NB_SITES // 2):
        sous_reseau[i*2+1] = ZN

def get_first_neighbours(i: int) -> list[int]:
    return [(i+1) % NB_SITES,
            (i-1) % NB_SITES,
            (i+L*2+1) % NB_SITES, 
            (i+L*2-1) % NB_SITES, 
            (i+2*L**2+1) % NB_SITES, 
            (i+2*L**2-1) % NB_SITES, 
            (i+2*L**2+2*L+1) % NB_SITES, 
            (i+2*L**2+2*L-1) % NB_SITES]


    

def get_second_neighbours(i: int) -> list[int]:
    return [(i+2) % NB_SITES,
            (i-2) % NB_SITES,
            (i+L*2) % NB_SITES, 
            (i-L*2) % NB_SITES, 
            (i+2*L**2) % NB_SITES, 
            (i-2*L**2) % NB_SITES ]

def choose_random_atom() -> int:
    return random.randint(0, NB_SITES-1)


FIRST_NEIGHBOURS = [get_first_neighbours(i) for i in range(NB_SITES)]
SECOND_NEIGHBOURS = [get_second_neighbours(i) for i in range(NB_SITES)]


def choose_random_neighbour(i: int) -> int:
    return FIRST_NEIGHBOURS[i][random.randint(0, len(FIRST_NEIGHBOURS[i]))]

def get_energy_neighbours(i: int) -> float:
    E = 0   
    for j in FIRST_NEIGHBOURS[i]:
        E += J[sous_reseau[i]][sous_reseau[j]]
    for j in SECOND_NEIGHBOURS[i]:
        E += I[sous_reseau[i]][sous_reseau[j]]
    return E


# Energy after virtually flipping the atom at index i
def get_energy_flipped_neighbours(i: int) -> float:
    E = 0   
    for j in FIRST_NEIGHBOURS[i]:
        E += J[(sous_reseau[i] + 1) % 2][sous_reseau[j]]
    for j in SECOND_NEIGHBOURS[i]:
        E += I[(sous_reseau[i] + 1) % 2][sous_reseau[j]]
    return E


def get_full_energy() -> float:
    E = 0
    for i in range(0, NB_SITES):
        E += get_energy_neighbours(i)
    return E/2

def energy_after_flip(E: float, i: int, j: int) -> float:
    delta_before = get_energy_neighbours(i) + get_energy_neighbours(j) # We don't remove the double interaction between i and j as it is also present after flipping.
    delta_after = get_energy_flipped_neighbours(i) + get_energy_flipped_neighbours(j)
    return E - delta_before + delta_after
        
init_at_t0()
print(sous_reseau)
print(get_first_neighbours(431))

E = get_full_energy()
print(E)

print(get_energy_neighbours(300))

print(energy_after_flip(E, 300, 301))

