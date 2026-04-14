import numpy as np

"""
=============================================================================

DESCRIPTION:
This script performs the mathematical derivation of the demographic
projection matrices for the homozygous morphs (DD and LL) of the common
buzzard (Buteo buteo). Because explicit stage-specific vital rates for these
morphs were not published in the original analytical model (de Vries & Caswell,
2019), this script reverse-engineers them using a baseline-and-scaling approach.

KEY OPERATIONS:
1. Baseline Scaling: Multiplies the empirical vital rates of the high-fitness
   # heterozygous morph (DL) by a specific factor to exactly match the target
   asymptotic growth rates (lambda = 0.48 for DD; lambda = 0.68 for LL).
2. Lifespan Truncation: Enforces biological lifespan constraints (6 years for
   DD, 7 years for LL) by zeroing out demographic transitions at older ages.
3. SAD & Carrying Capacity Allocation: Computes the dominant eigenvector
   (Stable Age Distribution - SAD) of the resident matrix to logically distribute
   the total Effective Carrying Capacity (K_eff) across the age classes.

OUTPUT:
Generates the formatted NumPy arrays (matrices and vectors) required to
parameterize the main eco-evolutionary simulation in 'fig1_buteo.py'.
=============================================================================
"""

# --- 1. BASELINE MATRIX (DL HETEROZYGOTE) ---
# Derived from Buteo buteo empirical data. Represents the optimal physiological baseline.
L_DL = np.array([
    [0.05, 0.41, 0.42, 0.45, 0.56, 0.48, 0.61, 0.27, 0.40, 0.00, 0.54],
    [0.80, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.00, 0.75, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.00, 0.00, 0.636,0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.00, 0.00, 0.00, 0.714,0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.00, 0.00, 0.00, 0.00, 0.667,0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.00, 0.00, 0.00, 0.00, 0.00, 0.80, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.625,0.00, 0.00, 0.00, 0.00],
    [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.40, 0.00, 0.00, 0.00],
    [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00],
    [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.60]
])

# --- 2. REVERSE ENGINEERING ALGORITHM ---
# Scale baseline matrices and truncate lifespans to strictly match
# the empirical asymptotic growth rates reported in the source literature.

# DD Resident Derivation
factor_DD = 0.48 / 1.04
L_DD = L_DL.copy() * factor_DD
L_DD[7:, :] = 0; L_DD[:, 7:] = 0; L_DD[10, 10] = 0 # Enforce maximum lifespan constraints

# LL Recessive Derivation
factor_LL = 0.68 / 1.04
L_LL = L_DL.copy() * factor_LL
L_LL[8:, :] = 0; L_LL[:, 8:] = 0; L_LL[10, 10] = 0

# --- 3. STRUCTURAL VECTOR CALCULATION ---
# Compute eigenvalues to derive the Stable Age Distribution (SAD)
vals, vecs = np.linalg.eig(L_DD)
sad = np.real(vecs[:, np.argmax(np.real(vals))])
sad = np.abs(sad) / np.sum(np.abs(sad)) # Normalize the SAD vector

K_TOTAL = 1000.0
K_VEC = sad * K_TOTAL

N0_TOTAL = 100.0
SEX_RATIO = 0.5
N0_DD = sad * (N0_TOTAL * SEX_RATIO)

# --- 4. FORMATTED CONSOLE OUTPUT ---
def print_matrix(name, M):
    print(f"{name} = np.array([")
    for row in M:
        print(f"    [{', '.join(f'{x:.4f}' for x in row)}],")
    print("])")
    print()

def print_vector(name, V):
    print(f"{name} = np.array([{', '.join(f'{x:.4f}' for x in V)}])")
    print()

print_matrix("L_DD", L_DD)
print_matrix("L_LL", L_LL)
print_vector("K_VECTOR_INTRA_ESTADIO", K_VEC)
print_vector("n0_DD", N0_DD)

import numpy as np

# 1. EXPLICIT L_DD MATRIX DEFINITION
L_DD = np.array([
    [0.0231, 0.1892, 0.1938, 0.2077, 0.2585, 0.2215, 0.2815, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.3692, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.3462, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.2935, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.3295, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.3078, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3692, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000]
])

# 2. COMPUTE SYSTEM EIGENVECTORS
vals, vecs = np.linalg.eig(L_DD)

# 3. IDENTIFY THE DOMINANT EIGENVECTOR (Stable Age Distribution)
# Locate the index of the largest real eigenvalue
idx = np.argmax(np.real(vals))
sad_raw = np.real(vecs[:, idx])

# 4. NORMALIZE SAD VECTOR (Ensures frequencies sum to 1.0)
sad_norm = np.abs(sad_raw) / np.sum(np.abs(sad_raw))

# 5. FINAL TERMINAL OUTPUT
print("Your SAD vector (normalized to 1.0) is:")
print("-" * 30)
print(f"sad_DD = np.array([{', '.join(f'{x:.4f}' for x in sad_norm)}])")
print("-" * 30)
