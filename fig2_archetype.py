#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ==========================================
# 1. THEORETICAL MATRICES (3x3 - Fast-lived Archetype)
# ==========================================

# --- L_aa (RESIDENT) ---
# Moderate fecundity. Parametrized to maintain a stable equilibrium near K.
L_aa = np.array([
    [0.0, 10, 12],
    [0.6, 0.0, 0.0],
    [0.0, 0.4, 0.0]
])

# --- L_AA (ADVANTAGEOUS INVADER) ---
# Explosive fecundity / survival advantage. Drives population growth beyond carrying capacity.
L_AA = np.array([
    [0.0, 10, 12],
    [0.75,  0.0,  0.0],
    [0.0,  0.4,  0.0]
])

# --- L_Aa (HETEROZYGOTE) ---
# Assuming codominance (intermediate advantage phenotype)
L_Aa = (L_AA + L_aa) / 2.0

# ==========================================
# 2. VECTORS & INITIALIZATION
# ==========================================

# Effective Carrying Capacity (K): 1000 individuals per stage
K_VECTOR_INTRA_ESTADIO = np.array([1000.0, 1000.0, 1000.0])

SEX_RATIO = 0.5

# 1. Residents (aa): Initialize near their demographic equilibrium
# Approximate Stable Age Distribution (SAD) = [0.7, 0.2, 0.1]
n0_aa = np.array([700.0, 200.0, 100.0])

# 2. Heterozygotes (Aa): Inject a small cohort of mutants/invaders into age class 0
n0_Aa = np.array([10.0, 0.0, 0.0])

# 3. Advantageous Homozygotes (AA): Initialize at 0 (will emerge via Mendelian segregation)
n0_AA = np.array([0.0, 0.0, 0.0])

# ==========================================
# 3. ECO-EVOLUTIONARY ENGINE
# ==========================================
def separar_F_y_S(Matriz):
    """Separates Projection Matrix into Fecundity (F) and Survival (S) components."""
    F = np.zeros_like(Matriz); S = np.zeros_like(Matriz)
    F[0, :] = Matriz[0, :]
    S[1:, :] = Matriz[1:, :]
    diag = np.diag(Matriz).copy(); diag[0] = 0.0
    np.fill_diagonal(S, diag)
    return F, S

def simular_dinamica(n_AA, n_Aa, n_aa, K_vec, generaciones, sex_ratio = SEX_RATIO):
    # Matrix decomposition
    F_AA, S_AA = separar_F_y_S(L_AA)
    F_Aa, S_Aa = separar_F_y_S(L_Aa)
    F_aa, S_aa = separar_F_y_S(L_aa)

    hist_AA = [n_AA.copy()]
    hist_Aa = [n_Aa.copy()]
    hist_aa = [n_aa.copy()]

    for t in range(generaciones):
        # A. Fecundity (Gamete Pool Aggregation)
        g_AA = F_AA @ n_AA; g_Aa = F_Aa @ n_Aa; g_aa = F_aa @ n_aa
        total_gametes = g_AA.sum() + g_Aa.sum() + g_aa.sum()

        # Calculate allele frequency (p = dominant allele A) assuming panmixia
        if total_gametes < 1e-12: p = 0
        else: p = (g_AA.sum() + 0.5 * g_Aa.sum()) / total_gametes
        q = 1.0 - p

        # B. Births (Hardy-Weinberg Principle)
        b_AA = total_gametes * (p**2)
        b_Aa = total_gametes * (2 * p * q)
        b_aa = total_gametes * (q**2)

        # C. Ecological Regulation (Ricker Scramble Competition)
        n_total_organismos = (n_AA + n_Aa + n_aa) / sex_ratio

        with np.errstate(divide='ignore', invalid='ignore'):
            exponent = - (n_total_organismos / K_vec)
            phi_vec = np.exp(exponent)

        # D. Transition (Density-regulated survival and aging)
        n_AA_next = S_AA @ (n_AA * phi_vec)
        n_Aa_next = S_Aa @ (n_Aa * phi_vec)
        n_aa_next = S_aa @ (n_aa * phi_vec)

        # E. Recruitment (Injection into class 0)
        # n_AA_next[0] = b_AA; n_Aa_next[0] = b_Aa; n_aa_next[0] = b_aa

        # Update system state
        n_AA, n_Aa, n_aa = n_AA_next, n_Aa_next, n_aa_next
        hist_AA.append(n_AA.copy()); hist_Aa.append(n_Aa.copy()); hist_aa.append(n_aa.copy())

    return np.array(hist_AA), np.array(hist_Aa), np.array(hist_aa)

# ==========================================
# 4. EXECUTION
# ==========================================
N_generations = 350
h_AA, h_Aa, h_aa = simular_dinamica(n0_AA, n0_Aa, n0_aa, K_VECTOR_INTRA_ESTADIO, N_generations, SEX_RATIO)
tiempo = np.arange(h_AA.shape[0])

# ==========================================
# 5. VISUALIZATION
# ==========================================
# Calculate Total Individuals (Adjusted by Sex Ratio)
N_AA = h_AA.sum(axis=1) / SEX_RATIO
N_Aa = h_Aa.sum(axis=1) / SEX_RATIO
N_aa = h_aa.sum(axis=1) / SEX_RATIO
N_total = N_AA + N_Aa + N_aa

# Genotype Frequencies
with np.errstate(divide='ignore', invalid='ignore'):
    freq_AA = N_AA / N_total
    freq_Aa = N_Aa / N_total
    freq_aa = N_aa / N_total

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Frequencies
ax1.plot(tiempo, freq_aa, label='Resident (aa)', color='#08519c', lw=3)
ax1.plot(tiempo, freq_Aa, label='Heterozygote (Aa)', color='#a1d99b', lw=3)
ax1.plot(tiempo, freq_AA, label='Invader (AA)', color='#de2d26', lw=3)
ax1.set_title('a) Genotype Frequencies')
ax1.set_ylim(0, 1.05)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Panel B: Abundance (Age-Class specific transients)
h_AA_tot = h_AA / SEX_RATIO
h_Aa_tot = h_Aa / SEX_RATIO
h_aa_tot = h_aa / SEX_RATIO

# Filter near-zero values to avoid log-scale rendering artifacts
limit = 1e-6
h_AA_tot[h_AA_tot < limit] = np.nan
h_Aa_tot[h_Aa_tot < limit] = np.nan
h_aa_tot[h_aa_tot < limit] = np.nan

# Plot dynamics for the 3 age classes
for i in range(3):
    ax2.plot(tiempo, h_aa_tot[:, i], color='#08519c', lw=2, alpha=0.8)
    ax2.plot(tiempo, h_Aa_tot[:, i], color='#a1d99b', lw=2, alpha=0.8)
    ax2.plot(tiempo, h_AA_tot[:, i], color='#de2d26', lw=2, alpha=0.8)

legend_elements = [Line2D([0],[0], color='#08519c', label='Resident (aa)'),
                   Line2D([0],[0], color='#a1d99b', label='Heterozygote (Aa)'),
                   Line2D([0],[0], color='#de2d26', label='Invader (AA)')]

ax2.set_yscale('log')
ax2.set_title('b) Age-Class Abundance (Transient Oscillations)')
ax2.set_xlabel('Generations')
ax2.set_ylabel('Abundance ($log_{10}$)')

# Fix lower limit; allow upper limit to scale automatically to visualize demographic overshoot
ax2.set_ylim(bottom=1e-1)

ax2.legend(handles=legend_elements)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
