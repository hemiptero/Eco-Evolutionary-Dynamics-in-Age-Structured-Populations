"""
appendix_B_bifurcation.py
-------------------------
Standalone script to compute the local stability (Jacobian matrix and
dominant eigenvalues) of the age-structured population under Ricker regulation.
Demonstrates the analytical period-2 bifurcation threshold using root-finding.
"""

import numpy as np
import scipy.linalg as la
import scipy.optimize as opt
import matplotlib.pyplot as plt

def evaluate_local_stability():
    K_eff = 1000.0
    advantages = np.linspace(0.0, 0.40, 100)
    dominant_eigenvalues = []

    print("--- LOCAL STABILITY ANALYSIS (JACOBIAN) ---")

    for adv in advantages:
        S_AA = np.array([[0, 0, 0], [0.60 + adv, 0, 0], [0, 0.40, 0]])
        F_AA = np.array([[0, 10, 12], [0, 0, 0], [0, 0, 0]])

        # 1. Encontrar el punto fijo numérico (N*) usando un solver algebraico
        # Definimos la función del sistema: f(N) = N_siguiente - N = 0
        def system_eq(N):
            # Prevenir valores negativos en el solver
            N = np.maximum(N, 0)
            phi = np.exp(-np.sum(N) / K_eff)
            N_next = S_AA @ (N * phi) + F_AA @ N
            return N_next - N

        # Adivinanza inicial (initial guess) basada en la capacidad de carga
        guess = np.array([K_eff, K_eff * 0.6, K_eff * 0.6 * 0.4])

        # Encontrar la raíz exacta
        sol = opt.root(system_eq, guess, method='hybr')
        N_star = np.maximum(sol.x, 0) # Asegurar que no haya poblaciones negativas

        total_N_star = np.sum(N_star)

        # 2. Calcular elementos del Jacobiano en el equilibrio exacto
        phi_star = np.exp(-total_N_star / K_eff)
        dphi_dN = -(1.0 / K_eff) * phi_star

        J = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                term_S = S_AA[i, j] * phi_star
                term_density = np.dot(S_AA[i, :], N_star) * dphi_dN
                term_F = F_AA[i, j]
                J[i, j] = term_S + term_density + term_F

        # 3. Calcular el eigenvalor dominante
        eigenvals = la.eigvals(J)
        dom_eig = np.max(np.abs(eigenvals))
        dominant_eigenvalues.append(dom_eig)

    # Identificar el punto de cruce exacto (donde el eigenvalor cruza 1.0)
    crossing_idx = np.where(np.array(dominant_eigenvalues) > 1.0)[0][0]
    critical_adv = advantages[crossing_idx] * 100

    print(f"Mathematical Bifurcation Point Found At: {critical_adv:.2f}% Advantage")

    # Visualización
    plt.figure(figsize=(8, 5))
    plt.plot(advantages * 100, dominant_eigenvalues, color='indigo', lw=2)
    plt.axhline(y=1.0, color='red', linestyle='--', label='Stability Threshold $(|\\lambda_{dom}| = 1)$')
    plt.axvline(x=critical_adv, color='gray', linestyle=':', label=f'Exact Bifurcation Onset (~{critical_adv:.1f}%)')

    plt.title('Analytical Local Stability Analysis\n(Root-Finding Decoupled Jacobian)', weight='bold')
    plt.xlabel('Early-Life Survival Advantage (%)')
    plt.ylabel('Dominant Eigenvalue $(|\\lambda_{dom}|)$')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('appendix_B_eigenvalues_final.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    evaluate_local_stability()
