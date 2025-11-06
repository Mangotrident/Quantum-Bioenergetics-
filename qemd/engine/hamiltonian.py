import numpy as np
def build_hamiltonian(site_energies, couplings, static_disorder_std=0.0, rng=None):
    n = len(site_energies)
    H = np.zeros((n, n), dtype=np.float64)
    if rng is None:
        rng = np.random.default_rng(42)
    eps = np.array(site_energies, dtype=np.float64)
    if static_disorder_std > 0:
        eps = eps + rng.normal(0.0, static_disorder_std, size=n)
    np.fill_diagonal(H, eps)
    for (i, j, Jij) in couplings:
        H[i, j] = Jij
        H[j, i] = Jij
    return H
