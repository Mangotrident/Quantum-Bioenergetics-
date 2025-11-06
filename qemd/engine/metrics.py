import numpy as np
def ete_from_sink_population(traj, n, sink_site, k_sink, t_grid):
    pop = []
    for v in traj:
        rho = v.reshape((n, n))
        pop.append(float(np.real(rho[sink_site, sink_site])))
    integral = np.trapz(np.array(pop), t_grid)
    return float(k_sink * integral)
def coherence_lifetime(traj, n, threshold=1/np.e, t_grid=None):
    if t_grid is None:
        t_grid = np.arange(len(traj))
    coh = []
    for v in traj:
        rho = v.reshape((n, n))
        off = rho.copy()
        np.fill_diagonal(off, 0.0)
        coh.append(float(np.mean(np.abs(off))))
    coh = np.array(coh)
    c0 = coh[0] if coh[0] > 1e-12 else np.max(coh) + 1e-12
    frac = coh / c0
    idx = np.where(frac < threshold)[0]
    return float(t_grid[idx[0]] - t_grid[0]) if len(idx) else float(t_grid[-1] - t_grid[0])
def qls(ete, tau_c, weights=(0.6, 0.4)):
    return float(weights[0]*ete + weights[1]*(np.arctan(tau_c)/(np.pi/2)))
