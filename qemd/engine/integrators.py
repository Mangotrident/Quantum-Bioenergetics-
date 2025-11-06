import numpy as np
from scipy.linalg import expm
def expm_krylov(L, dt, v):
    return expm(dt * L) @ v
def evolve_liouville(L, rho0, t_grid):
    v = rho0.reshape(-1, 1).astype(complex)
    traj = []
    last = t_grid[0]
    for t in t_grid:
        dt = t - last
        if dt > 0:
            v = expm_krylov(L, dt, v)
        traj.append(v.copy())
        last = t
    return traj
