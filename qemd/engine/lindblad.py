import numpy as np
def kron(a, b):
    return np.kron(a, b)
def liouvillian(H, dephasing_rates, k_sink, k_loss, sink_site):
    n = H.shape[0]
    I = np.eye(n, dtype=complex)
    Hc = H.astype(complex)
    Lh = -1j * (kron(Hc, I) - kron(I, Hc.T))
    Ld = np.zeros((n*n, n*n), dtype=complex)
    for i, gamma in enumerate(dephasing_rates):
        if gamma == 0:
            continue
        Pi = np.zeros((n, n), dtype=complex); Pi[i, i] = 1.0
        Ld += gamma * (kron(Pi, Pi) - 0.5 * (kron(Pi @ Pi, I) + kron(I, (Pi @ Pi).T)))
    Ps = np.zeros((n, n), dtype=complex); Ps[sink_site, sink_site] = 1.0
    Lsink = k_sink * (kron(Ps, Ps) - 0.5 * (kron(Ps, I) + kron(I, Ps.T)))
    if k_loss > 0:
        Psum = np.eye(n, dtype=complex)
        Lloss = k_loss * (kron(Psum, Psum) - 0.5 * (kron(Psum, I) + kron(I, Psum.T)))
    else:
        Lloss = 0.0
    return (Lh + Ld + Lsink + Lloss).astype(complex)
