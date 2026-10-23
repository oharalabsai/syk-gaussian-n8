"""vqa: GIST's variational TFD (TW:1176-1204) fitted to the exact TFD (Eq. S19, TW:1158), per Hamiltonian and JW map.

Conventions (sykg owns the Hamiltonian convention, this file owns the ansatz):
- Register: 8 qubits, qubit 1 = leftmost kron factor. Maps S2 (FIX) / S3 (SYM) read from the rendered PDF (TW:621-640):
  psi_L^1,2 = Z1 X2, Z1 Y2; psi_R^1,2 = X1, Y1; psi_L^j = Z^(e_j) X, psi_R^j = Z^(e_j) Y on qubit e_j + 1 (j >= 3),
  e_3..e_8 = S2 (5,6,4,2,7,3), S3 (6,7,3,4,2,5). Each Majorana carries 1/sqrt2 (TW:617-619).
- |I>: (psi_L^j + i psi_R^j)|I> = 0 for all j (TW:1162-1164) = Bell-type pair on qubits (1,2) x |0>^6.
- TFD = e^(-beta (H_L + H_R)/4)|I>, normalized (Eq. S19), beta = 3, H_R = H_L terms on psi_R.
- Ansatz (TW:1184-1190): |I> then 6 layers of Rx(thx) Rz(thz) per qubit (operator order: Rz acts first, `order="zx"`)
  with a linear CNOT chain cx(q,q+1), q = 1..7, after each of the first 5 layers: 96 parameters, 35 CNOT.
  Variants: start "H" (|+>^8), order "xz", chain "up".
- F = |<TFD|psi(theta)>|^2 (Eq. S21); gradient by the adjoint method (exact, one backward sweep).
Public: tfd(terms, J, jw_map), fit(terms, J, jw_map, ...) -> (theta, F), fidelity(theta, target, ...)."""
from __future__ import annotations

import os
from functools import lru_cache, reduce

import numpy as np
from scipy.optimize import minimize

NQ, NL, NCX = 8, 6, 5
BETA = 3.0
ORDER, CHAIN = "xz", "rdown"    # pinned: qiskit TwoLocal([rx,rz], cx, linear), reversed index
MAPS = {"S2": (5, 6, 4, 2, 7, 3), "S3": (6, 7, 3, 4, 2, 5)}      # e_j for j = 3..8 (TW:621-640, PDF p.12)
MAPS.update(fix=MAPS["S2"], sym=MAPS["S3"])                     # sykg protocol names
_X = np.array([[0, 1], [1, 0]], complex)
_Y = np.array([[0, -1j], [1j, 0]], complex)
_Z = np.diag([1.0 + 0j, -1.0])
_I = np.eye(2, dtype=complex)


# ============================ Majoranas, |I>, exact TFD ===========================================================
def _string(nz: int, P: np.ndarray) -> np.ndarray:
    return reduce(np.kron, [_Z] * nz + [P] + [_I] * (NQ - nz - 1)) / np.sqrt(2)


@lru_cache(maxsize=None)
def majoranas(jw_map: str) -> tuple:
    """(psi_L[1..8], psi_R[1..8]) as 256x256 arrays (index 0 = Majorana 1)."""
    e = (1, 1) + MAPS[jw_map]
    L = [_string(e[0], _X), _string(e[1], _Y)] + [_string(k, _X) for k in e[2:]]
    R = [_string(0, _X), _string(0, _Y)] + [_string(k, _Y) for k in e[2:]]
    return tuple(L), tuple(R)


@lru_cache(maxsize=None)
def i_state(jw_map: str) -> np.ndarray:
    """|I>: the unique null vector of sum_j A_j^dag A_j, A_j = psi_L^j + i psi_R^j (TW:1162-1164)."""
    L, R = majoranas(jw_map)
    M = sum((l + 1j * r).conj().T @ (l + 1j * r) for l, r in zip(L, R))
    w, v = np.linalg.eigh(M)
    assert w[0] < 1e-10 < w[1], f"|I> not unique: {w[:3]}"
    s = v[:, 0]
    return s * np.exp(-1j * np.angle(s[np.argmax(np.abs(s))]))


def hamiltonians(terms, J, jw_map: str):
    """(H_L, H_R): sum_T J_T psi_a psi_b psi_c psi_d, i^(q/2) absorbed in J (sykg convention)."""
    L, R = majoranas(jw_map)
    HL = sum(c * L[a - 1] @ L[b - 1] @ L[cc - 1] @ L[d - 1] for c, (a, b, cc, d) in zip(J, terms))
    HR = sum(c * R[a - 1] @ R[b - 1] @ R[cc - 1] @ R[d - 1] for c, (a, b, cc, d) in zip(J, terms))
    return HL, HR


def tfd(terms, J, jw_map: str, beta: float = BETA) -> np.ndarray:
    """Exact TFD, Eq. S19: e^(-beta H_tot/4)|I> normalized (the 2^(N/2)/sqrt Z prefactor is the norm)."""
    HL, HR = hamiltonians(terms, J, jw_map)
    w, v = np.linalg.eigh(HL + HR)
    s = v @ (np.exp(-beta * (w - w.min()) / 4) * (v.conj().T @ i_state(jw_map)))
    return s / np.linalg.norm(s)


# ============================ ansatz simulator (adjoint gradient) =================================================
_ZS = np.array([[1 - 2 * ((i >> (NQ - 1 - q)) & 1) for i in range(2 ** NQ)] for q in range(NQ)], float)  # Z_q signs


CHAINS = {"down": [(q, q + 1) for q in range(NQ - 1)],               # cx(1,2), cx(2,3), .., cx(7,8)
          "up": [(q, q + 1) for q in range(NQ - 2, -1, -1)],         # cx(7,8), .., cx(1,2)
          "rdown": [(q + 1, q) for q in range(NQ - 2, -1, -1)],      # cx(8,7), .., cx(2,1): qiskit 'linear', reversed index
          "rup": [(q + 1, q) for q in range(NQ - 1)]}                # cx(2,1), .., cx(8,7): qiskit 'reverse_linear'


@lru_cache(maxsize=None)
def _chain_perm(chain: str) -> np.ndarray:
    """psi_new = psi[perm] for one linear CNOT layer CHAINS[chain] = (control, target) pairs (0-based) in time order."""
    idx = np.arange(2 ** NQ)
    for c, t in reversed(CHAINS[chain]):      # U = P_last..P_first: perm = f_first o .. o f_last, built inside out
        idx = idx ^ (((idx >> (NQ - 1 - c)) & 1) << (NQ - 1 - t))     # qubit 1 = MSB
    return idx


def start_state(start: str, jw_map: str) -> np.ndarray:
    if start == "I":
        return i_state(jw_map).copy()
    return np.full(2 ** NQ, 2 ** (-NQ / 2), complex)       # "H": |+>^8 (variant start)


def _rx_layer(psi: np.ndarray, th: np.ndarray) -> np.ndarray:
    t = psi.reshape((2,) * NQ)
    for q in range(NQ):
        c, s = np.cos(th[q] / 2), -1j * np.sin(th[q] / 2)
        a, b = np.take(t, 0, axis=q), np.take(t, 1, axis=q)
        t = np.stack([c * a + s * b, s * a + c * b], axis=q)
    return t.reshape(-1)


def _xflip(psi: np.ndarray, q: int) -> np.ndarray:
    return np.flip(psi.reshape((2,) * NQ), axis=q).reshape(-1)


def _rz_phase(th: np.ndarray) -> np.ndarray:
    return np.exp(-0.5j * (th @ _ZS))


def simulate(theta: np.ndarray, psi0: np.ndarray, order: str = ORDER, chain: str = CHAIN) -> np.ndarray:
    th = np.asarray(theta, float).reshape(NL, NQ, 2)        # [..., 0] = thx, [..., 1] = thz
    perm, psi = _chain_perm(chain), psi0.copy()
    for l in range(NL):
        if order == "zx":
            psi = _rx_layer(psi * _rz_phase(th[l, :, 1]), th[l, :, 0])
        else:
            psi = _rx_layer(psi, th[l, :, 0]) * _rz_phase(th[l, :, 1])
        if l < NCX:
            psi = psi[perm]
    return psi


def infidelity_grad(theta, psi0, target, order: str = ORDER, chain: str = CHAIN):
    """(1 - F, dF/dtheta * -1) with F = |<target|psi(theta)>|^2, adjoint sweep."""
    th = np.asarray(theta, float).reshape(NL, NQ, 2)
    perm = _chain_perm(chain)
    inv = np.argsort(perm)
    phi = simulate(th, psi0, order, chain)
    a = np.vdot(target, phi)
    lam = target.copy()
    g = np.zeros((NL, NQ, 2), complex)

    def gx(l):
        g[l, :, 0] = [-0.5j * np.vdot(lam, _xflip(phi, q)) for q in range(NQ)]

    def gz(l):
        g[l, :, 1] = -0.5j * (_ZS @ (lam.conj() * phi))

    for l in range(NL - 1, -1, -1):
        if l < NCX:
            phi, lam = phi[inv], lam[inv]
        if order == "zx":
            gx(l)
            phi, lam = _rx_layer(phi, -th[l, :, 0]), _rx_layer(lam, -th[l, :, 0])
            gz(l)
            c = _rz_phase(-th[l, :, 1])
            phi, lam = phi * c, lam * c
        else:
            gz(l)
            c = _rz_phase(-th[l, :, 1])
            phi, lam = phi * c, lam * c
            gx(l)
            phi, lam = _rx_layer(phi, -th[l, :, 0]), _rx_layer(lam, -th[l, :, 0])
    dF = 2 * np.real(np.conj(a) * g)
    return 1.0 - abs(a) ** 2, -dF.reshape(-1)


def fidelity(theta, target, jw_map: str, start: str = "I", order: str = ORDER, chain: str = CHAIN) -> float:
    return float(abs(np.vdot(target, simulate(theta, start_state(start, jw_map), order, chain))) ** 2)


def fit_one(theta0, target, jw_map: str, start: str = "I", order: str = ORDER, chain: str = CHAIN,
            maxiter: int = 5000) -> tuple:
    """One L-BFGS-B run from theta0 -> (theta, F, info)."""
    psi0 = start_state(start, jw_map)
    r = minimize(infidelity_grad, np.asarray(theta0, float).reshape(-1), args=(psi0, target, order, chain),
                 jac=True, method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-15, "gtol": 1e-9})
    return r.x, 1.0 - float(r.fun), {"nit": int(r.nit), "status": int(r.status), "gmax": float(np.abs(r.jac).max())}


def fit_homotopy(theta0, terms, J, jw_map: str, beta: float = BETA, steps: int = 7, start: str = "I",
                 order: str = ORDER, chain: str = CHAIN, maxiter: int = 3000) -> tuple:
    """Track the optimum along beta' = 0 .. beta (target |I> at beta' = 0), warm-starting each step."""
    th = np.asarray(theta0, float)
    for b in np.linspace(0.0, beta, steps):
        th, F, info = fit_one(th, tfd(terms, J, jw_map, b), jw_map, start, order, chain, maxiter)
    return th, F, info


def fit_path(theta0, terms0, J0, terms1, J1, jw_map: str, steps: int = 5, maxiter: int = 3000) -> tuple:
    """Track the optimum along H(s) = (1-s) H0 + s H1, s = 0 .. 1 (H0 = Eq. 11 with theta0 its optimum)."""
    u = sorted(set(map(tuple, terms0)) | set(map(tuple, terms1)))
    a = np.array([dict(zip(map(tuple, terms0), J0)).get(t, 0.0) for t in u])
    b = np.array([dict(zip(map(tuple, terms1), J1)).get(t, 0.0) for t in u])
    th, F, info = np.asarray(theta0, float), 0.0, {}
    for s in np.linspace(0.0, 1.0, steps + 1)[1:]:
        th, F, info = fit_one(th, tfd(u, (1 - s) * a + s * b, jw_map), jw_map, maxiter=maxiter)
    return th, F, info


THETA_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vqa_eq11_theta.json")


def eq11_theta(jw_map: str):
    """Eq. 11's converged optimum per map (from the VQA Eq. 11 job; None until committed)."""
    if not os.path.exists(THETA_FILE):
        return None
    import json
    return np.array(json.load(open(THETA_FILE))[jw_map]["theta"], float)


def _fit_worker(task):
    """task = (kind, theta0 or seed, terms, J, jw_map, maxiter, path_from) -> (kind, theta, F)."""
    kind, x, terms, J, jw_map, maxiter, path_from = task
    if kind == "path":
        th, F, _ = fit_path(x, *path_from, terms, J, jw_map, maxiter=maxiter)
    else:
        rng = np.random.default_rng(x if kind != "warm" else 0)
        th0 = x if kind == "warm" else (rng.normal(0, 0.1, NL * NQ * 2) if kind == "small" else
                                        rng.uniform(-np.pi, np.pi, NL * NQ * 2))
        th, F, _ = fit_one(th0, tfd(terms, J, jw_map), jw_map, maxiter=maxiter)
    return kind, [float(v) for v in th], float(F)


def tasks(terms, J, jw_map: str, theta0="eq11", n_starts: int = 14, seed: int = 0, maxiter: int = 3000,
          n_small: int = 4) -> list:
    """Start list for one (Hamiltonian, map): Eq. 11 optimum as a warm start, the same along the Eq. 11 -> H path,
    and n_starts uniform random starts. Picklable, so a batch job can pmap tasks of many candidates at once."""
    if isinstance(theta0, str):
        theta0 = eq11_theta(jw_map)
    out = []
    if theta0 is not None:
        import sykg
        out += [("warm", np.asarray(theta0, float), terms, J, jw_map, maxiter, None),
                ("path", np.asarray(theta0, float), terms, J, jw_map, maxiter, sykg.eq11())]
    out += [("small" if k < n_small else "rand", 1000 * seed + k, terms, J, jw_map, maxiter, None)
            for k in range(n_starts)]
    return out


def fit(terms, J, jw_map: str, theta0="eq11", n_starts: int = 14, seed: int = 0, maxiter: int = 3000,
        workers: int | None = None) -> tuple:
    """Fast fitter for one candidate and map -> (theta, F): warm + path starts from Eq. 11's optimum plus n_starts
    random starts, L-BFGS-B, spawn pool over the box's cores ($LH_CPUS)."""
    jobs = tasks(terms, J, jw_map, theta0, n_starts, seed, maxiter)
    n = workers or int(os.environ.get("LH_CPUS") or os.cpu_count() or 1)
    if n <= 1:
        res = [_fit_worker(j) for j in jobs]
    else:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(min(n, len(jobs)), mp_context=mp.get_context("spawn")) as ex:
            res = list(ex.map(_fit_worker, jobs))
    _, th, F = max(res, key=lambda x: x[2])
    return np.array(th), F


def state(theta, jw_map: str) -> np.ndarray:
    """|psi_VQA> in sykg's register basis (same kron order: pass as sykg.TwoSided(..., tfd_state=state(...)))."""
    return simulate(theta, start_state("I", jw_map))
