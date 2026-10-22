"""diag: single-Hamiltonian chaos diagnostics at beta = 3 (OTOC, two-point function, size winding), one-sided H in map S4 from sykg.

TW definitions (TW = arXiv:2604.10090 .txt lines): two-point G_i(t) = Re<psi_i(0) psi_i(t)>_beta (Eq. 12, S13);
OTOC F^(i,j)(t) = -Re<[psi_i(0), psi_j(t)]^2>_beta (Eq. 13, TW:252-255); <O> = Tr(rho O), rho = e^{-beta H}/Z;
size winding (S40-S44, S49): c_P(t) = Tr(psi_P^dag rho^(1/2) psi_i(t)) / Tr 1, psi_P = 2^(|P|/2) i^(|P|(|P|-1)/2)
prod_{j in P} psi_j, q(l) = sum_{|P|=l} c_P^2, theta(l) = arg q(l) - arg q(1), R(l) = |q(l)| / sum |c_P|^2;
after the kick (S50-S51): e^{i mu V}|P> = e^{i mu |P|/(qN)}|P>, so theta_after(l) = theta(l) + 2 mu (l-1)/(qN).
psi_i(t) = e^{iHt} psi_i e^{-iHt}. {psi_a, psi_b} = delta_ab (sykg.chi_s4).

Frozen criterion: src/sykn8/l5_thresholds.json. passes(terms, J) is the chaos-diagnostics entry point."""
from __future__ import annotations

import itertools
import json
import os
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sykg  # noqa: E402

BETA, N, Q = sykg.BETA, sykg.N8, 4
TS = np.round(np.arange(0.0, 20.0 + 1e-9, 0.02), 10)        # PDF axes run to t = 20
T_WIND = 2.21                                                # PDF panel (6)
LATE = (8.0, 20.0)
THRESH = Path(__file__).resolve().parent / "l5_thresholds.json"   # the frozen diagnostic thresholds


# ---------------------------------------------------------------- spectral pieces
def eig(terms, J):
    E, V = np.linalg.eigh(sykg.h_one_sided(terms, J))
    w = np.exp(-BETA * (E - E.min()))
    return E, V, w / w.sum()


def psi_eig(V) -> np.ndarray:
    """(N, d, d) Majoranas in the H eigenbasis."""
    C = sykg.chi_s4(N)
    return np.einsum("im,kmn,nj->kij", V.conj().T, C, V)


def evolve(A, E, ts):
    """A(t) = e^{iHt} A e^{-iHt} in the eigenbasis: A_mn e^{i(E_m-E_n)t}. A (..., d, d) -> (T, ..., d, d)."""
    ph = np.exp(1j * np.multiply.outer(ts, E[:, None] - E[None, :]))
    return A[None] * ph.reshape((len(ts),) + (1,) * (A.ndim - 2) + ph.shape[1:])


# ---------------------------------------------------------------- correlators
def two_point(E, P, w, ts=TS) -> np.ndarray:
    """(N, T): G_i(t) = Re Tr(rho psi_i psi_i(t)) = Re sum_mn w_m |psi_mn|^2 e^{i(E_n - E_m)t}."""
    a2 = np.abs(P) ** 2                                       # (N, d, d)
    ph = np.exp(1j * np.multiply.outer(ts, E[None, :] - E[:, None]))   # (T, m, n): e^{i(E_n-E_m)t}
    return np.einsum("m,kmn,tmn->kt", w, a2, ph).real


def otoc(E, P, w, ts=TS):
    """F_ii (N, T) and F_ij (28, T) in itertools.combinations order. Uses psi^2 = 1/2:
    F = 1/2 - Re[Tr(rho A B A B) + Tr(rho B A B A)], A = psi_i, B = psi_j(t)."""
    Pt = evolve(P, E, ts)                                     # (T, N, d, d)
    def F(i, j):
        B = Pt[:, j]
        AB = P[i][None] @ B
        BA = B @ P[i][None]
        x = np.einsum("m,tmm->t", w, AB @ AB) + np.einsum("m,tmm->t", w, BA @ BA)
        return 0.5 - x.real
    Fii = np.array([F(i, i) for i in range(N)])
    Fij = np.array([F(i, j) for i, j in itertools.combinations(range(N), 2)])
    return Fii, Fij


@lru_cache(maxsize=None)
def strings():
    """(256, d, d) psi_P in map S4 (sykg basis) and sizes |P|."""
    C = sykg.chi_s4(N); d = C.shape[1]; mats, sizes = [], []
    for bits in range(2 ** N):
        S = [j for j in range(N) if bits >> j & 1]; L = len(S)
        M = np.eye(d, dtype=complex)
        for j in S:
            M = M @ C[j]
        mats.append(2 ** (L / 2) * (1j) ** (L * (L - 1) // 2) * M); sizes.append(L)
    return np.array(mats), np.array(sizes)


def winding(E, V, w, t=T_WIND, mu=-sykg.MU):
    """per insertion i: theta(l), theta_after(l), R(l) at l = 1,3,5,7; c_P from S49 (lab basis)."""
    S, sz = strings()
    C = sykg.chi_s4(N); d = C.shape[0 + 1]
    U = V @ np.diag(np.exp(-1j * E * t)) @ V.conj().T          # e^{-iHt}
    sq = V @ np.diag(np.sqrt(w)) @ V.conj().T                  # rho^(1/2)
    ls = np.array([1, 3, 5, 7]); out = []
    for i in range(N):
        A = sq @ U.conj().T @ C[i] @ U                         # rho^(1/2) psi_i(t)
        c = np.einsum("pji,ji->p", S.conj(), A) / d            # Tr(psi_P^dag A)/d
        q = np.array([np.sum(c[sz == l] ** 2) for l in ls])
        Pl = np.array([np.sum(np.abs(c[sz == l]) ** 2) for l in ls])
        out.append({"theta": unwrap_rel(q), "R": np.abs(q) / Pl, "Pl": Pl,
                    "theta_after": unwrap_rel(q) + 2 * mu * (ls - 1) / (Q * N)})
    return ls, out


def unwrap_rel(q) -> np.ndarray:
    """arg q(l) - arg q(1), unwrapped across l with each step in (-pi, pi] (reproduces the PDF t = 0 line 0, pi, 2pi, 3pi)."""
    a = np.angle(q)
    st = np.diff(a); st = np.pi - np.mod(np.pi - st, 2 * np.pi)
    return np.concatenate([[0.0], np.cumsum(st)])


# ---------------------------------------------------------------- metrics (pre-registered)
def first_cross(t, y, thr, below: bool) -> float:
    idx = np.nonzero(y <= thr if below else y >= thr)[0]
    return float(t[idx[0]]) if len(idx) else float(t[-1])


def metrics(terms, J, ts=TS, curves: bool = False) -> dict:
    E, V, w = eig(terms, J)
    P = psi_eig(V)
    Fii, Fij = otoc(E, P, w, ts)
    G = two_point(E, P, w, ts)
    late = (ts >= LATE[0]) & (ts <= LATE[1])
    tij = np.array([first_cross(ts, f, 0.75, True) for f in Fij])
    tii = np.array([first_cross(ts, f, 0.25, False) for f in Fii])
    tg = np.array([first_cross(ts, g, 0.25, True) for g in G])
    m = {"E1": np.ptp(tij), "E2": np.ptp(tii), "E3": tij.max(), "E4": tij.std() / tij.mean(),
         "L1": Fij[:, late].std(axis=1).mean(), "L2": Fii[:, late].std(axis=1).mean(),
         "L3": np.ptp(Fij[:, late], axis=1).max(), "L4": np.ptp(Fii[:, late], axis=1).max(),
         "L5": Fij[:, late].mean(axis=0).std(), "L6": Fii[:, late].mean(axis=0).std(),
         "otoc_ij_late_mean": Fij[:, late].mean(), "otoc_ii_late_mean": Fii[:, late].mean(),
         "otoc_ij_late_min": Fij[:, late].min(), "otoc_ii_early_max": Fii[:, ts <= 5].max(),
         "g_half_spread": np.ptp(tg), "g_half_max": tg.max(), "g_late_absmax": np.abs(G[:, late]).max(),
         "g_late_absmean": np.abs(G[:, late]).mean(), "g0": G[:, 0].mean(),
         "F_ii0": np.abs(Fii[:, 0]).max(), "F_ij0": Fij[:, 0].mean()}
    ls, wd = winding(E, V, w)
    slope = lambda th: np.polyfit(ls, th, 1)[0]
    lin = lambda th: np.abs(th - np.polyval(np.polyfit(ls, th, 1), ls)).max()
    Rbar = [np.sum(x["R"] * x["Pl"]) / np.sum(x["Pl"]) for x in wd]
    m.update({"wind_Rbar_min": min(Rbar), "wind_Rbar_mean": float(np.mean(Rbar)),
              "wind_Rmin": min(x["R"].min() for x in wd),
              "wind_slope_min": min(slope(x["theta"]) for x in wd), "wind_slope_max": max(slope(x["theta"]) for x in wd),
              "wind_after_slope_max": max(slope(x["theta_after"]) for x in wd),
              "wind_lin_max": max(lin(x["theta"]) for x in wd),
              "wind_theta7_mean": float(np.mean([x["theta"][-1] for x in wd])),
              "wind_after7_mean": float(np.mean([x["theta_after"][-1] for x in wd])),
              "wind_reversed_all": float(all(slope(x["theta"]) > 0 > slope(x["theta_after"]) for x in wd))})
    m = {k: float(v) for k, v in m.items()}
    if curves:
        m["_curves"] = {"t": ts.tolist(), "Fii": Fii.tolist(), "Fij": Fij.tolist(), "G": G.tolist(),
                        "wind": [{k: v.tolist() for k, v in x.items()} for x in wd]}
    return m


# ---------------------------------------------------------------- frozen criterion (chaos diagnostics)
def load_thresholds(path: Path = THRESH) -> dict:
    return json.loads(Path(path).read_text())


def verdict(m: dict, th: dict) -> dict:
    """th['rules']: list of {metric, max} or {metric, min}; every rule must hold. Returns per-rule booleans."""
    res = {}
    for r in th["rules"]:
        v = m[r["metric"]]
        res[r["name"]] = bool(("max" not in r or v <= r["max"]) and ("min" not in r or v >= r["min"]))
    return res


def passes(terms, J, th: dict | None = None) -> dict:
    """chaos-diagnostics entry: {'pass': bool, 'otoc': bool, 'winding': bool, 'two_point': bool, 'rules': {...}, 'metrics': {...}}."""
    th = th or load_thresholds()
    m = metrics(terms, J)
    v = verdict(m, th)
    grp = {g: all(v[r["name"]] for r in th["rules"] if r["group"] == g) for g in ("otoc", "winding", "two_point")}
    return {"pass": all(grp.values()), **grp, "rules": v, "metrics": m}
