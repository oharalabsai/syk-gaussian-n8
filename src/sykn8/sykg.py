"""sykg: the shared N=8 numpy library in GIST's conventions. Import, never copy.

Sources (TW = arXiv:2604.10090 .txt line numbers). Conventions pinned here:
- H = sum_T J_T chi_a chi_b chi_c chi_d over sorted quartets T = (a<b<c<d), labels 1..N, i^(q/2) absorbed into J as
  Eq. 11 is written (TW:240-247). {chi_a, chi_b} = delta_ab (each Majorana carries 1/sqrt2, TW:617-619).
- Gaussian ensemble: fixed K, uniform pruning (TW:203-206), Var J_T = J^2 (q-1)!/(p N^(q-1)), p = K/C(N,4), J = sqrt2
  (TW:197-199, 128-129). Binary: J/sqrt(K) random signs (Eq. 10, TW:218); Eq. 11 = 1/sqrt5 (TW:246).
- One-sided spectra: map S4 (TW:645-654), 16 levels; sectors by Z^4: TW 'odd' = Z^4 = +1, 'even' = Z^4 = -1
  (pinned by Eq. 11: r_odd 0.537, r_even 0.479, TW:258-259). r: raw levels, no trim.
- Two-sided: register of N qubits (leftmost kron factor = qubit 1), maps S2 (FIX) / S3 (SYM) (TW:621-640, read
  from the rendered PDF page 11). State axes (P, Q, register); T = register qubit 1 measured directly (TW:296-300).
"""
from __future__ import annotations

import itertools
from functools import lru_cache
from math import comb, factorial, sqrt

import numpy as np

N8 = 8
Q4 = 4
J_SCALE = sqrt(2.0)          # TW:128-129
BETA = 3.0                   # TW:128-129
MU = 12.0                    # |mu|, TW:1095
T0_FIX = 1.8                 # TW:346-350
R_ODD_EQ11, R_EVEN_EQ11, R_TOL = 0.537163, 0.478812, 0.0005   # TW:258-259 (3-decimal rounding)


# ============================ term sets and ensembles ============================================================
@lru_cache(maxsize=None)
def _terms_all(N: int) -> tuple:
    return tuple(itertools.combinations(range(1, N + 1), Q4))


def terms_all(N: int = N8) -> list:
    """All C(N,4) sorted quartets (1-based labels), lexicographic. N=8: 70."""
    return list(_terms_all(N))


def term_index(N: int = N8) -> dict:
    return {t: i for i, t in enumerate(_terms_all(N))}


def eq11():
    """Eq. 11 (TW:240-247): binary K=10, J = 1/sqrt5, signs as printed."""
    terms = [(1, 2, 4, 7), (1, 2, 6, 8), (1, 3, 4, 5), (1, 3, 7, 8), (1, 4, 5, 6),
             (1, 5, 6, 7), (2, 3, 4, 8), (2, 3, 5, 7), (2, 4, 6, 7), (3, 5, 6, 8)]
    J = np.array([-1, 1, -1, 1, -1, 1, -1, 1, -1, 1], float) / sqrt(5.0)
    return terms, J


# Pinned Lie-1 term order for Eq. 11 (indices into eq11()'s printed list): our reconstruction of the Eq. 11 compile order,
# the unique one of all 1568 Lie-1 product classes within 0.001 of six reference Trotter numbers.
EQ11_PDF_ORDER = (9, 3, 5, 1, 4, 8, 7, 6, 2, 0)

# Term orders arXiv:2604.10090 used for Eq. 11 (sent by its authors, 2026-10-23): indices into eq11()'s printed
# list, order[0] applied first, (H_L, H_R). The backward Lie-1 step (same_order or adjoint) was not stated.
EQ11_GIST_ORDERS = {"fix": ((8, 5, 6, 1, 9, 3, 4, 2, 7, 0), (8, 5, 6, 4, 9, 3, 1, 7, 2, 0)),
                    "sym": ((3, 8, 7, 5, 0, 1, 2, 9, 6, 4), (3, 7, 8, 1, 2, 5, 6, 0, 4, 9))}


def default_term_order(terms) -> list:
    """pinned Lie-1 / Suzuki-2 term order: Eq. 11's term set -> EQ11_PDF_ORDER (matched by quartet, so any listing
    of Eq. 11 gets it); any other term set -> its list order (term 0 applied first)."""
    tt = [tuple(sorted(x)) for x in terms]; e = eq11()[0]
    if sorted(tt) == sorted(e) and len(tt) == len(e):
        pos = {q: i for i, q in enumerate(tt)}
        return [pos[e[k]] for k in EQ11_PDF_ORDER]
    return list(range(len(tt)))


def gauss_sigma(K: int, N: int = N8, Jscale: float = J_SCALE) -> float:
    """std of J_T at fixed K: Var = J^2 (q-1)!/(p N^(q-1)), p = K/C(N,4) (TW:197-199, 203-206)."""
    p = K / comb(N, Q4)
    return sqrt(Jscale ** 2 * factorial(Q4 - 1) / (p * N ** (Q4 - 1)))


def _pick(K: int, rng, N: int) -> list:
    allt = _terms_all(N)
    idx = np.sort(rng.choice(len(allt), size=K, replace=False))
    return [allt[i] for i in idx]


def draw_gaussian(K: int, rng, N: int = N8, Jscale: float = J_SCALE):
    """Fixed-K uniform pruning, Gaussian couplings. Returns (terms sorted, J array)."""
    return _pick(K, rng, N), rng.normal(0.0, gauss_sigma(K, N, Jscale), size=K)


def draw_binary(K: int, rng, N: int = N8, Jscale: float = J_SCALE):
    """Eq. 10: J/sqrt(K) with random signs, fixed-K uniform pruning."""
    return _pick(K, rng, N), rng.choice([-1.0, 1.0], size=K) * Jscale / sqrt(K)


def draw_bernoulli(p: float, rng, N: int = N8, Jscale: float = J_SCALE):
    """Each quartet kept with prob p (JH Eq. 2.6-2.7 style); variance uses p (GIST normalization). K varies."""
    allt = _terms_all(N)
    keep = np.nonzero(rng.random(len(allt)) < p)[0]
    sig = sqrt(Jscale ** 2 * factorial(Q4 - 1) / (p * N ** (Q4 - 1)))
    return [allt[i] for i in keep], rng.normal(0.0, sig, size=len(keep))


def to_vec(terms, J, N: int = N8) -> np.ndarray:
    """(terms, J) -> length-C(N,4) coupling vector (zeros for absent quartets)."""
    ix = term_index(N); v = np.zeros(comb(N, Q4))
    for t, j in zip(terms, J):
        v[ix[tuple(sorted(t))]] += j
    return v


def from_vec(v, N: int = N8):
    allt = _terms_all(N); nz = np.nonzero(v)[0]
    return [allt[i] for i in nz], np.asarray(v)[nz]


# ============================ one-sided map S4 and spectra =========================================================
_X = np.array([[0, 1], [1, 0]], complex); _Y = np.array([[0, -1j], [1j, 0]], complex)
_Z = np.diag([1.0 + 0j, -1.0]); _I2 = np.eye(2, dtype=complex)


def _kron(ops):
    out = np.ones((1, 1), complex)
    for o in ops:
        out = np.kron(out, o)
    return out


@lru_cache(maxsize=None)
def chi_s4(N: int = N8) -> np.ndarray:
    """Eq. S4: chi_{2m-1} = Z^{m-1} X / sqrt2, chi_{2m} = Z^{m-1} Y / sqrt2 on N/2 qubits. Shape (N, d, d), index 0 = chi_1."""
    n = N // 2; out = []
    for j in range(1, N + 1):
        m = (j + 1) // 2
        out.append(_kron([_Z] * (m - 1) + [_X if j % 2 else _Y] + [_I2] * (n - m)) / sqrt(2))
    return np.array(out)


@lru_cache(maxsize=None)
def term_mats(N: int = N8) -> np.ndarray:
    """(C(N,4), d, d): chi_a chi_b chi_c chi_d for every sorted quartet (Hermitian, square = 1/16)."""
    C = chi_s4(N)
    return np.array([C[a - 1] @ C[b - 1] @ C[c - 1] @ C[d - 1] for a, b, c, d in _terms_all(N)])


@lru_cache(maxsize=None)
def _parity_idx(N: int = N8):
    """indices of the Z^{N/2} = +1 (TW 'odd') and -1 (TW 'even') sectors in the S4 basis."""
    zz = np.diagonal(_kron([_Z] * (N // 2))).real
    return np.nonzero(zz > 0)[0], np.nonzero(zz < 0)[0]


@lru_cache(maxsize=None)
def _sector_mats(N: int = N8):
    T = term_mats(N); po, pe = _parity_idx(N)
    return np.ascontiguousarray(T[:, po][:, :, po]), np.ascontiguousarray(T[:, pe][:, :, pe])


def h_one_sided(terms, J, N: int = N8) -> np.ndarray:
    """2^(N/2) x 2^(N/2) one-sided H in map S4."""
    return np.einsum("k,kij->ij", to_vec(terms, J, N), term_mats(N))


def h_batch(V: np.ndarray, N: int = N8) -> np.ndarray:
    """V: (B, C(N,4)) coupling vectors -> (B, d, d) Hamiltonians."""
    return np.einsum("bk,kij->bij", np.atleast_2d(V), term_mats(N), optimize=True)


def spectra_batch(V: np.ndarray, N: int = N8) -> np.ndarray:
    """(B, C(N,4)) -> (B, 2^(N/2)) full one-sided spectra, ascending."""
    return np.linalg.eigvalsh(h_batch(V, N))


def sector_spectra_batch(V: np.ndarray, N: int = N8):
    """(B, C(N,4)) -> (E_odd (B, d/2), E_even (B, d/2)); odd = Z^{N/2} = +1 sector (TW label)."""
    V = np.atleast_2d(V); To, Te = _sector_mats(N)
    Ho = np.einsum("bk,kij->bij", V, To, optimize=True); He = np.einsum("bk,kij->bij", V, Te, optimize=True)
    return np.linalg.eigvalsh(Ho), np.linalg.eigvalsh(He)


def sector_spectra(terms, J, N: int = N8):
    Eo, Ee = sector_spectra_batch(to_vec(terms, J, N)[None], N)
    return Eo[0], Ee[0]


def spectrum(terms, J, N: int = N8) -> np.ndarray:
    return np.linalg.eigvalsh(h_one_sided(terms, J, N))


def r_values(E, tol: float = 1e-12) -> np.ndarray:
    """gap ratios min(s_i, s_{i+1})/max(...) on raw sorted levels (last axis). 0/0 (double degeneracy) -> 0."""
    E = np.sort(np.asarray(E, float), axis=-1); s = np.diff(E, axis=-1)
    a, b = s[..., :-1], s[..., 1:]
    hi = np.maximum(a, b); lo = np.minimum(a, b)
    return np.where(hi > tol, lo / np.where(hi > tol, hi, 1.0), 0.0)


def r_mean(E, tol: float = 1e-12):
    """<r> of Eq. 9 / S12 on raw levels, no trim, no unfolding. Works on the last axis of a stack."""
    return r_values(E, tol).mean(axis=-1)


def eq11_check() -> dict:
    """V1 r gate: every lh job writes this into metrics.json (known_ok, eq11_r_err)."""
    Eo, Ee = sector_spectra(*eq11())
    ro, re_ = float(r_mean(Eo)), float(r_mean(Ee))
    err = max(abs(ro - R_ODD_EQ11), abs(re_ - R_EVEN_EQ11))
    return {"known_ok": int(err <= R_TOL), "eq11_r_err": err, "r_odd": ro, "r_even": re_}


# ============================ SFF (Eq. 8 / S10 / S11; Fig. 4 windows) ==============================================
def sff_Y(E, t, alpha: float = 3.0) -> np.ndarray:
    """Y(alpha, t) = sum_i exp(-alpha E_i^2 - i E_i t). E (..., n), t (T,) -> (..., T) complex."""
    E = np.asarray(E, float)[..., None]; t = np.asarray(t, float)
    return (np.exp(-alpha * E ** 2) * np.exp(-1j * E * t)).sum(axis=-2)


def sff_filtered(E, t, alpha: float = 3.0, normalize: bool = False) -> np.ndarray:
    """|Y(alpha,t)|^2: UNNORMALIZED (the quantity disorder-averaged in Eq. S11, TW:788-795) unless normalize=True
    (h(alpha,t) of Eq. S10)."""
    Y = sff_Y(E, t, alpha); y2 = np.abs(Y) ** 2
    if normalize:
        y2 = y2 / (np.abs(sff_Y(E, [0.0], alpha)) ** 2)
    return y2


def sff_windows(t_lo: float = 0.1, t_hi: float = 1000.0, n_win: int = 18, per_win: int = 64):
    """[OUR CHOICE, TW unstated] 18 log-spaced windows on [t_lo, t_hi] (Fig. 4, TW:227-232), per_win points each."""
    edges = np.geomspace(t_lo, t_hi, n_win + 1)
    return np.array([np.linspace(edges[i], edges[i + 1], per_win, endpoint=False) for i in range(n_win)])


def sff_time_avg(E, windows=None) -> np.ndarray:
    """raw (alpha=0) SFF |Y(0,t)|^2/|Y(0,0)|^2 averaged inside each window. E (..., n) -> (..., n_win)."""
    W = sff_windows() if windows is None else np.asarray(windows)
    E = np.asarray(E, float); n = E.shape[-1]
    g = np.abs(sff_Y(E, W.ravel(), 0.0)) ** 2 / n ** 2
    return g.reshape(E.shape[:-1] + W.shape).mean(axis=-1)


def sff_scalar(E, windows=None) -> np.ndarray:
    """chaos gate: min over windows of the time-averaged single-sample SFF. Divide by the DENSE mean (v2_bars/l1.json) yourself."""
    return sff_time_avg(E, windows).min(axis=-1)


# ============================ term-set structure: commutation, mirror, d12 =========================================
def _mask(T) -> int:
    m = 0
    for a in T:
        m |= 1 << (a - 1)
    return m


def n_anticommuting_pairs(terms) -> int:
    """quartets anticommute iff they share an odd number of Majoranas. Eq. 11: 11 (TW:248-249)."""
    ms = [_mask(t) for t in terms]
    return sum(bin(a & b).count("1") % 2 for a, b in itertools.combinations(ms, 2))


@lru_cache(maxsize=None)
def _subset_parity(N: int = N8) -> np.ndarray:
    """(2^N, C(N,4)) bool: |S cap T| odd."""
    tm = np.array([_mask(t) for t in _terms_all(N)])
    S = np.arange(2 ** N)[:, None]
    x = S & tm[None, :]
    pc = np.zeros_like(x)
    for b in range(N):
        pc += (x >> b) & 1
    return (pc % 2).astype(bool)


def mirror_sets(terms, N: int = N8) -> list:
    """all S subset [N] (bitmask) with |S cap T| odd for every T: Gamma_S anticommutes with H -> exact E -> -E."""
    ix = term_index(N); cols = [ix[tuple(sorted(t))] for t in terms]
    ok = _subset_parity(N)[:, cols].all(axis=1)
    return [int(s) for s in np.nonzero(ok)[0]]


def mirror_gf2(terms, N: int = N8) -> bool:
    """some Gamma_S anticommutes with every term (GF(2) solvable) -> exact E -> -E of the FULL spectrum.
    SUFFICIENT only: brute force found K=7 term sets with an exact mirror spectrum and no such S, so gates
    must use mirror_exact (GF(2) OR the spectral residual)."""
    return len(mirror_sets(terms, N)) > 0


def mirror_residual(E) -> np.ndarray:
    """max_i |E_i + E_{n-1-i}| on the sorted FULL spectrum (last axis); 0 = exact E -> -E."""
    E = np.sort(np.asarray(E, float), axis=-1)
    return np.abs(E + E[..., ::-1]).max(axis=-1)


def mirror_exact(terms, J, N: int = N8, tol: float = 1e-9) -> bool:
    """V4 gate: exact E -> -E on the full 16-level spectrum by either mechanism."""
    return mirror_gf2(terms, N) or bool(mirror_residual(spectrum(terms, J, N)) < tol)


def near_mirror(terms, J, N: int = N8) -> float:
    """m(H) = min_S sum_{T: |S cap T| even} J_T^2 / sum_T J_T^2 (0 = exact mirror)."""
    ix = term_index(N); cols = [ix[tuple(sorted(t))] for t in terms]
    w = np.asarray(J, float) ** 2
    even = ~_subset_parity(N)[:, cols]
    return float((even @ w).min() / w.sum())


def d12(terms, pair=(1, 2)) -> int:
    """number of quartets NOT containing exactly one of the pair (|S cap T| even, S = pair). Eq. 11: 3."""
    return sum(1 for t in terms if (pair[0] in t) + (pair[1] in t) != 1)


def dmin(terms, N: int = N8) -> int:
    return min(d12(terms, pr) for pr in itertools.combinations(range(1, N + 1), 2))


@lru_cache(maxsize=None)
def _hidden_candidates(N: int = N8) -> np.ndarray:
    """bitmasks S with |S| even, S != empty, S != [N]: Gamma_S commutes with the fermion parity Gamma_[N]."""
    full = 2 ** N - 1
    return np.array([S for S in range(1, full) if bin(S).count("1") % 2 == 0])


def hidden_z2_sets(terms, N: int = N8) -> list:
    """all nontrivial even S with |S cap T| even for every term T: Gamma_S commutes with H and the parity, so it is a
    hidden Z2 that splits each parity sector 4 + 4 (S and its complement give the same symmetry inside a sector)."""
    ix = term_index(N); cols = [ix[tuple(sorted(t))] for t in terms]
    S = _hidden_candidates(N)
    ok = ~_subset_parity(N)[S][:, cols].any(axis=1)
    return [int(x) for x in S[ok]]


def hidden_z2(terms, N: int = N8) -> bool:
    """True if a hidden commuting Z2 Gamma_S exists (GF(2): |S cap T| even for all T, S nontrivial even)."""
    return len(hidden_z2_sets(terms, N)) > 0


def hidden_z2_residual(terms, J, N: int = N8) -> float:
    """min over nontrivial even S of sum_{T: |S cap T| odd} J_T^2 / sum_T J_T^2 (0 = exact hidden Z2)."""
    ix = term_index(N); cols = [ix[tuple(sorted(t))] for t in terms]
    w = np.asarray(J, float) ** 2
    odd = _subset_parity(N)[_hidden_candidates(N)][:, cols]
    return float((odd @ w).min() / w.sum())


def hidden_parity_flip_sets(terms, N: int = N8) -> list:
    """all S with |S| ODD and |S cap T| even for every term: Gamma_S commutes with H and anticommutes with the parity,
    so it maps one parity sector onto the other and the two sector spectra are identical."""
    ix = term_index(N); cols = [ix[tuple(sorted(t))] for t in terms]
    S = np.array([x for x in range(1, 2 ** N) if bin(x).count("1") % 2 == 1])
    ok = ~_subset_parity(N)[S][:, cols].any(axis=1)
    return [int(x) for x in S[ok]]


def hidden_parity_flip(terms, N: int = N8) -> bool:
    return len(hidden_parity_flip_sets(terms, N)) > 0


def sector_spectra_equal(terms, J, N: int = N8, tol: float = 1e-9) -> bool:
    """spectral check: odd- and even-sector spectra identical. hidden_parity_flip implies it (theorem); the converse
    FAILS at small K (19 of 3000 random K 3..12 sets, all K 5..8, mostly odd-S mirror sets), so gates use
    parity_flip_exact."""
    Eo, Ee = sector_spectra(terms, J, N)
    return bool(np.abs(np.sort(Eo) - np.sort(Ee)).max() < tol)


def parity_flip_exact(terms, J, N: int = N8) -> bool:
    """gate: identical parity-sector spectra by either mechanism (GF(2) odd S, or spectral)."""
    return hidden_parity_flip(terms, N) or sector_spectra_equal(terms, J, N)


def sector_reducible(terms, J, N: int = N8, tol: float = 1e-8) -> tuple:
    """brute-force block test per parity sector (odd, even): the sector Hamiltonian is block-diagonal in some fixed
    basis independent of the couplings iff the term matrices have a nontrivial common commutant. Computed as the
    null-space dimension of X -> ([P_k, X])_k over 8x8 complex X for the sector-projected term matrices P_k (couplings
    J only select which terms enter). Returns the commutant dimensions; 1 = irreducible, >= 2 = blocks."""
    To, Te = _sector_mats(N); ix = term_index(N); cols = [ix[tuple(sorted(t))] for t in terms]
    dims = []
    for Tm in (To, Te):
        d = Tm.shape[1]; I = np.eye(d)
        A = np.vstack([np.kron(Tm[k], I) - np.kron(I, Tm[k].T) for k in cols])
        sv = np.linalg.svd(A, compute_uv=False)
        dims.append(int((sv < tol * max(sv.max(), 1.0)).sum() + max(0, d * d - len(sv))))
    return tuple(dims)


# ============================ two-sided: JW maps S2 / S3, TFD, protocol, I_PT ======================================
# Register = 8 qubits, qubit 1 = leftmost kron factor = MSB of the basis index. Eq. S2/S3 (rendered PDF p. 11-12):
# psi_L^1 = Z X, psi_L^2 = Z Y (qubit 2), psi_R^1 = X, psi_R^2 = Y (qubit 1); for j >= 3 psi_L^j = Z^p X, psi_R^j = Z^p Y
# (Z on qubits 1..p, X/Y on qubit p+1), with p below; every Majorana / sqrt2.
JW_EXP = {"fix": (5, 6, 4, 2, 7, 3), "sym": (6, 7, 3, 4, 2, 5)}   # p for j = 3..8
NREG = 8
DREG = 2 ** NREG
PROTOCOLS = ("sym", "fix")


def _mono(label: str):
    """Pauli string (one char per qubit, qubit 1 first) as a monomial: op|x> = phase[x] |perm[x]>."""
    n = len(label); x = np.arange(2 ** n); perm = x.copy(); ph = np.ones(2 ** n, complex)
    for q, c in enumerate(label, start=1):
        b = (x >> (n - q)) & 1
        if c in "XY":
            perm ^= 1 << (n - q)
        if c == "Y":
            ph *= 1j * (1 - 2 * b)
        elif c == "Z":
            ph *= (1 - 2 * b)
    return perm, ph


def _mono_mul(A, B):
    """monomial of A @ B."""
    return A[0][B[0]], B[1] * A[1][B[0]]


def _mono_dense(M):
    d = len(M[0]); out = np.zeros((d, d), complex); out[M[0], np.arange(d)] = M[1]
    return out


def _mono_apply(M, S):
    """M @ s along the last axis of S."""
    out = np.empty_like(S); out[..., M[0]] = S * M[1]
    return out


def jw_labels(protocol: str) -> dict:
    """Pauli labels of psi_L^j, psi_R^j (without the 1/sqrt2) for map S2 ('fix') or S3 ('sym')."""
    lab = {("L", 1): "ZX" + "I" * 6, ("L", 2): "ZY" + "I" * 6, ("R", 1): "X" + "I" * 7, ("R", 2): "Y" + "I" * 7}
    for j, p in zip(range(3, 9), JW_EXP[protocol]):
        lab[("L", j)] = "Z" * p + "X" + "I" * (NREG - p - 1)
        lab[("R", j)] = "Z" * p + "Y" + "I" * (NREG - p - 1)
    return lab


@lru_cache(maxsize=None)
def majoranas(protocol: str) -> dict:
    """{('L'|'R', j): monomial} with the 1/sqrt2 included."""
    return {k: (lambda m: (m[0], m[1] / sqrt(2)))(_mono(v)) for k, v in jw_labels(protocol).items()}


@lru_cache(maxsize=None)
def majoranas_dense(protocol: str) -> dict:
    return {k: _mono_dense(m) for k, m in majoranas(protocol).items()}


def _term_mono(protocol: str, side: str, T):
    M = majoranas(protocol); out = M[(side, T[0])]
    for a in T[1:]:
        out = _mono_mul(out, M[(side, a)])
    return out


@lru_cache(maxsize=None)
def _term_monos(protocol: str, side: str):
    """for all 70 quartets: (perm (70,256), 4*phase (70,256)); 4*chi_a..chi_d is a Hermitian involution."""
    ms = [_term_mono(protocol, side, T) for T in _terms_all(N8)]
    return np.array([m[0] for m in ms]), np.array([4 * m[1] for m in ms])


@lru_cache(maxsize=None)
def i_state(protocol: str) -> np.ndarray:
    """|I>: (psi_L^j + i psi_R^j)|I> = 0 for all j (TW S6). Unique up to phase."""
    D = majoranas_dense(protocol)
    A = np.vstack([D[("L", j)] + 1j * D[("R", j)] for j in range(1, 9)])
    w, v = np.linalg.eigh(A.conj().T @ A)
    if w[0] > 1e-9 or w[1] < 1e-6:
        raise RuntimeError("|I> not unique")
    vec = v[:, 0]; k = int(np.argmax(np.abs(vec)))
    return vec * np.exp(-1j * np.angle(vec[k]))


@lru_cache(maxsize=None)
def v_eig(protocol: str):
    """V = H_int/(qN), H_int = i sum_j psi_L^j psi_R^j (Eq. 3, TW:138-143). <I|V|I> = -1/8 (ground state)."""
    D = majoranas_dense(protocol)
    V = sum(1j * D[("L", j)] @ D[("R", j)] for j in range(1, 9)) / (Q4 * N8)
    w, v = np.linalg.eigh(V)
    I = i_state(protocol)
    if abs(np.real(I.conj() @ V @ I) + 0.125) > 1e-12:
        raise RuntimeError("|I> is not the H_int ground state")
    return w, v


@lru_cache(maxsize=None)
def swap_ql(protocol: str, pair=(1, 2)) -> np.ndarray:
    """Eq. S6 on (Q, register), Q = left kron factor; Dirac mode psi_L^- = (psi_L^i + i psi_L^j)/sqrt2 (Eq. S5, S9)."""
    D = majoranas_dense(protocol); a, b = pair
    m = (D[("L", a)] + 1j * D[("L", b)]) / sqrt(2); p = m.conj().T
    P0 = np.diag([1.0, 0]); P1 = np.diag([0, 1.0]); s01 = np.array([[0, 1.0], [0, 0]]); s10 = s01.T
    SW = np.kron(P0, m @ p) + np.kron(s01, p) + np.kron(s10, m) + np.kron(P1, p @ m)
    if not np.allclose(SW @ SW.conj().T, np.eye(2 * DREG)):
        raise RuntimeError("SWAP_QL not unitary")
    return SW


def _readout_ops(protocol: str, pair=(1, 2)):
    """T readout: the R Dirac mode of the pair as a qubit: X -> sqrt2 psi_R^i, Y -> sqrt2 psi_R^j, Z -> -iXY.
    For pair (1,2) these are the Paulis of register qubit 1 (TW:296-300, Eq. S8)."""
    M = majoranas(protocol); a, b = pair
    X = (M[("R", a)][0], M[("R", a)][1] * sqrt(2)); Y = (M[("R", b)][0], M[("R", b)][1] * sqrt(2))
    XY = _mono_mul(X, Y)
    return [X, Y, (XY[0], -1j * XY[1])]


_PAULI4 = [np.eye(2, dtype=complex), _X, _Y, _Z]


def entropy2(rho) -> np.ndarray:
    """base-2 von Neumann entropy over a stack (..., d, d)."""
    w = np.linalg.eigvalsh(0.5 * (rho + np.conj(np.swapaxes(rho, -1, -2))))
    w = np.clip(w, 0, None)
    return -(np.where(w > 1e-14, w * np.log2(np.where(w > 1e-14, w, 1.0)), 0.0)).sum(axis=-1)


def _apply_eig(w, v, t, S):
    """exp(-i H t) on S (n, 4, 256), t (n,) or scalar, H = v diag(w) v^dag."""
    t = np.broadcast_to(np.asarray(t, float), (S.shape[0],))
    X = (S @ v.conj()) * np.exp(-1j * np.outer(t, w))[:, None, :]
    return X @ v.T


class TwoSided:
    """Eq. 11-style two-sided model for one Hamiltonian (terms, J), one JW map, one injection pair.
    Protocol (TW Fig. 2, Eq. 4-5, S6, S8): |Bell>_PQ |TFD>; e^{+iH_L t0}; SWAP_QL; e^{-iH_L t0}; e^{i mu V};
    e^{-iH_R t1}; I_PT. H_R terms during the t0 stage cancel exactly (they commute with SWAP_QL) and H_L after the
    kick commutes with the readout, so both are omitted, as on hardware."""

    def __init__(self, terms, J, protocol: str = "fix", pair=(1, 2), order=None, beta: float = BETA,
                 tfd_state=None, r_in_t0: bool = False, order_r=None, back: str = "same_order"):
        """order = Trotter term order (None -> default_term_order). r_in_t0=True also Trotterizes H_R in the t0 stage
        (e^{+-iH_tot t0}, same order both sides; then the R products do not cancel under same_order). Default False:
        R absent from the t0 stage, which with EQ11_PDF_ORDER reproduces the PDF Lie-1 curves to 3e-4 (
        True misses by 2.2e-3 in that order). Exact and Suzuki-2 are identical either way."""
        if protocol not in PROTOCOLS:
            raise ValueError(protocol)
        self.terms = [tuple(sorted(t)) for t in terms]; self.J = np.asarray(J, float)
        self.protocol, self.pair, self.beta = protocol, tuple(pair), beta
        self.order = default_term_order(self.terms) if order is None else list(order)
        self.order_r = self.order if order_r is None else list(order_r)   # H_R's own term order (GIST, Oct 23)
        self.back = back                       # Lie-1 backward step: "same_order" or "adjoint" (reversed order)
        self.r_in_t0 = r_in_t0
        ix = term_index(N8); cols = np.array([ix[t] for t in self.terms])
        self.mL = tuple(a[cols] for a in _term_monos(protocol, "L"))
        self.mR = tuple(a[cols] for a in _term_monos(protocol, "R"))
        self.HL = self._dense(self.mL); self.HR = self._dense(self.mR)
        self.wL, self.vL = np.linalg.eigh(self.HL); self.wR, self.vR = np.linalg.eigh(self.HR)
        if tfd_state is None:
            wt, vt = np.linalg.eigh(self.HL + self.HR)
            tfd = vt @ (np.exp(-beta * wt / 4) * (vt.conj().T @ i_state(protocol)))   # Eq. S19
        else:
            tfd = np.asarray(tfd_state, complex)
        self.tfd = tfd / np.linalg.norm(tfd)
        self.SW = swap_ql(protocol, self.pair); self.wV, self.vV = v_eig(protocol)
        self.ro = _readout_ops(protocol, self.pair)

    def set_pair(self, pair):
        """switch the injection/readout pair in place (H, TFD unchanged)."""
        self.pair = tuple(pair); self.SW = swap_ql(self.protocol, self.pair)
        self.ro = _readout_ops(self.protocol, self.pair)
        return self

    def _dense(self, m):
        H = np.zeros((DREG, DREG), complex); cols = np.arange(DREG)
        for k in range(len(self.J)):
            H[m[0][k], cols] += self.J[k] * m[1][k] / 4
        return H

    # ---- evolutions on S (n, 4, 256) ----
    def _product(self, m, t, S, seq):
        """prod over seq of exp(-i J_k (P_k/4) t), seq[0] applied first; P_k = 4 chi..chi, P_k^2 = 1."""
        t = np.broadcast_to(np.asarray(t, float), (S.shape[0],))
        for k in seq:
            th = (self.J[k] * t / 4)[:, None, None]
            S = np.cos(th) * S - 1j * np.sin(th) * _mono_apply((m[0][k], m[1][k]), S)
        return S

    def evolve(self, side: str, t, S, trotter=None):
        """exp(-i H_side t). trotter None = exact; 'lie1' = one Lie step in self.order; 'suz2' = one Suzuki-2 step.
        Backward (t < 0) keeps the SAME term order with negated angles ('same_order'); for 'suz2' that equals
        the exact inverse of the forward step (palindromic). back='adjoint' reverses the Lie-1 order when every t <= 0."""
        m = self.mL if side == "L" else self.mR
        if trotter is None:
            w, v = (self.wL, self.vL) if side == "L" else (self.wR, self.vR)
            return _apply_eig(w, v, t, S)
        o = self.order if side == "L" else self.order_r
        if trotter == "lie1":
            neg = self.back == "adjoint" and bool(np.all(np.asarray(t, float) <= 0))
            return self._product(m, t, S, o[::-1] if neg else o)
        if trotter == "suz2":
            return self._product(m, np.asarray(t, float) / 2, S, o + o[::-1])
        raise ValueError(trotter)

    def init_state(self, n: int = 1) -> np.ndarray:
        S = np.zeros((4, DREG), complex); S[0] = S[3] = self.tfd / sqrt(2)      # |Bell>_PQ (x) |TFD>, index 2p+q
        return np.broadcast_to(S, (n, 4, DREG)).copy()

    def swap(self, S):
        n = S.shape[0]
        return (S.reshape(n, 2, 2 * DREG) @ self.SW.T).reshape(n, 4, DREG)

    def kick(self, mu: float, S):
        return _apply_eig(self.wV, self.vV, -mu, S)                            # e^{i mu V}

    def rho_pt(self, S) -> np.ndarray:
        """(n, 4, 4) rho_PT, P (x) T, from <sigma_P (x) O_b> with O_b the readout-mode Paulis."""
        n = S.shape[0]; Sp = S.reshape(n, 2, 2 * DREG)
        rho = np.zeros((n, 4, 4), complex)
        for b, O in enumerate([None] + self.ro):
            OS = Sp if O is None else _mono_apply(O, S).reshape(n, 2, 2 * DREG)
            Xb = np.einsum("npk,nPk->npP", OS, Sp.conj())                         # (X_b)_{p,P}
            rho += 0.5 * np.einsum("npP,ab->npaPb", Xb, _PAULI4[b]).reshape(n, 4, 4)
        return rho

    def mi(self, S) -> np.ndarray:
        rho = self.rho_pt(S); r4 = rho.reshape(-1, 2, 2, 2, 2)
        rP = np.einsum("npaPa->npP", r4); rT = np.einsum("npapA->naA", r4)
        return entropy2(rP) + entropy2(rT) - entropy2(rho)

    def ipt(self, mu: float, ts, trotter=None, t0: float = T0_FIX) -> np.ndarray:
        """I_PT (Eq. 4, base 2). 'sym': t0 = t1 = t for t in ts; 'fix': t0 fixed, t1 = ts."""
        ts = np.atleast_1d(np.asarray(ts, float))
        if self.protocol == "sym":
            S = self.init_state(len(ts)); tb = ts
        else:
            S = self.init_state(1); tb = np.array([t0])
        S = self.evolve("L", -tb, S, trotter)
        if self.r_in_t0:
            S = self.evolve("R", -tb, S, trotter)
        S = self.swap(S)
        S = self.evolve("L", tb, S, trotter)
        if self.r_in_t0:
            S = self.evolve("R", tb, S, trotter)
        S = self.kick(mu, S)
        if self.protocol == "fix":
            S = np.broadcast_to(S, (len(ts),) + S.shape[1:]).copy()
        S = self.evolve("R", ts, S, trotter)
        return self.mi(S)

    def delta(self, ts, trotter=None, t0: float = T0_FIX, mu: float = MU) -> np.ndarray:
        """dI = I(mu < 0) - I(mu > 0) (Eq. 5)."""
        return self.ipt(-abs(mu), ts, trotter, t0) - self.ipt(abs(mu), ts, trotter, t0)


def ipt_curve(terms, J, protocol, mu, ts, trotter=None, t0=T0_FIX, beta=BETA, pair=(1, 2), tfd_state=None,
              order=None, term_order=None) -> np.ndarray:
    o = term_order if term_order is not None else order
    return TwoSided(terms, J, protocol, pair, o, beta, tfd_state).ipt(mu, ts, trotter, t0)


def delta_ipt(terms, J, protocol, ts, trotter=None, t0=T0_FIX, beta=BETA, pair=(1, 2), tfd_state=None,
              order=None, mu=MU, term_order=None) -> np.ndarray:
    o = term_order if term_order is not None else order
    return TwoSided(terms, J, protocol, pair, o, beta, tfd_state).delta(ts, trotter, t0, mu)


def all_pairs_delta(terms, J, protocol, ts, trotter=None, t0=T0_FIX, beta=BETA, tfd_state=None, order=None,
                    mu=MU, term_order=None) -> dict:
    """{(i, j): dI(ts)} for all 28 injection pairs i < j (Eq. S9); readout mode = the same pair on R
    ([INFERRED] from Eq. S9: TW replaces the Dirac mode of Eq. S5; the teleported mode is the same pair on R)."""
    m = TwoSided(terms, J, protocol, (1, 2), term_order if term_order is not None else order, beta, tfd_state)
    return {pr: m.set_pair(pr).delta(ts, trotter, t0, mu) for pr in itertools.combinations(range(1, 9), 2)}


# ============================ signal-gate helpers (chaos, signal, VQA gates) =====================================================
def peak_window(ts, d, lo: float = 1.0, hi: float = 4.0):
    """(max dI, t at max) over t in [lo, hi] (signal window)."""
    ts = np.asarray(ts, float); w = (ts >= lo) & (ts <= hi)
    i = int(np.argmax(np.where(w, d, -np.inf)))
    return float(d[i]), float(ts[i])


def revival(ts, d, lo: float = 8.0, hi: float = 20.0) -> float:
    """max |dI| over t in [lo, hi] (revival gate, after TW:1027-1028)."""
    ts = np.asarray(ts, float); w = (ts >= lo) & (ts <= hi)
    return float(np.abs(np.asarray(d)[w]).max())


def trotter_eta(ts, d_exact, d_trotter, lo: float = 1.0, hi: float = 4.0):
    """SI S8 (Eq. S24-S25, TW:1292-1301): delta = max A_ex - max A_tr over [lo, hi]; eta = |delta| / max A_ex
    (the abs is ours, the VQA gate: Trotter above exact cannot pass)."""
    a = peak_window(ts, d_exact, lo, hi)[0]; b = peak_window(ts, d_trotter, lo, hi)[0]
    return a - b, abs(a - b) / a


def coupling_typicality(J, K: int, N: int = N8):
    """(sum J^2, KS p-value of J vs N(0, sigma_K^2)) for the chaos gate."""
    from scipy.stats import kstest
    J = np.asarray(J, float)
    return float((J ** 2).sum()), float(kstest(J / gauss_sigma(K, N), "norm").pvalue)


def sumJ2_percentile(J, K: int, N: int = N8) -> float:
    """percentile of sum J^2 inside its fixed-K Gaussian ensemble (sum J^2 / sigma^2 ~ chi^2_K); chaos-gate typical =
    inside [0.05, 0.95]."""
    from scipy.stats import chi2
    return float(chi2.cdf((np.asarray(J, float) ** 2).sum() / gauss_sigma(K, N) ** 2, K))


def jw_map_ok(labels: dict) -> bool:
    """a JW map is usable iff the 16 Majoranas satisfy {psi_a, psi_b} = delta_ab, the R Dirac mode psi_R^1,2 is
    register qubit 1 (T readout) and |I> is unique. Mislabelled maps fail (V4)."""
    D = {k: _mono_dense(_mono(v)) / sqrt(2) for k, v in labels.items()}
    ops = list(D.values())
    for a in range(len(ops)):
        for b in range(a, len(ops)):
            ac = ops[a] @ ops[b] + ops[b] @ ops[a]
            if not np.allclose(ac, np.eye(DREG) if a == b else 0, atol=1e-12):
                return False
    if labels[("R", 1)] != "X" + "I" * 7 or labels[("R", 2)] != "Y" + "I" * 7:
        return False
    A = np.vstack([D[("L", j)] + 1j * D[("R", j)] for j in range(1, 9)])
    w = np.linalg.eigvalsh(A.conj().T @ A)
    return bool(w[0] < 1e-9 and w[1] > 1e-6)


# ============================ rung-0 self tests and planted faults ==============================================
def v4_planted_faults() -> dict:
    """planted faults that sykg-level gates must reject: {name: rejected (True = good)}. Faults about routing
    (SWAPs dropped, seed cherry-pick, approximation_degree < 1) belong to verify.py / the routing code, not here."""
    out = {}; rng = np.random.default_rng(7)
    t11, J11 = eq11()
    # 1. E -> -E symmetric H: every quartet holds exactly one of (1,2) (DG:34); test on the FULL spectrum, and the
    #    baseline Eq. 11 (whose Z^4 = +1 sector alone is +/-E symmetric) must NOT be rejected.
    one = [t for t in terms_all() if (1 in t) + (2 in t) == 1]
    tm = [one[i] for i in sorted(rng.choice(len(one), 24, replace=False))]
    Jm = rng.normal(0, gauss_sigma(24), 24)
    Eo, _ = sector_spectra(t11, J11)
    out["mirror_symmetric_H"] = mirror_exact(tm, Jm) and not mirror_exact(t11, J11) \
        and bool(mirror_residual(Eo) < 1e-9)
    # 2. commuting term set (all pairs commute): 0 anticommuting pairs and degenerate / non-GOE sectors.
    comm = [(1, 2, 3, 4), (1, 2, 5, 6), (1, 2, 7, 8), (3, 4, 5, 6), (3, 4, 7, 8), (5, 6, 7, 8), (1, 3, 5, 7),
            (1, 3, 6, 8), (1, 4, 5, 8), (1, 4, 6, 7)]
    Jc = rng.normal(0, gauss_sigma(10), 10); Eo, Ee = sector_spectra(comm, Jc)
    out["commuting_set"] = n_anticommuting_pairs(comm) == 0 and min(r_mean(Eo), r_mean(Ee)) < 0.45
    # 3. sign-flipped coupling: Eq. 11 with one sign flipped must fail the eq11 r check.
    Jf = J11.copy(); Jf[0] *= -1; Eo, Ee = sector_spectra(t11, Jf)
    out["sign_flipped_coupling"] = max(abs(r_mean(Eo) - R_ODD_EQ11), abs(r_mean(Ee) - R_EVEN_EQ11)) > R_TOL
    # 4. mislabelled JW map: psi_R^3 given X instead of Y (breaks the algebra) and L/R pair (1,2) swapped (wrong T).
    bad = dict(jw_labels("fix")); bad[("R", 3)] = bad[("L", 3)]
    bad2 = dict(jw_labels("fix"))
    bad2[("L", 1)], bad2[("R", 1)] = bad2[("R", 1)], bad2[("L", 1)]
    bad2[("L", 2)], bad2[("R", 2)] = bad2[("R", 2)], bad2[("L", 2)]
    out["mislabelled_jw_map"] = jw_map_ok(jw_labels("fix")) and jw_map_ok(jw_labels("sym")) \
        and not jw_map_ok(bad) and not jw_map_ok(bad2)
    # 5. revival-only signal: a curve whose peak sits at t = 10 fails the [1, 4] peak and the revival gate.
    ts = np.linspace(0, 20, 201); d = 0.3 * np.exp(-((ts - 10) / 0.7) ** 2)
    out["revival_only_signal"] = peak_window(ts, d)[0] < 0.05 and revival(ts, d) > 0.25
    # 6. high-sum J^2 draw: K = 24 couplings scaled by 2 sit above the ensemble's 95th sum-J^2 percentile (the KS
    #    test alone catches only ~50% of such draws at K = 24: measured 0.507 over 2000, so the sum gate carries it).
    Jt = rng.normal(0, gauss_sigma(24), 24) * 2.0
    out["high_sumJ2_draw"] = sumJ2_percentile(Jt, 24) > 0.95
    # 7. plain qubit SWAP(Q, register qubit 2) instead of Eq. S6 (no JW string): must kill the Eq. 11 signal.
    m = TwoSided(t11, J11, "fix"); ts = np.linspace(0, 6, 61); d_ok = m.delta(ts)
    x = np.arange(2 * DREG); q, r = x // DREG, x % DREG; b2 = (r >> (NREG - 2)) & 1
    dst = b2 * DREG + ((r & ~(1 << (NREG - 2))) | (q << (NREG - 2)))
    m.SW = np.zeros((2 * DREG, 2 * DREG)); m.SW[dst, x] = 1.0
    out["plain_swap_ql"] = d_ok.max() > 0.25 and m.delta(ts).max() < 0.05
    # 8. Trotter dI far above exact: eta uses |delta|, so A_tr = 1.5 A_ex gives eta = 0.5.
    ts = np.linspace(0, 6, 61); de = 0.2 * np.exp(-((ts - 2.5) / 0.8) ** 2)
    out["trotter_above_exact"] = trotter_eta(ts, de, 1.5 * de)[1] > 0.4
    return out


def selftest(verbose: bool = True) -> dict:
    """rung 0 (seconds): eq11 r, structure counts, TFD thermal, JW-map invariance, Trotter limits, V4 faults."""
    res = dict(eq11_check()); t, J = eq11(); m = TwoSided(t, J, "fix")
    res["eq11_anticommuting_pairs"] = n_anticommuting_pairs(t); res["eq11_d12"] = d12(t)
    E = spectrum(t, J); w = np.exp(-BETA * E)
    eth = float((E * w).sum() / w.sum())
    res["tfd_energy_err"] = abs(float(np.real(m.tfd.conj() @ m.HL @ m.tfd)) - eth)
    ts = np.linspace(0, 6, 61)
    res["map_invariance_err"] = max(
        float(np.abs(TwoSided(t, J, "fix").delta(ts, tr) - _as(TwoSided(t, J, "sym"), "fix").delta(ts, tr)).max())
        for tr in (None, "lie1", "suz2"))
    res["hidden_z2_check"] = _hidden_z2_brute(200)
    res["v4"] = v4_planted_faults()
    res["all_ok"] = bool(res["known_ok"] and res["eq11_anticommuting_pairs"] == 11 and res["eq11_d12"] == 3
                         and res["tfd_energy_err"] < 1e-10 and res["map_invariance_err"] < 1e-10
                         and res["hidden_z2_check"]["mismatch"] == 0 and not hidden_z2(t)
                         and res["hidden_z2_check"]["mismatch_flip"] == 0 and not hidden_parity_flip(t)
                         and all(res["v4"].values()))
    if verbose:
        print(res)
    return res


def _hidden_z2_brute(n: int = 200, seed: int = 3) -> dict:
    """random Gaussian term sets (K 4..20, plus planted hidden-Z2 sets): hidden_z2 must equal 'some parity sector is
    reducible' (sector_reducible > 1). Planted: S = {1,2} with every term holding both or neither of 1, 2."""
    rng = np.random.default_rng(seed); mism = 0; n_h = 0; mism_f = 0; n_f = 0; n_conv = 0
    both_or_none = [t for t in terms_all() if (1 in t) == (2 in t)]
    no8 = [t for t in terms_all() if 8 not in t]       # S = {8}: |S cap T| = 0 for all T -> planted parity flip
    for k in range(n):
        if k % 4 == 1:
            K = int(rng.integers(4, len(no8) + 1))
            tt = [no8[i] for i in sorted(rng.choice(len(no8), K, replace=False))]
            JJ = rng.normal(0, gauss_sigma(max(K, 1)), K)
        elif k % 4 == 0:
            K = int(rng.integers(4, len(both_or_none) + 1))
            tt = [both_or_none[i] for i in sorted(rng.choice(len(both_or_none), K, replace=False))]
            JJ = rng.normal(0, gauss_sigma(max(K, 1)), K)
        else:
            tt, JJ = draw_gaussian(int(rng.integers(4, 21)), rng)
        h = hidden_z2(tt); red = max(sector_reducible(tt, JJ)) > 1
        n_h += h; mism += h != red
        f = hidden_parity_flip(tt); eq = sector_spectra_equal(tt, JJ)
        n_f += f; mism_f += f and not eq; n_conv += eq and not f
    return {"n": n, "n_hidden": n_h, "mismatch": mism, "n_flip": n_f, "mismatch_flip": mism_f,
            "equal_without_flip": n_conv}


def _as(m: TwoSided, protocol: str) -> TwoSided:
    """run a model built in one map with the other protocol's schedule (map-invariance test only)."""
    m.protocol = protocol
    return m
