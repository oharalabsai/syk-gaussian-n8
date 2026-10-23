"""verify: the frozen verifier (gates: chaos, signal, vqa, diag, cost).

ONE entry: verify(listing) -> Verdict. The listing is the candidate's data only (terms, couplings, per protocol: JW map,
Trotter, term order, operating times, VQA theta, seed lists, shipped routed ISA QASM 3 with its final layout and claimed
2q). Every metric is rebuilt here from that data: physics by sykg / vqa / diag, circuits by circuits (imported, never
copied); no number claimed in the listing is used as a metric. Bars come from frozen JSON files (BAR_FILES), pinned
with the shared code by sha256 in PINS at freeze; a Verdict passes only if every gate passes AND the pins match.
After freeze this file never changes: a needed change is a new frozen version that re-runs on ALL candidates.

Repository copy: the frozen verifier (version 3.1) of the run with its batch-job entry points removed (no job plumbing).
Paths point at this package (v2_bars/ and the modules next to this file); PINS are the sha256 values of the run's
frozen files, seven re-pinned to this copy after comment-only edits (originals noted at PINS); the keys keep
the run's paths as a record, so check_pins reports a mismatch here. Certification ran as a job on the
run's x86_64 platform; `uv run verify` (reproduce.py) recomputes the numbers from the listings."""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from functools import lru_cache
from importlib.metadata import version
from pathlib import Path
from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, ValidationError, model_validator

HERE = Path(__file__).resolve().parent
R = HERE   # package root
sys.path.insert(0, str(HERE))
import diag  # noqa: E402
import sykg  # noqa: E402
import vqa  # noqa: E402

VERSION = "3.1"                                               # 2.0 = mixed maps (S2 FIX / S3 SYM), secondary record
K_RANGE = (24, 30)                                           # 24 = Kc, 30 = our cap
PROTO = {"FIX": ("fix", "S2", "lie1", "lie1"), "SYM": ("sym", "S2", "suzuki2", "suz2")}   # sykg, map, circuits, sykg trot
# one JW map, S2, for both protocols, baselines and candidates alike
MAP_OF = {PROTO[P][0]: PROTO[P][1] for P in PROTO}           # sykg protocol -> its JW map
SYKG_MAP = {"S2": "fix", "S3": "sym"}                        # sykg builds each map under the protocol that owns it


def two_sided(terms, J, p: str, mp: str, order=None, tfd_state=None):
    """sykg model in JW map mp with protocol p's schedule (map invariance; sykg._as, as sykg.selftest uses it)."""
    return sykg._as(sykg.TwoSided(terms, J, SYKG_MAP[mp], order=order, tfd_state=tfd_state), p)
T_PEAK = np.round(np.arange(1.0, 4.0001, 0.05), 4)           # signal window [1, 4] (same grid as the l2 bars)
T_REV = np.round(np.arange(8.0, 20.0001, 0.1), 4)            # signal revival window [8, 20]
SIG8 = ("sym_ex_12", "fix_ex_12", "sym_tr_12", "fix_tr_12", "sym_ex_all", "fix_ex_all", "sym_tr_all", "fix_tr_all")
F_VQA_MIN, COST_FACTOR, N_SEEDS = 0.96, 0.70, 20             # vqa, cost
EQUIV_TOL, IPT_TOL, MAX_ACTIVE, THETA_TOL = 1e-6, 1e-8, 22, 1e-3
SHIP_IPT_TOL = 1e-5      # routed ISA state: opt-3 resynthesis leaves ~1e-7 numerical error (box: 7.2e-7 on Eq. 11 FIX);
                         # consistent with EQUIV_TOL (1 - F <= 1e-6); the logical circuit tie keeps 1e-8
BAR_FILES = {k: f"v2_bars/{k}.json" for k in ("l1", "l2", "l5", "m2e", "vqa_eq11")}   # frozen copies of the run's bar files
# (l1 ensemble bands, l2 typical draws, l5 diagnostic thresholds, m2e cost bars, Eq. 11 VQA fit)
CODE_FILES = ("sykg.py", "circuits.py", "vqa.py", "diag.py")
VERSIONS = {"qiskit": "2.5.2", "qiskit-ibm-runtime": "0.49.0", "numpy": "2.5.2", "scipy": "1.18.1",
            "qiskit_qasm3_import": "0.6.0", "pytket": "2.18.5"}   # pytket: cost-gate tket block synthesis
MACHINE = "x86_64"                                                 # platform.machine() of the lh box image, set at freeze:
                                                             # 2q counts differ across architectures (Mac arm64 579 vs
                                                             # box 577 on one seed), so certification = lh job only
PINS: dict[str, str] = {  # 3.0 freeze; keys are the run's paths (a record)
    # Re-pinned by this repository after comment / docstring-only edits (code unchanged, AST equal; bar files: label
    # strings only, numbers and keys unchanged; m2e.json settings.circuits_py_blob / vqa_py_blob re-pinned to the
    # edited files, run values 266dfb615c9d983e4d969e0456325bc082b93a33 / b3189bebfc552fcd7ad8d001b85bcebd9dc9f057).
    # circuits.py / sykg.py re-pinned 2026-10-23 after adding the opt-in order_r / back arguments (default path unchanged,
    # every reproduce check equal); 22 Oct values 14c1a509... / 77f46c27..., m2e blob d47cf3d.
    # The run's original sha256 values:
    #   scripts/circuits.py: 1be605e9af57d8b2619278ee68b2c6cd79a9ea400cfb9a2221a200addc1d7090
    #   scripts/diag.py: 63f5ff4e5596c72ffc436fe419146dc80e2256e31d38a88fe8bd153a44e2dfc6
    #   scripts/sykg.py: 1cb012260c7aff85e99faa840fb9ec973d29aa8bae8213e587c7bb324ea5da73
    #   scripts/vqa.py: 337f82808406134aea1aa5addafb912db68f4e53c79dcb0d2562c66c82682955
    #   scripts/v2_bars/l5.json: 23b2d54e3e1488c595182bd1f18b3087ea259fcf7565319dcf9d125c6362d719
    #   scripts/v2_bars/l1.json: 772540a8c1d94ba04173ef3ab42b5c37c7a105ba473ce3d019c08c248bfb30b0
    #   scripts/v2_bars/m2e.json: cd8a6b419b944d831eb544e489973b12ac8db7519c52d8096910c7db8509a298
    "scripts/circuits.py": "739b1b6757fd9ca06ae50a4471d30c20878dd5608c9cf0ea13348d1d458857a5",
    "scripts/diag.py": "3f70e6e4558bd61092beae8a273c85bf0fdb829f12a937c40af384458febf44d",
    "scripts/lhio.py": "36dcbc7c414f233e9837afc23f4ce3b904e4c974f4d4bebf54e8b3cbfe1a56a6",
    "scripts/sykg.py": "c2f8644811a0d8f505f349e4f04ce4c48f7de379480ccaa37d520e703efd5bcf",
    "scripts/v2_bars/l1.json": "7bbe1bbbfb434c0c72ff88bbdf2397f0cefc25f2e0be77b95e596eac2f50a986",
    "scripts/v2_bars/l2.json": "31d9d569fa715e67ad09bc507987f5f7e895fd16e5b7ad2b1f858e406f15b668",
    "scripts/v2_bars/l5.json": "59c949018fbd69a1dfa0dd556b9c1eef373e502392e2ba48b707f3b727bfadcb",
    "scripts/v2_bars/m2e.json": "e05c00f5cd34c799cc40a91c6d4665704bad7d9028d0a0d69daaf7e3d7ff8acc",
    "scripts/v2_bars/vqa_eq11.json": "86e100acc7244eb9ae8103369c2b9ac6cb3ab8e6914cc95ad9964440c767d532",
    "scripts/vqa.py": "97eceb3db33dfa53cf5edd369841e4331a5fabcc7d57c53e62b358822aa58276",
}                                                            # sha256 of EVERY bar + code file (relative to R)
VQA_BUDGET = {"n_starts": 720, "seed": 1, "n_small": 64, "maxiter": 3000}   # vqa.tasks: seeds 1000..1719, 1000..1063
# small-angle, no warm start (720-start rule); the SAME budget for Eq. 11's bar and every candidate.
# The vqa gate is judged on the verifier's OWN fit, never on the listed theta (theta is a free knob that moves eta).
J_FLOOR = 1e-3                                               # min |J_T| / sigma_K (near-zero terms fake K)
UNITARY_OPS = {"cz", "rz", "sx", "x", "id", "delay", "barrier"}
BACKEND = "nighthawk"
SYNTH_BASIS = ["cx", "cz", "swap", "rz", "sx", "x", "h", "s", "sdg"]   # simulate the SYNTHESISED circuit (routing trap)


# ============================ schemas =============================================================================
class RouteSettings(BaseModel):
    model_config = ConfigDict(extra="forbid")                # an unknown setting is never silently ignored
    backend: str                 # descriptive; routing always on circuits.backend(BACKEND) (FakeNighthawk)
    opt_level: int
    layout_method: str | None
    routing_method: str | None
    hls: str
    approximation_degree: float
    synth: str
    backward: str
    prep: str
    vqa_seed: int
    term_order: str = ""         # the bars' rule for baselines; a candidate's own order is ProtocolListing.order
    tket: str = ""               # tket block-synthesis descriptor (synth "tket"); must name the pinned pytket
    circuits_py_blob: str = ""   # git blob of circuits.py / vqa.py the bars were routed with: must equal the code
    circuits_py_last_commit: str = ""
    vqa_py_blob: str = ""
    rule: str = ""               # the bars' map rule (reported)
    jw: dict[str, str] = {}      # v3: the bars' JW map per protocol; must equal PROTO's (S2 for both)


class Ship(BaseModel):
    seed: int                    # the routing-seed list A seed whose transpile this is
    final_layout: list[int]      # physical qubit holding logical qubit 0..9 at the end (tq.layout.final_index_layout())
    qasm3: str
    claimed_2q: int


class ProtocolListing(BaseModel):
    jw: str
    trotter: Literal["lie1", "suzuki2"]
    order: list[int]             # Trotter term order (indices into terms), the one the circuit uses
    t0: float
    t1: float
    mu: float
    vqa_theta: list[float]       # 96, vqa.py layout, fitted in this protocol's map
    ship: Ship                   # seeds and settings are never listed: the verifier uses the pinned ones only


class Draw(BaseModel):
    """term_seed: sykg.draw_gaussian(K, default_rng(term_seed)) gives the term set (uniform draw, the run's screening and bank code,
    seed [seed0, K, chunk, i]). coupling_seed None: J is that same draw's J. Else Stage B (the run's screening code,
    draws): J = default_rng(coupling_seed).normal(0, sigma_K, K) on the sorted term set, seed [seed0, K, b, chunk, i]."""
    term_seed: list[int] | None = None     # None only for term_set "optimized"
    coupling_seed: list[int] | None = None


class OptRecord(BaseModel):
    """optimized term set: the optimizer script that produced it, pinned by sha, and its job."""
    script: str                  # path relative to R
    script_sha256: str
    job: str


class Selection(BaseModel):
    """pipeline files (paths relative to R) from which the verifier DERIVES the pooled 1 in M:
    bank = l4_termbank.json (Stage A counts.n_drawn / n_bank, and its sets), stats = Stage B
    bank_stats.json lists ([{bank_id, n, survivors}]) pooled over every screened set; sha = expected sha256 per file."""
    bank: str
    stats: list[str]
    sha: dict[str, str]


class BankSet(BaseModel):
    terms: list[tuple[int, int, int, int]]
    seed: list[int]


class Bank(BaseModel):
    """term-set bank: counts either top level (v2/v3 banks) or under "counts" (v1); one schema."""
    n_drawn: float
    n_bank: int
    sets: list[BankSet]

    @model_validator(mode="before")
    @classmethod
    def _flat(cls, d):
        if isinstance(d, dict) and "n_drawn" not in d and "counts" in d:
            d = {**d, "n_drawn": d["counts"]["n_drawn"]}
        return d


class Listing(BaseModel):
    id: str
    terms: list[tuple[int, int, int, int]]
    J: list[float]
    protocols: dict[Literal["FIX", "SYM"], ProtocolListing]
    term_set: Literal["uniform", "optimized"] = "uniform"   # optimized: verified + reported, never counted
    optimizer: OptRecord | None = None
    draw: Draw | None = None     # ensemble provenance: (terms, J) must rebuild from it exactly


class Metric(BaseModel):
    name: str
    group: str                   # listing, chaos, signal, vqa, diag, cost, pins
    value: float | str | None = None
    bar: str = ""
    passed: bool | None = None   # None = reported only


class Verdict(BaseModel):
    id: str
    version: str
    verify_sha256: str
    passed: bool
    m2_countable: bool           # passed AND term_set "uniform"; distinct term sets = aggregator
    failed: list[str]
    metrics: list[Metric]
    error: str | None = None

    def m(self, name: str) -> Metric:
        return next(x for x in self.metrics if x.name == name)


class Bars(BaseModel):
    files: dict[str, dict]       # BAR_FILES key -> parsed JSON ({} if missing)
    sha: dict[str, str]          # path relative to R -> sha256 (bar files and CODE_FILES)


# ============================ helpers ============================================================================
def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load_bars(files: dict[str, str] = BAR_FILES, root: Path = R) -> Bars:
    got, sha = {}, {}
    for k, rel in files.items():
        p = root / rel
        got[k] = json.loads(p.read_text()) if p.exists() else {}
        if p.exists():
            sha[rel] = sha256(p)
    for rel in CODE_FILES:
        sha[rel] = sha256(root / rel)
    return Bars(files=got, sha=sha)


def get(d, *keys):
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def M(name, group, value=None, bar="", passed=None) -> Metric:
    if isinstance(value, (np.floating, np.integer)):
        value = float(value)
    return Metric(name=name, group=group, value=value, bar=bar, passed=None if passed is None else bool(passed))


def gate(name, group, value, lo=None, hi=None, strict=False) -> Metric:
    """value inside [lo, hi] (None = open side; strict: lo < value). A missing bar or value fails."""
    if value is None or (lo is None and hi is None) or (isinstance(value, float) and math.isnan(value)):
        return M(name, group, value, f"[{lo}, {hi}] (bar or value missing)", False)
    ok = (lo is None or (value > lo if strict else value >= lo)) and (hi is None or value <= hi)
    return M(name, group, value, f"{'(' if strict else '['}{lo}, {hi}]", ok)


def pct(q, v) -> float | None:
    """percentile (0..100) of v in a 101-entry quantile table."""
    return None if q is None else float(np.interp(v, q, np.arange(len(q))))


# ============================ physics values (sykg) ==============================================================
def signal_values(terms, J, orders: dict[str, list[int] | None]) -> dict[str, float]:
    """signal gate: max over [1, 4] of dI, pair (1,2) and the 28-pair mean, exact TFD, exact and hardware Trotter, both
    protocols (FIX t0 = 1.8); revival max |dI| over [8, 20] (exact); t_op = argmax of the exact mu < 0 I_PT."""
    v = {}
    for P, (p, _, _, tr) in PROTO.items():
        m = sykg.TwoSided(terms, J, p, order=orders[p])
        Im = m.ipt(-sykg.MU, T_PEAK)
        v[f"{p}_ex_12"] = float(np.max(Im - m.ipt(sykg.MU, T_PEAK)))
        v[f"{p}_t_op"] = float(T_PEAK[int(np.argmax(Im))])
        v[f"{p}_tr_12"] = float(np.max(m.delta(T_PEAK, tr)))
        for kind, t in (("ex", None), ("tr", tr)):
            ap = sykg.all_pairs_delta(terms, J, p, T_PEAK, t, order=orders[p])
            v[f"{p}_{kind}_all"] = float(np.mean([np.max(d) for d in ap.values()]))
        v[f"{p}_rev"] = sykg.revival(T_REV, m.delta(T_REV))
    return v


class FitShard(BaseModel):
    """one slice [lo, hi) of the verifier's own VQA_BUDGET fit (a fit job, `verify.py --fit`): its best theta and F."""
    key: str                     # sha256 of the fit key (terms, J, map, budget)
    lo: int
    hi: int
    theta: list[float]
    F: float
    job: str = ""


FITS: dict[str, list[FitShard]] = {}                         # supplied shards (main: $V2_FITS), keyed by fit key sha


def fit_key(terms, J, mp: str) -> str:
    return json.dumps({"terms": [list(map(int, t)) for t in terms], "J": [float(j) for j in J], "map": mp,
                       "budget": VQA_BUDGET})


@lru_cache(maxsize=8)
def vqa_best(key: str, lo: int = 0, hi: int | None = None) -> tuple[tuple[float, ...], float]:
    """the verifier's own VQA fit: max-F over the pinned VQA_BUDGET starts [lo, hi) (vqa.tasks / vqa._fit_worker), spawn pool
    over $LH_CPUS. key = fit_key(...) (cached: Eq. 11 as a candidate reuses its bar fit)."""
    d = json.loads(key)
    ts = vqa.tasks([tuple(t) for t in d["terms"]], np.asarray(d["J"], float), d["map"], theta0=None, **d["budget"])[lo:hi]
    n = int(os.environ.get("LH_CPUS") or os.cpu_count() or 1)
    if n > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(min(n, len(ts)), mp_context=mp.get_context("spawn")) as ex:
            res = list(ex.map(vqa._fit_worker, ts, chunksize=4))
    else:
        res = [vqa._fit_worker(t) for t in ts]
    _, th, F = max(res, key=lambda r: r[2])
    return tuple(th), float(F)


def own_fit(terms, J, mp: str) -> tuple[tuple[float, ...], float]:
    """max F over the whole budget: from supplied fit shards if they tile [0, n_starts) exactly, else computed here."""
    k = fit_key(terms, J, mp)
    sh = sorted(FITS.get(hashlib.sha256(k.encode()).hexdigest(), []), key=lambda x: x.lo)
    if sh and sh[0].lo == 0 and sh[-1].hi == VQA_BUDGET["n_starts"] and all(
            a.hi == b.lo for a, b in zip(sh, sh[1:])):
        b = max(sh, key=lambda x: x.F)
        return tuple(b.theta), float(b.F)
    return vqa_best(k)


def eta_value(terms, J, p: str, order, theta) -> tuple[float, float]:
    """vqa gate: (F of the VQA TFD in the protocol's map, eta = |max A_ex - max A_tr| / max A_ex) from that TFD."""
    mp, tr = MAP_OF[p], {"fix": "lie1", "sym": "suz2"}[p]
    st = vqa.state(np.asarray(theta, float), mp)
    F = float(abs(np.vdot(vqa.tfd(terms, J, mp), st)) ** 2)
    m = two_sided(terms, J, p, mp, order=order, tfd_state=st)
    d_ex = m.delta(T_PEAK)
    eta = sykg.trotter_eta(T_PEAK, d_ex, m.delta(T_PEAK, tr))[1]
    return F, (eta if np.max(d_ex) > 0 else float("nan"))       # no signal = no eta (fails closed)


@lru_cache(maxsize=4)
def eq11_ref(theta_key: str, fit: bool) -> dict:
    """Eq. 11's own values with the same settings (PDF term order). vqa-gate bar F / eta from the pinned same-budget Eq. 11
    fit (v2_bars/vqa_eq11.json = vqa_eq11_theta_budget.json: the identical VQA_BUDGET procedure; the box test asserts
    own_fit(Eq. 11) reproduces it), so a certification job does not refit Eq. 11 each time."""
    t, J = sykg.eq11()
    v = signal_values(t, J, {"fix": None, "sym": None})
    th = json.loads(theta_key) if theta_key else {}
    for p, mp in MAP_OF.items():
        if fit and get(th, mp, "theta") is not None:
            v[f"{p}_F"], v[f"{p}_eta"] = eta_value(t, J, p, None, th[mp]["theta"])
    return v


# ============================ checks (each returns list[Metric]) =================================================
def check_pins(B: Bars) -> list[Metric]:
    need = set(BAR_FILES.values()) | set(CODE_FILES)
    bad = sorted(k for k in need if k not in PINS or B.sha.get(k) != PINS[k]) + sorted(set(PINS) - need)
    import platform
    def ver(k):
        try:
            return version(k)
        except Exception:                                    # noqa: BLE001
            return "missing"
    vers = {k: ver(k) for k in VERSIONS}
    return [M("machine", "pins", platform.machine(), MACHINE or "UNFROZEN", platform.machine() == MACHINE),
            M("pins", "pins", "UNFROZEN" if not PINS else (",".join(bad) or "ok"),
              "every bar + code file present, sha256 = PINS",
              bool(PINS) and not bad),
            M("versions", "pins", json.dumps(vers), json.dumps(VERSIONS), vers == VERSIONS)]


def check_listing(L: Listing, B: Bars) -> list[Metric]:
    K, out = len(L.terms), []
    out.append(gate("K", "listing", K, *K_RANGE))
    import circuits as C
    ok = (len(L.J) == K and all(math.isfinite(j) and j != 0 for j in L.J) and list(L.terms) == sorted(set(L.terms))
          and all(list(t) == sorted(set(t)) and 1 <= t[0] and t[3] <= 8 for t in L.terms))
    out.append(M("term_set_id", "listing", hashlib.sha256(json.dumps(sorted(map(list, L.terms))).encode()).hexdigest()[:16],
                 "reported: distinct term sets = distinct ids"))
    out.append(M("terms_wellformed", "listing", K, "sorted list of distinct sorted quartets in 1..8, K finite nonzero J",
                 ok))
    canon = C.greedy_order(sorted(L.terms))                   # the one canonical Trotter order
    out.append(M("ensemble_draw", "listing", L.draw.model_dump_json() if L.draw else None,
                 "(terms, J) rebuild exactly from the draw seeds (else tuned: L6 flagged list)", ensemble_ok(L)))
    seeds_a = get(B.files["m2e"], "seeds_a") or []
    for P, (p, mp, ctr, _) in PROTO.items():
        x = L.protocols.get(P)
        if x is None:
            out.append(M(f"{P}_present", "listing", None, "both protocols listed", False))
            continue
        out += [M(f"{P}_jw", "listing", x.jw, mp, x.jw == mp),
                M(f"{P}_trotter", "listing", x.trotter, ctr, x.trotter == ctr),
                M(f"{P}_order", "listing", str(x.order), f"canonical greedy_order(sorted terms) = {canon}",
                  ok and x.order == canon),
                M(f"{P}_mu", "listing", x.mu, f"{-sykg.MU}", x.mu == -sykg.MU),
                M(f"{P}_t0", "listing", x.t0, "1.8 (FIX) / = t1 (SYM)",
                  x.t0 == sykg.T0_FIX if P == "FIX" else x.t0 == x.t1),
                M(f"{P}_theta", "listing", len(x.vqa_theta), "96 finite",
                  len(x.vqa_theta) == 96 and all(map(math.isfinite, x.vqa_theta))),
                M(f"{P}_ship_seed", "listing", x.ship.seed, "in the pinned routing-seed list A", x.ship.seed in seeds_a)]
    return out


def ensemble_ok(L: Listing) -> bool:
    d, K = L.draw, len(L.terms)
    if d is None:
        return False
    if L.term_set == "optimized":                             # terms from the pinned optimizer; J must still be i.i.d.
        o = L.optimizer
        if o is None or d.coupling_seed is None or not (R / o.script).exists() or sha256(R / o.script) != o.script_sha256:
            return False
    else:
        if d.term_seed is None or len(d.term_seed) != 4 or d.term_seed[1] != K:
            return False
        tA, JA = sykg.draw_gaussian(K, np.random.default_rng(d.term_seed))
        if [tuple(int(a) for a in t) for t in tA] != list(L.terms):
            return False
    if d.coupling_seed is None:
        return bool(np.array_equal(JA, np.asarray(L.J, float)))
    if len(d.coupling_seed) != 5 or d.coupling_seed[1] != K:
        return False
    JB = np.random.default_rng(d.coupling_seed).normal(0.0, sykg.gauss_sigma(K), size=K)
    return bool(np.array_equal(JB, np.asarray(L.J, float)))


def check_selection(L: Listing, sel: Selection | None) -> list[Metric]:
    """reported pooled 1 in M; GATE for a Stage B candidate (coupling_seed set): its term set and Stage A seed are a
    set of the pinned-by-sha bank (no private term-set search)."""
    stage_b = L.term_set == "uniform" and L.draw is not None and L.draw.coupling_seed is not None
    if sel is None:
        return [M("selection_1_in_M", "chaos", None, "no pipeline files given"),
                M("bank_membership", "chaos", None, "Stage B needs the bank", False if stage_b else None)]
    files = [sel.bank, *sel.stats]
    ok = all((R / f).exists() and sha256(R / f) == sel.sha.get(f) for f in files)
    if not ok:
        return [M("selection_sha", "chaos", "MISMATCH or missing", "every selection file = its sha", False)]
    bank = Bank.model_validate_json((R / sel.bank).read_text())
    rows = [r for f in sel.stats for r in json.loads((R / f).read_text())]
    nb, sv = sum(int(r["n"]) for r in rows), sum(int(r["survivors"]) for r in rows)
    n_sets = len({r["bank_id"] for r in rows}) if rows else float(bank.n_bank)   # Stage B screened only these sets
    m_terms = bank.n_drawn / max(float(n_sets), 1.0)
    member = L.term_set == "uniform" and any(
        list(st.terms) == list(L.terms) and L.draw is not None and L.draw.term_seed is not None
        and st.seed == list(L.draw.term_seed) for st in bank.sets)
    return [M("selection_1_in_M", "chaos", m_terms * (nb / max(sv, 1) if rows else 1.0), f"reported: (Stage A {bank.n_drawn:.0f} "
              f"/ {n_sets} sets screened in Stage B) x (Stage B {nb} / {sv}), files {json.dumps(sel.sha)}"),
            M("bank_membership", "chaos", member, "term set + Stage A seed in the bank", member if stage_b else None)]


def git_blob(p: Path) -> str:
    b = p.read_bytes()
    return hashlib.sha1(b"blob %d\0" % len(b) + b).hexdigest()


def check_chaos(L: Listing, B: Bars) -> list[Metric]:
    l1 = B.files["l1"]; D = get(l1, "K", "70"); Kq = get(l1, "K", str(len(L.terms)), "q") or {}
    Eo, Ee = sykg.sector_spectra(L.terms, L.J)
    ro, re_ = float(sykg.r_mean(Eo)), float(sykg.r_mean(Ee))
    dm = get(l1, "dense_sff_scalar_mean")
    sff = float(sykg.sff_scalar(sykg.spectrum(L.terms, L.J))) / dm if dm else None
    out = []
    for k, v in (("r_odd", ro), ("r_even", re_), ("sff_scalar", sff)):
        out.append(gate(k, "chaos", v, *(get(D, f"{k}_c90") or (None, None))))
        out.append(M(f"{k}_pct_ownK", "chaos", pct(Kq.get(k), v) if v is not None else None, "reported"))
    return out


def check_typicality(L: Listing, B: Bars) -> list[Metric]:
    K = len(L.terms); J = np.asarray(L.J, float)
    s, ks = sykg.coupling_typicality(J, K)
    Kb = get(B.files["l1"], "K", str(K)); q = get(Kb, "q") or {}
    nac = float(sykg.n_anticommuting_pairs(L.terms))
    return [M("sumJ2_pct_ownK", "chaos", pct(q.get("sumJ2"), s), "reported (l1 table)"),
            M("anticomm_pairs_pct_ownK", "chaos", pct(q.get("anticomm_pairs"), nac), "reported (l1 table)"),
            gate("min_absJ_over_sigma", "chaos", float(np.abs(J).min() / sykg.gauss_sigma(K)), J_FLOOR),
            gate("sumJ2_pct", "chaos", sykg.sumJ2_percentile(J, K), 0.05, 0.95),
            M("sumJ2", "chaos", s, f"ensemble c90 {get(Kb, 'sumJ2_c90')}"),
            gate("ks_p", "chaos", ks, 0.05, None, strict=True),
            gate("anticomm_pairs", "chaos", nac, get(Kb, "anticomm_pairs_p5")),
            ]


def signal_metrics(v: dict, eq: dict, B: Bars) -> list[Metric]:
    """pure: signal gates on precomputed values (8 comparisons vs max(DENSE median, Eq. 11), revival vs DENSE p95)."""
    F = get(B.files["l2"], "F") or {}
    out = []
    for k in SIG8:
        dm, e = get(F, "dense", "dI_median", k), eq.get(k)
        out.append(gate(k, "signal", v.get(k), None if dm is None or e is None else max(dm, e)))
    for p in ("sym", "fix"):
        out.append(gate(f"{p}_rev", "signal", v.get(f"{p}_rev"), None, get(F, "dense", "rev_p95", p)))
    fe = get(F, "eq11") or {}
    ks = SIG8 + ("sym_rev", "fix_rev")
    diff = max((abs(fe[k] - eq[k]) for k in ks if k in fe and k in eq), default=None)
    out.append(M("eq11_bars_consistent", "signal", diff, "|l2 eq11 - live eq11| <= 1e-6 on the 8 peaks + 2 revivals",
                 diff is not None and diff <= 1e-6 and all(k in fe for k in ks)))
    return out


def eta_metrics(vals: dict, eq: dict) -> list[Metric]:
    """pure: vqa gates; vals/eq carry {p}_F and {p}_eta."""
    out = []
    for p in ("fix", "sym"):
        f = vals.get(f"{p}_F")
        out.append(M(f"{p}_vqa_F_marginal", "vqa", None if f is None else
                     ("MARGINAL" if abs(f - F_VQA_MIN) < 0.002 else "clear"), "reported: |F - 0.96| < 0.002 = MARGINAL"))
        out.append(gate(f"{p}_vqa_F", "vqa", vals.get(f"{p}_F"), F_VQA_MIN))
        out.append(gate(f"{p}_eq11_vqa_F", "vqa", eq.get(f"{p}_F"), F_VQA_MIN))
        out.append(gate(f"{p}_eta", "vqa", vals.get(f"{p}_eta"), None, eq.get(f"{p}_eta")))
    return out


def check_signal_and_eta(L: Listing, B: Bars, fit: bool = True) -> list[Metric]:
    orders = {PROTO[P][0]: (L.protocols[P].order if P in L.protocols else None) for P in PROTO}
    v = signal_values(L.terms, L.J, orders)
    eq = eq11_ref(json.dumps(B.files["vqa_eq11"], sort_keys=True), fit)
    out = signal_metrics(v, eq, B)
    for P, x in L.protocols.items():
        p = PROTO[P][0]
        tread = x.t1 if P == "FIX" else x.t0
        out.append(M(f"{P}_t_operating", "signal", tread, f"exact mu<0 I_PT argmax on T_PEAK = {v[f'{p}_t_op']}",
                     abs(tread - v[f"{p}_t_op"]) < 1e-9))
    ev = {}
    for P, x in L.protocols.items():
        p, mp = PROTO[P][0], PROTO[P][1]
        Fl = eta_value(L.terms, L.J, p, x.order, x.vqa_theta)[0]
        if fit:
            th, _ = own_fit(L.terms, L.J, mp)
            jobs = sorted({x.job for x in FITS.get(hashlib.sha256(fit_key(L.terms, L.J, mp).encode()).hexdigest(), [])})
            out.append(M(f"{P}_vqa_fit_source", "vqa", ",".join(jobs) or "computed in this job",
                         "reported: the verifier's own VQA_BUDGET fit (shards tile [0, n_starts))"))
            ev[f"{p}_F"], ev[f"{p}_eta"] = eta_value(L.terms, L.J, p, x.order, th)
            out.append(gate(f"{P}_listed_theta_F", "vqa", Fl, max(F_VQA_MIN, ev[f"{p}_F"] - THETA_TOL)))
    return out + (eta_metrics(ev, eq) if fit else [M("m2c", "vqa", None, "not evaluated (fit=False)", False)])


def check_diag(L: Listing, B: Bars) -> list[Metric]:
    th = B.files["l5"]
    out = []
    if th:
        d = diag.passes(L.terms, L.J, th)
        out += [M(f"diag_{g}", "diag", None, "l5_thresholds.json", d[g]) for g in ("otoc", "winding", "two_point")]
        out.append(M("otoc_E3", "diag", d["metrics"]["E3"], "gate in l5 rules; DENSE-anchored E3 <= 3.50 is SECONDARY, "
                     "reported only"))
    else:
        out.append(M("diag", "diag", None, "l5_thresholds.json missing", False))
    out.append(M("mirror_exact", "diag", None, "no exact E -> -E (GF(2) or full 16-level spectrum)",
                 not sykg.mirror_exact(L.terms, L.J)))
    out.append(M("no_hidden_z2", "diag", sykg.hidden_z2_residual(L.terms, L.J), "no hidden commuting Z2 Gamma_S (GF(2), "
                 "sykg.hidden_z2; value = residual, 0 = exact)", not sykg.hidden_z2(L.terms)))
    out.append(M("no_parity_flip", "diag", None, "no odd-|S| hidden symmetry and no identical parity-sector spectra "
                 "(sykg.parity_flip_exact)", not sykg.parity_flip_exact(L.terms, L.J)))
    out.append(M("sector_commutant_dims", "diag", str(sykg.sector_reducible(L.terms, L.J)), "reported (1, 1) = irreducible"))
    m = sykg.near_mirror(L.terms, L.J)
    out.append(gate("near_mirror_m", "diag", m, get(B.files["l1"], "K", "70", "near_mirror_p5"), None, strict=True))
    out.append(M("near_mirror_pct_ownK", "diag", pct(get(B.files["l1"], "K", str(len(L.terms)), "q", "near_mirror"), m),
                 "reported (l1 table)"))
    out += [M("d12", "diag", float(sykg.d12(L.terms)), "reported"), M("dmin", "diag", float(sykg.dmin(L.terms)), "reported")]
    return out


# ============================ cost and the shipped circuit (circuits / qiskit) ===================================
def _cand(L: Listing, P: str, x: ProtocolListing, S: RouteSettings, mu: float | None = None, fitted: bool = False):
    """circuits.Candidate from the listing; fitted=False: generic VQA angles (pinned vqa_seed) for counts."""
    import circuits as C
    return C.Candidate(id=L.id, terms=L.terms, J=L.J, jw=x.jw, protocol=P, trotter=x.trotter, t0=x.t0, t1=x.t1,
                       mu=x.mu if mu is None else mu, order=x.order, backward=S.backward, prep=S.prep,
                       vqa_seed=S.vqa_seed, synth=S.synth, theta=list(x.vqa_theta) if fitted else None)


def lower_median(rs: list[tuple[int, int, int]]) -> tuple[int, int, int]:
    """(2q, depth, seed) of sorted index (n-1)//2, ties by depth then seed (routing statistic)."""
    return sorted(rs)[(len(rs) - 1) // 2]


def _ipt(sv, keep=(7, 9)) -> float:
    from qiskit.quantum_info import entropy, partial_trace
    rho = partial_trace(sv, [k for k in range(10) if k not in keep])
    return float(entropy(partial_trace(rho, [0]), base=2) + entropy(partial_trace(rho, [1]), base=2)
                 - entropy(rho, base=2))


def _build(spec: str):
    """spawn-pool worker: one circuits build (the tket block synthesis is the slow part, ~20 s on a box)."""
    import circuits as C
    d = json.loads(spec)
    c = C.Candidate.model_validate(d["cand"])
    return C.step_circuit(c) if d["kind"] == "step" else C.full_circuit(c, with_prep=d["kind"] == "full")


def build_all(specs: dict) -> dict:
    """{key: (Candidate, "full" | "noprep" | "step")} -> {key: QuantumCircuit}, built in parallel."""
    js = [json.dumps({"cand": c.model_dump(mode="json"), "kind": k}) for c, k in specs.values()]
    n = min(len(js), int(os.environ.get("LH_CPUS") or os.cpu_count() or 1))
    if n > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(n, mp_context=mp.get_context("spawn")) as ex:
            res = list(ex.map(_build, js))
    else:
        res = [_build(j) for j in js]
    return dict(zip(specs, res))


def logical_state(qc_noprep, x: ProtocolListing, S: RouteSettings):
    """|Bell>_PQ |VQA(theta)> evolved by the protocol circuit (no prep) built by circuits from the listing, SYNTHESISED
    with the pinned plugin first (an unsynthesised PauliEvolutionGate simulates as the exact exponential)."""
    from qiskit import transpile
    from qiskit.quantum_info import Statevector
    from qiskit.transpiler.passes import HLSConfig
    qc = transpile(qc_noprep, basis_gates=SYNTH_BASIS, optimization_level=0,
                   hls_config=None if S.hls == "default" else HLSConfig(PauliEvolution=[S.hls]))
    init = np.kron(np.array([1, 0, 0, 1]) / np.sqrt(2), vqa.state(np.asarray(x.vqa_theta, float), x.jw))
    return Statevector(init).evolve(qc)


def shipped_checks(L: Listing, P: str, x: ProtocolListing, S: RouteSettings, be, ship_count: int | None,
                   circ: dict) -> list[Metric]:
    """recount, ISA validity on the backend, state equivalence to the claimed protocol, physics tie to sykg."""
    import circuits as C
    from qiskit import QuantumCircuit, qasm3
    from qiskit.quantum_info import Statevector
    out = []
    try:
        q = qasm3.loads(x.ship.qasm3)
    except Exception as e:                                   # noqa: BLE001
        return [M(f"{P}_ship_parse", "cost", str(e)[:200], "QASM 3 parses", False)]
    n2 = C.count2q(q)
    out.append(M(f"{P}_ship_2q", "cost", float(n2), f"= claimed {x.ship.claimed_2q}", n2 == x.ship.claimed_2q))
    out.append(M(f"{P}_ship_2q_reproduced", "cost", float(n2), f"= pinned transpile of the FITTED circuit at seed "
                 f"{x.ship.seed}: {ship_count}", n2 == ship_count))
    from qiskit.circuit.library import get_standard_gate_name_mapping
    std, tg, bad = get_standard_gate_name_mapping(), be.target, []
    for i in q.data:                                         # a QASM `gate cz ...` redefinition is not the native cz
        nm, qa = i.operation.name, tuple(q.find_bit(b).index for b in i.qubits)
        if nm != "barrier" and (nm not in UNITARY_OPS or type(i.operation) is not type(std[nm])
                                or not tg.instruction_supported(nm, qa)):
            bad.append((nm, qa))
    out.append(M(f"{P}_ship_isa", "cost", str(bad[:3]), f"unitary ISA ops on {BACKEND} couplers", not bad))
    if bad:
        return out
    act = sorted({q.find_bit(b).index for i in q.data for b in i.qubits})
    lay = x.ship.final_layout
    order = list(lay) + [a for a in act if a not in lay]
    if len(lay) != 10 or len(set(lay)) != 10 or len(order) > MAX_ACTIVE or any(not 0 <= k < q.num_qubits for k in lay):
        return out + [M(f"{P}_ship_equiv", "cost", len(order), f"final_layout 10 distinct, <= {MAX_ACTIVE} active", False)]
    pos = {k: i for i, k in enumerate(order)}
    qc = QuantumCircuit(len(order))
    qc.global_phase = q.global_phase
    for i in q.data:
        if i.operation.name not in ("barrier", "delay", "id"):
            qc.append(i.operation, [pos[q.find_bit(b).index] for b in i.qubits])
    sv = Statevector(qc).data
    lg = logical_state(circ[(P, "np", x.mu)], x, S)
    fid = float(abs(np.vdot(lg.data, sv[:1024])) ** 2)
    out.append(gate(f"{P}_ship_equiv", "cost", 1.0 - fid, None, EQUIV_TOL))
    v10 = sv[:1024] / max(np.linalg.norm(sv[:1024]), 1e-300)  # the ROUTED circuit's own state on its final layout
    p, tr = PROTO[P][0], PROTO[P][3]
    m = two_sided(L.terms, L.J, p, x.jw, order=x.order, tfd_state=vqa.state(np.asarray(x.vqa_theta, float), x.jw))
    tread = x.t1 if P == "FIX" else x.t0
    err = max(abs(_ipt(lg if mu == x.mu else logical_state(circ[(P, "np", mu)], x, S))
                  - float(m.ipt(mu, [tread], tr, x.t0)[0])) for mu in (-sykg.MU, sykg.MU))
    out.append(gate(f"{P}_circuit_physics_tie", "cost", err, None, IPT_TOL))
    from qiskit.quantum_info import Statevector as SV
    out.append(gate(f"{P}_ship_ipt_tie", "cost", abs(_ipt(SV(v10)) - float(m.ipt(x.mu, [tread], tr, x.t0)[0])), None,
                    SHIP_IPT_TOL))                                # I_PT of the routed shipped state = sykg (mu = listing mu)
    return out


def check_cost(L: Listing, B: Bars) -> list[Metric]:
    import circuits as C
    m2e, K, out = B.files["m2e"], len(L.terms), []
    pin = get(m2e, "settings")
    if pin is None:
        return [M("m2e_bars", "cost", None, "v2_bars/m2e.json missing", False)]
    S = RouteSettings.model_validate(pin)
    code_ok = (S.circuits_py_blob == git_blob(HERE / "circuits.py") and S.vqa_py_blob == git_blob(HERE / "vqa.py")
               and (S.synth != "tket" or f"pytket=={VERSIONS['pytket']}" in S.tket)
               and S.jw == {P: PROTO[P][1] for P in PROTO})
    if S.approximation_degree != 1.0 or not S.backend.startswith("FakeNighthawk") or not code_ok:
        return [M("m2e_settings", "cost", S.model_dump_json(), "approximation_degree 1, FakeNighthawk, bars routed with "
                  "THIS circuits.py / vqa.py (git blobs), the pinned pytket and the v3 map (S2 both)", False)]
    be = C.backend(BACKEND)
    route = lambda qc, s: C.route(qc, be, s, S.opt_level, S.layout_method, S.routing_method, S.hls)
    specs = {}
    for P, x in L.protocols.items():
        if x.jw == PROTO[P][1]:
            specs[(P, "gen")] = (_cand(L, P, x, S), "full")
            specs[(P, "fit")] = (_cand(L, P, x, S, fitted=True), "full")
            for mu in {x.mu, -sykg.MU, sykg.MU}:
                specs[(P, "np", mu)] = (_cand(L, P, x, S, mu), "noprep")
    if "FIX" in L.protocols:
        specs["step"] = (_cand(L, "FIX", L.protocols["FIX"], S), "step")
    circ = build_all(specs)
    for P, x in L.protocols.items():
        if x.jw != PROTO[P][1]:
            out.append(M(f"{P}_cost", "cost", None, "map not supported by this version", False))
            continue
        qc = circ[(P, "gen")]
        bar = get(m2e, "typical_med_of_med", P, str(K))
        bar = None if bar is None else COST_FACTOR * bar
        basis, npass = get(m2e, "basis", P, str(K)), get(m2e, "ad_pass_n", P, str(K))
        out.append(M(f"{P}_bar_basis", "cost", f"{basis} (n_pass_a_d {npass})", "typical samples passing "
                     "a.-d. (>= 50), else all with the pass count recorded",
                     (basis == "pass_a_d" and (npass or 0) >= 50) or (basis == "all" and npass is not None and npass < 50)))
        sa, sb = get(m2e, "seeds_a"), get(m2e, "seeds_b")
        if len(sa) != N_SEEDS or len(sb) != N_SEEDS or set(sa) & set(sb):
            out.append(M(f"{P}_seed_lists", "cost", None, "pinned 20 + 20 disjoint", False))
            continue
        for tag, seeds in (("A", sa), ("B", sb)):
            rs = [(lambda r: (r[1], r[2], s))(route(qc, s)) for s in seeds]
            med = lower_median(rs)
            out.append(gate(f"{P}_2q_median_{tag}", "cost", float(med[0]), None, bar))
            out.append(M(f"{P}_2q_min_{tag}", "cost", float(min(rs)[0]), "reported"))
            out.append(M(f"{P}_depth_median_{tag}", "cost", float(med[1]), "reported"))
        ship_count = route(circ[(P, "fit")], x.ship.seed)[1] if x.ship.seed in sa else None
        out += shipped_checks(L, P, x, S, be, ship_count, circ)
        out.append(M(f"{P}_angles_note", "cost", None, "the bar statistic routes generic VQA angles (vqa_seed, as the "
                     "baselines); the ship is the fitted-theta circuit, recounted on its own. 2q measured angle independent "
                     "for generic angles (4 sets x 6 seeds identical, test_angle_independence); fitted Eq. 11 "
                     "SYM differed on 2/20 seeds, so the shipped 2q is reported, never the statistic"))
    if "FIX" in L.protocols:                                  # Table S1-style one Lie-1 step of e^(-i H_L t0), reported
        st = lower_median([(lambda r: (r[1], r[2], s))(route(circ["step"], s)) for s in get(m2e, "seeds_a")])
        out += [M("step_2q_median", "cost", float(st[0]), "reported (Table S1 style)"),
                M("step_depth_median", "cost", float(st[1]), "reported")]
    return out


# ============================ entry ==============================================================================
def verify(listing, bars: Bars | None = None, cost: bool = True, selection: Selection | None = None,
           fit: bool = True) -> Verdict:
    """listing: Listing, dict, JSON text or a path. cost=False skips the cost gate, fit=False skips the verifier's own VQA fit (vqa gate);
    either makes the Verdict fail."""
    me = sha256(Path(__file__))
    try:
        if isinstance(listing, Listing):
            L = listing
        elif isinstance(listing, dict):
            L = Listing.model_validate(listing)
        else:
            s = str(listing)
            L = Listing.model_validate_json(Path(s).read_text() if not s.lstrip().startswith("{") else s)
    except (ValidationError, OSError, ValueError) as e:
        return Verdict(id="?", version=VERSION, verify_sha256=me, passed=False, m2_countable=False,
                       failed=["listing_schema"], metrics=[],
                       error=str(e)[:500])
    B = bars or load_bars()
    ms: list[Metric] = []
    try:
        ms += check_pins(B)
        for f in (check_listing, check_chaos, check_typicality, check_diag):
            ms += f(L, B)
        ms += check_signal_and_eta(L, B, fit)
        ms += check_selection(L, selection)
        ms += check_cost(L, B) if cost else [M("m2e", "cost", None, "not evaluated (cost=False)", False)]
        err = None
    except Exception as e:                                   # noqa: BLE001  any crash = not certified
        err = f"{type(e).__name__}: {e}"[:500]
    failed = [x.name for x in ms if x.passed is False]
    ok = not failed and err is None
    return Verdict(id=L.id, version=VERSION, verify_sha256=me, passed=ok, m2_countable=ok and L.term_set == "uniform",
                   failed=failed, metrics=ms, error=err)


def _one(x: dict) -> str:
    sel = Selection.model_validate(x["sel"]) if x["sel"] else None
    try:
        return verify(x["listing"], selection=sel).model_dump_json()
    except Exception as e:                                   # noqa: BLE001  never lose the other verdicts
        return Verdict(id=str(json.loads(x["listing"]).get("id", "?")), version=VERSION, verify_sha256="", passed=False,
                       m2_countable=False, failed=["crash"], metrics=[], error=str(e)[:500]).model_dump_json()


def load_fits(path: Path) -> None:
    for line in path.read_text().splitlines():
        if line.strip():
            x = FitShard.model_validate_json(line)
            FITS.setdefault(x.key, []).append(x)
