"""`uv run verify`: Eq. 11 known answers + every listing's noiseless numbers, recomputed from its hamiltonian.json.

Known answers: Eq. 11 level-spacing ratios (r_odd 0.537163, r_even 0.478812), sykg.selftest (TFD energy, JW-map
invariance, Trotter limits), and Eq. 11 exact pair-(1,2) Delta I_PT peaks vs Fig. S1 of arXiv:2604.10090 (SYM 0.212,
FIX 0.276, +/- 0.003 bits read off the figure).
Listings (candidates/<id>/ and kit/<id>/): r_odd, r_even, SFF / DENSE, the 8 Delta I_PT peaks and 2 revivals, VQA F and
Trotter eta from the listed angles, OTOC E3, the noiseless dI of the shipped circuit model at its readout time and the
two-qubit count of every shipped QASM file, each against the number the run reported (row.json or predictions.json)."""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from pydantic import BaseModel

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import diag  # noqa: E402
import sykg  # noqa: E402
import verify as V  # noqa: E402
import vqa  # noqa: E402

TOL = 1e-8                   # recomputed vs reported, same code and pins (observed ~1e-15)
TOL_I = 0.003                #  Eq. 11 vs Fig. S1 of arXiv:2604.10090 (bits, figure read-off)
TS = np.round(np.arange(0, 20.0001, 0.1), 1)
FIG_S1 = {"sym": 0.212, "fix": 0.276}   # Eq. 11 exact pair-(1, 2) Delta I_PT peak, arXiv:2604.10090 Fig. S1


class Check(BaseModel):
    group: str               # eq11 | <listing id>
    name: str
    value: float | None
    ref: float | None
    tol: float
    ok: bool


class Report(BaseModel):
    known_ok: int            # 1 iff every check passes
    eq11_r_err: float
    n_checks: int
    n_failed: int
    max_err: float           # largest |value - ref| over the float listing checks (tolerance TOL or SHIP_IPT_TOL)
    failed: list[str]
    seconds: float
    checks: list[Check]


def chk(group: str, name: str, value, ref, tol: float = TOL) -> Check:
    v = None if value is None else float(value)
    r = None if ref is None else float(ref)
    ok = v is not None and r is not None and np.isfinite(v) and abs(v - r) <= tol
    return Check(group=group, name=name, value=v, ref=r, tol=tol, ok=bool(ok))


# ------------------------------------------------------------------ Eq. 11 known answers
def eq11_checks() -> list[Check]:
    st = sykg.selftest(False)
    out = [chk("eq11", "r_odd", st["r_odd"], sykg.R_ODD_EQ11, sykg.R_TOL),
           chk("eq11", "r_even", st["r_even"], sykg.R_EVEN_EQ11, sykg.R_TOL),
           chk("eq11", "selftest_all_ok", float(st["all_ok"]), 1.0, 0.0)]
    T, J = sykg.eq11()
    for p, ref in FIG_S1.items():
        m = sykg.TwoSided(T, J, p)
        out.append(chk("eq11", f"figS1_{p}_dI", sykg.peak_window(TS, m.ipt(-sykg.MU, TS) - m.ipt(sykg.MU, TS))[0], ref, TOL_I))
    return out


# ------------------------------------------------------------------ listings
def qasm_files(d: Path, P: str) -> tuple[Path, Path]:
    """(mu = -12, mu = +12) shipped QASM of protocol P: kit/ naming P_mu-12.qasm, candidates/ naming P.qasm."""
    m = d / f"{P}_mu-12.qasm"
    return (m if m.exists() else d / f"{P}.qasm"), d / f"{P}_mu+12.qasm"


def count2q_file(p: Path) -> int:
    import circuits as C
    from qiskit import qasm3
    return C.count2q(qasm3.loads(p.read_text()))


def model_dI(terms, J, P: str, x: dict) -> float:
    """noiseless dI = I_PT(mu -12) - I_PT(mu +12) of the protocol's one-step Trotter circuit with the listed VQA TFD,
    at the listed readout time (the shipped circuit's model; verify.py shipped_checks tie)."""
    p, tr = V.PROTO[P][0], V.PROTO[P][3]
    m = V.two_sided(terms, J, p, x["jw"], order=x["order"], tfd_state=vqa.state(np.asarray(x["vqa_theta"], float), x["jw"]))
    t = x["t1"] if P == "FIX" else x["t0"]
    return float(m.ipt(-sykg.MU, [t], tr, x["t0"])[0] - m.ipt(sykg.MU, [t], tr, x["t0"])[0])


def listing_values(d: Path) -> dict:
    """every recomputable noiseless number of one listing directory."""
    h = json.loads((d / "hamiltonian.json").read_text())
    terms, J, pr = [tuple(t) for t in h["terms"]], h["J"], h["protocols"]
    l1 = json.loads((HERE / "v2_bars" / "l1.json").read_text())
    Eo, Ee = sykg.sector_spectra(terms, J)
    v = {"r_odd": float(sykg.r_mean(Eo)), "r_even": float(sykg.r_mean(Ee)),
         "sff_rel": float(sykg.sff_scalar(sykg.spectrum(terms, J))) / l1["dense_sff_scalar_mean"],
         "otoc_E3": float(diag.passes(terms, J, diag.load_thresholds())["metrics"]["E3"])}
    v |= V.signal_values(terms, J, {"fix": pr["FIX"]["order"], "sym": pr["SYM"]["order"]})
    for P, x in pr.items():
        p = P.lower()
        v[f"F_{p}"], v[f"eta_{p}"] = V.eta_value(terms, J, p, x["order"], x["vqa_theta"])
        v[f"dI0_{p}"] = model_dI(terms, J, P, x)
        qm, qp = qasm_files(d, P)
        v[f"2q_{p}_m"] = float(count2q_file(qm)) if qm.exists() else None
        v[f"2q_{p}_p"] = float(count2q_file(qp)) if qp.exists() else None
        v[f"claimed_2q_{p}"] = float(x["ship"]["claimed_2q"])
    return v


def references(d: Path) -> dict:
    """the reported numbers, keyed like listing_values: row.json (candidates/) or predictions.json (kit/)."""
    h = json.loads((d / "hamiltonian.json").read_text())
    ref: dict = {}
    for P, x in h["protocols"].items():
        ref[f"F_{P.lower()}"] = x["vqa_F"]
        ref[f"2q_{P.lower()}_m"] = x["ship"]["claimed_2q"]
    if (d / "row.json").exists():
        r = json.loads((d / "row.json").read_text())
        ref |= {k: r[k] for k in ("r_odd", "r_even", "sff_rel", "eta_fix", "eta_sym")}
        ref |= {k: r["dI"][k] for k in V.SIG8} | {"sym_rev": r["rev_sym"], "fix_rev": r["rev_fix"]}
        ref |= {"otoc_E3": r["extra"].get("otoc_E3")} | {f"2q_{P.lower()}_m": r["ship"][P]["2q"] for P in r["ship"]}
    if (d / "verdict.json").exists():          # a listing whose row.json is a placeholder: its frozen-verifier verdict metrics
        vm = json.loads((d / "verdict.json").read_text())["metrics"]
        names = {"sff_scalar": "sff_rel", "fix_vqa_F": "F_fix", "sym_vqa_F": "F_sym", "fix_eta": "eta_fix",
                 "sym_eta": "eta_sym"}
        ref |= {names.get(k, k): x for k, x in vm.items() if ref.get(names.get(k, k)) is None}
    if (d / "predictions.json").exists():
        s = json.loads((d / "predictions.json").read_text())
        st = s.get("stats", {})
        ref |= {k2: st[k1] for k1, k2 in (("r_odd", "r_odd"), ("r_even", "r_even"), ("sff_rel_DENSE", "sff_rel"),
                                          ("otoc_E3", "otoc_E3"), ("eta_FIX", "eta_fix"), ("eta_SYM", "eta_sym"))
                if k1 in st}
        for P, x in s["protos"].items():
            p = P.lower()
            ref |= {f"dI0_{p}": x["dI0"], f"{p}_ex_12": x["dI_exact_peak"], f"{p}_tr_12": x["dI_trotter_peak"],
                    f"2q_{p}_m": x["twoq"][0], f"2q_{p}_p": x["twoq"][1]}
    return {k: v for k, v in ref.items() if v is not None}


def listing_checks(d: str) -> list[dict]:
    """spawn-pool worker: one listing directory -> its checks (as dicts)."""
    d = Path(d)
    gid = f"{d.parent.name}/{d.name}"
    v, ref = listing_values(d), references(d)
    out = []
    for k in sorted(ref):
        tol = (0.05 if k == "otoc_E3" else 0.0 if k.startswith("2q_") else V.SHIP_IPT_TOL if k.startswith("dI0_")
               else TOL)   # E3 is reported to 0.01; the reported dI0 is the routed circuits' noiseless simulation
        # (median over five routing seeds), which carries the routed ISA circuit's resynthesis error (<= ~1e-6)
        out.append(chk(gid, k, v.get(k), ref[k], tol).model_dump())
    for p in ("fix", "sym"):                 # every shipped file carries the claimed count
        if v.get(f"2q_{p}_m") is not None:
            out.append(chk(gid, f"2q_{p}_m_eq_claimed", v[f"2q_{p}_m"], v[f"claimed_2q_{p}"], 0.0).model_dump())
    return out


def listing_dirs(root: Path) -> list[Path]:
    return sorted(p.parent for base in ("candidates", "kit") for p in (root / base).glob("*/hamiltonian.json"))


def run(root: Path, only: list[str] | None = None, workers: int | None = None) -> Report:
    t0 = time.time()
    dirs = [d for d in listing_dirs(root) if not only or any(o in str(d) for o in only)]
    n = max(1, min(len(dirs), workers or os.cpu_count() or 1))
    with ProcessPoolExecutor(n, mp_context=mp.get_context("spawn")) as ex:
        fut = ex.map(listing_checks, [str(d) for d in dirs])
        eq = eq11_checks()
        cs = eq + [Check.model_validate(c) for r in fut for c in r]
    fail = [f"{c.group}:{c.name}" for c in cs if not c.ok]
    errs = [abs(c.value - c.ref) for c in cs if c.tol in (TOL, V.SHIP_IPT_TOL) and c.value is not None and c.ref is not None]
    return Report(known_ok=int(not fail), eq11_r_err=float(sykg.eq11_check()["eq11_r_err"]), n_checks=len(cs),
                  n_failed=len(fail), max_err=max(errs, default=0.0), failed=fail, seconds=time.time() - t0, checks=cs)


def default_root() -> Path:
    for p in (Path.cwd(), HERE.parents[1]):
        if (p / "candidates").is_dir():
            return p
    return HERE.parents[1]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="verify", description=__doc__.split("\n")[0])
    ap.add_argument("--only", nargs="*", help="substring(s) of listing directories to check (default: all)")
    ap.add_argument("--out", default="out/verify.json", help="report JSON path")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--root", default=None, help="repository root (default: cwd or the checkout)")
    a = ap.parse_args(argv)
    rep = run(Path(a.root) if a.root else default_root(), a.only, a.workers)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(rep.model_dump_json(indent=1))
    groups = sorted({c.group for c in rep.checks}, key=lambda g: (g != "eq11", g))
    for g in groups:
        cs = [c for c in rep.checks if c.group == g]
        bad = [c.name for c in cs if not c.ok]
        e = max((abs(c.value - c.ref) for c in cs if c.value is not None and c.ref is not None and not c.name.endswith(
                 "_t") and not c.name.startswith(("2q_", "otoc", "selftest"))),
                default=0.0)
        print(f"{'ok  ' if not bad else 'FAIL'} {g:58s} {len(cs):3d} checks  max|diff| {e:.1e}"
              + (f"  failed: {', '.join(bad)}" if bad else ""))
    print(f"known_ok {rep.known_ok}  checks {rep.n_checks}  failed {rep.n_failed}  max|diff| {rep.max_err:.2e}  "
          f"eq11_r_err {rep.eq11_r_err:.2e}  {rep.seconds:.0f} s  -> {a.out}")
    return 0 if rep.known_ok else 1


if __name__ == "__main__":
    sys.exit(main())
