"""`uv run route <listing dir>`: rebuild the fitted circuit of a listing from its hamiltonian.json (tket block synthesis,
S2 map, mu = -12 and +12), route it on FakeNighthawk at the shipped transpiler seed with the pinned setting
(optimization level 3, approximation degree 1, default layout / routing / HLS), and compare the two-qubit count and the
QASM 3 body with the shipped files. Routing is platform sensitive (the run's numbers come from x86_64 Linux; an arm64 Mac
gave 199 instead of 200 on one FIX circuit); the physics is not."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from pydantic import BaseModel

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import circuits as C  # noqa: E402
from reproduce import default_root, qasm_files  # noqa: E402

SETTINGS = dict(opt_level=3, layout_method=None, routing_method=None, hls="default")   # v2_bars/m2e.json settings


class Routed(BaseModel):
    id: str
    P: str
    mu: float
    seed: int
    twoq: int
    oneq: int
    depth: int
    shipped_2q: int | None
    match_2q: bool | None
    match_qasm: bool | None        # routed QASM 3 = shipped file, header comments stripped
    final_layout: list[int]


def candidate(ham: dict, P: str, mu: float | None = None) -> C.Candidate:
    """the shipped (fitted) circuit's Candidate (the run's candidate builder)."""
    x = ham["protocols"][P]
    return C.Candidate(id=ham["id"], terms=[tuple(q) for q in ham["terms"]], J=ham["J"], jw=x["jw"], protocol=P,
                       trotter=x["trotter"], t0=x["t0"], t1=x["t1"], mu=x["mu"] if mu is None else mu, order=x["order"],
                       backward="same_order", prep="I", vqa_seed=0, synth="tket", theta=list(x["vqa_theta"]))


def body(q: str) -> str:
    return "\n".join(ln for ln in q.strip().splitlines() if not ln.startswith("//"))


def route_listing(d: Path, protocols=("FIX", "SYM"), seed: int | None = None) -> list[Routed]:
    from qiskit import qasm3
    ham = json.loads((d / "hamiltonian.json").read_text())
    be, out = C.backend("nighthawk"), []
    for P in protocols:
        x = ham["protocols"][P]
        s = x["ship"]["seed"] if seed is None else seed
        for mu, f in zip((-12.0, 12.0), qasm_files(d, P)):
            tq, n2, dep = C.route(C.full_circuit(candidate(ham, P, mu)), be, s, **SETTINGS)
            ops = dict(tq.count_ops())
            ship = f.read_text() if f.exists() and seed is None else None
            n_ship = C.count2q(qasm3.loads(ship)) if ship else None
            out.append(Routed(id=ham["id"], P=P, mu=mu, seed=s, twoq=n2, depth=dep,
                              oneq=sum(v for k, v in ops.items() if k in ("rz", "sx", "x", "id")),
                              shipped_2q=n_ship, match_2q=None if ship is None else n2 == n_ship,
                              match_qasm=None if ship is None else body(C.to_qasm3(tq)) == body(ship),
                              final_layout=[int(p) for p in tq.layout.final_index_layout()]))
    return out


def resolve(arg: str, root: Path) -> Path:
    for p in (Path(arg), root / "kit" / arg, root / "candidates" / arg):
        if (p / "hamiltonian.json").exists():
            return p
    hits = [p.parent for p in root.glob("*/*/hamiltonian.json") if arg in p.parent.name]
    kit = [h for h in hits if h.parent.name == "kit"]
    if len(hits) == 1 or len(kit) == 1:            # a label in both kit/ and candidates/ resolves to kit/
        return (kit or hits)[0]
    raise SystemExit(f"no unique listing matches {arg!r}: {[str(h) for h in hits]}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="route", description=__doc__.split("\n")[0])
    ap.add_argument("listing", nargs="?", default="p325_tb_23_24_44300310_1429_19277_201",
                    help="listing directory, its name, or a unique substring (default: the kit's primary listing)")
    ap.add_argument("--protocol", choices=("FIX", "SYM", "both"), default="both")
    ap.add_argument("--seed", type=int, default=None, help="transpiler seed (default: the shipped seed; then compare)")
    ap.add_argument("--out", default="out/route.json")
    a = ap.parse_args(argv)
    t0 = time.time()
    rs = route_listing(resolve(a.listing, default_root()), ("FIX", "SYM") if a.protocol == "both" else (a.protocol,),
                       a.seed)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps([r.model_dump() for r in rs], indent=1))
    for r in rs:
        print(f"{r.id} {r.P} mu {r.mu:+.0f} seed {r.seed}: 2q {r.twoq} (shipped {r.shipped_2q}, match {r.match_2q}), "
              f"1q {r.oneq}, depth {r.depth}, QASM identical {r.match_qasm}")
    ok = all(r.match_2q is not False for r in rs)
    print(f"route_ok {int(ok)}  {time.time() - t0:.0f} s  -> {a.out}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
