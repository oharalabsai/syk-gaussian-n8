"""`uv run noise <listing dir>`: FakeNighthawk noisy pair-(1,2) dI = I_PT(mu -12) - I_PT(mu +12) of one listing's
shipped circuit pair (one protocol, the shipped routing seed), as in the run's noise simulation.

Each routed circuit is restricted to its active qubits (each keeps its own FakeNighthawk errors) and simulated as a
density matrix (qiskit-aer): `ideal` (no noise), `gate` (NoiseModel gate errors incl. gate-time T1 / T2 relaxation),
`gate_ro` (+ readout as a confusion map on naive two-qubit Pauli tomography of P and register qubit 1, PSD clipped;
the run's PRIMARY variant). Known answer: the ideal dI equals the sykg model of the same circuit (1e-5). Reference: the
kit's predictions.json `<P>_shipped_seed_gate_ro` when present."""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
from pydantic import BaseModel

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import circuits as C  # noqa: E402
from reproduce import default_root, model_dI, qasm_files  # noqa: E402
from route import resolve  # noqa: E402

KEEP = (7, 9)                        # logical register qubit 1 (replaces T) and P
MODEL_TOL, REF_TOL = 1e-5, 1e-6
_BE: dict = {}


class Noisy(BaseModel):
    id: str
    P: str
    seed: int
    twoq: list[int]                  # mu -12, mu +12
    n_active: list[int]
    dI_model: float                  # sykg, noiseless, listed theta
    dI_ideal: float
    dI_gate: float
    dI_gate_ro: float
    ref_gate_ro: float | None
    model_ok: bool
    ref_ok: bool | None


def be():
    if "nh" not in _BE:
        _BE["nh"] = C.backend("nighthawk")
    return _BE["nh"]


# ------------------------------------------------------------------ simulation
def compact(tq, final_layout):
    """active qubits only: final layout of logical 0..9 first, then the other active qubits (verify.py order)."""
    from qiskit import QuantumCircuit
    lay = list(final_layout)
    act = sorted({tq.find_bit(b).index for i in tq.data if i.operation.name not in ("barrier", "delay")
                  for b in i.qubits})
    order = lay + [a for a in act if a not in lay]
    pos = {k: i for i, k in enumerate(order)}
    qc = QuantumCircuit(len(order))
    qc.global_phase = tq.global_phase
    for i in tq.data:
        if i.operation.name == "barrier":
            continue
        qs = [tq.find_bit(b).index for b in i.qubits]
        if all(q in pos for q in qs):
            qc.append(i.operation, [pos[q] for q in qs])
    return qc, order


def noise_model(order):
    """FakeNighthawk gate errors (depolarizing + thermal relaxation on gate durations) on the compact register."""
    from qiskit_aer.noise import NoiseModel
    from qiskit_aer.noise.device.models import basic_device_gate_errors
    if "ge" not in _BE:
        _BE["ge"] = basic_device_gate_errors(target=be().target)
    pos = {k: i for i, k in enumerate(order)}
    nm = NoiseModel(basis_gates=["cz", "id", "rz", "sx", "x"])
    for name, qs, err in _BE["ge"]:
        if name != "reset" and all(q in pos for q in qs):
            nm.add_quantum_error(err, name, [pos[q] for q in qs])
    return nm


def readout(order, keep) -> list[tuple[float, float]]:
    """(e0 = P(1|0), e1 = P(0|1)) of the physical qubits behind compact indices keep."""
    from qiskit_aer.noise.device.models import basic_device_readout_errors
    ro = {qs[0]: np.asarray(e.probabilities) for qs, e in basic_device_readout_errors(target=be().target)}
    return [(float(ro[order[k]][0, 1]), float(ro[order[k]][1, 0])) for k in keep]


_P = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])]


def ro_map(rho: np.ndarray, e: list[tuple[float, float]]) -> np.ndarray:
    """naive 2-qubit Pauli tomography under independent readout confusion: per qubit, c_s -> a c_s + b c_I (s = X,Y,Z),
    a = 1 - e0 - e1, b = e1 - e0; rho index little endian (qubit 0 = keep[0]); then clipped to PSD, trace 1."""
    ops = {(i, j): np.kron(_P[j], _P[i]) for i in range(4) for j in range(4)}
    c = np.array([[np.trace(rho @ ops[(i, j)]).real for j in range(4)] for i in range(4)])
    (a0, b0), (a1, b1) = [(1 - e0 - e1, e1 - e0) for e0, e1 in e]
    c[1:, :] = a0 * c[1:, :] + b0 * c[0:1, :]
    c[:, 1:] = a1 * c[:, 1:] + b1 * c[:, 0:1]
    r = sum(c[i, j] * ops[(i, j)] for i in range(4) for j in range(4)) / 4
    w, v = np.linalg.eigh((r + r.conj().T) / 2)
    w = np.clip(w, 0, None)
    return (v * (w / w.sum())) @ v.conj().T


def _S(r) -> float:
    w = np.clip(np.linalg.eigvalsh(r), 0, None)
    w = w[w > 1e-15]
    return float(-(w * np.log2(w)).sum())


def ipt2(rho: np.ndarray) -> float:
    r = rho.reshape(2, 2, 2, 2)                      # (q1, q0, q1', q0')
    return _S(np.einsum("abcb->ac", r)) + _S(np.einsum("abad->bd", r)) - _S(rho)


def dm(qc, keep, nm=None) -> np.ndarray:
    from qiskit_aer import AerSimulator
    q = qc.copy()
    q.save_density_matrix(qubits=list(keep), label="rho")
    sim = AerSimulator(method="density_matrix", noise_model=nm)
    return np.asarray(sim.run(q, shots=1).result().data(0)["rho"])


def simulate(tq, final_layout) -> dict:
    """one routed circuit -> I_PT of logical (7, 9): ideal, gate, gate_ro."""
    qc, order = compact(tq, final_layout)
    nm = noise_model(order)
    rg = dm(qc, KEEP, nm)
    return {"n_active": len(order), "ideal": ipt2(dm(qc, KEEP)), "gate": ipt2(rg),
            "gate_ro": ipt2(ro_map(rg, readout(order, KEEP)))}


# ------------------------------------------------------------------ shipped files
def header_layout(text: str) -> list[int]:
    """final layout (physical qubit of logical 0..9) from the QASM header line 'logical->physical at the end: P->$n,
    Q->$n, reg1->$n, ..., reg8->$n' (P = logical 9, Q = 8, register qubit k = logical 8 - k)."""
    m = dict(re.findall(r"(P|Q|reg\d)->\$(\d+)", text.split("include", 1)[0]))
    lay = [0] * 10
    lay[9], lay[8] = int(m["P"]), int(m["Q"])
    for k in range(1, 9):
        lay[C.rq(k)] = int(m[f"reg{k}"])
    return lay


def noisy_listing(d: Path, P: str) -> Noisy:
    from qiskit import qasm3
    ham = json.loads((d / "hamiltonian.json").read_text())
    x = ham["protocols"][P]
    sims = []
    for f in qasm_files(d, P):
        t = f.read_text()
        sims.append(simulate(qasm3.loads(t), header_layout(t)) | {"2q": C.count2q(qasm3.loads(t))})
    m, p = sims
    pred = d / "predictions.json"
    ref = json.loads(pred.read_text()).get("stats", {}).get(f"{P}_shipped_seed_gate_ro") if pred.exists() else None
    dmod = model_dI([tuple(q) for q in ham["terms"]], ham["J"], P, x)
    gro = m["gate_ro"] - p["gate_ro"]
    return Noisy(id=ham["id"], P=P, seed=x["ship"]["seed"], twoq=[m["2q"], p["2q"]],
                 n_active=[m["n_active"], p["n_active"]], dI_model=dmod, dI_ideal=m["ideal"] - p["ideal"],
                 dI_gate=m["gate"] - p["gate"], dI_gate_ro=gro, ref_gate_ro=ref,
                 model_ok=abs(m["ideal"] - p["ideal"] - dmod) <= MODEL_TOL,
                 ref_ok=None if ref is None else abs(gro - ref) <= REF_TOL)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="noise", description=__doc__.split("\n")[0])
    ap.add_argument("listing", nargs="?", default="p325_tb_23_24_44300310_1429_19277_201",
                    help="listing directory, its name, or a unique substring (default: the kit's primary listing)")
    ap.add_argument("--protocol", choices=("FIX", "SYM"), default="FIX")
    ap.add_argument("--out", default="out/noise.json")
    a = ap.parse_args(argv)
    t0 = time.time()
    r = noisy_listing(resolve(a.listing, default_root()), a.protocol)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(r.model_dump_json(indent=1))
    print(f"{r.id} {r.P} seed {r.seed} 2q {r.twoq} active {r.n_active}: dI model {r.dI_model:.6f} ideal "
          f"{r.dI_ideal:.6f} gate {r.dI_gate:.6f} gate_ro {r.dI_gate_ro:.6f} (kit {r.ref_gate_ro}) "
          f"model_ok {r.model_ok} ref_ok {r.ref_ok}  {time.time() - t0:.0f} s -> {a.out}")
    return 0 if r.model_ok and r.ref_ok is not False else 1


if __name__ == "__main__":
    sys.exit(main())
