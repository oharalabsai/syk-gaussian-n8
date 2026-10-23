"""GIST traversable-wormhole circuits for N = 8 sparse SYK, and their routing on IBM fake backends.

Conventions (TW = arXiv:2604.10090 .txt line numbers):
- H_{L,R} = sum_T c_T psi_a psi_b psi_c psi_d (ordered, i^{q/2} absorbed in c_T, as Eq. 11), {psi_a^i, psi_b^j} = delta.
- JW maps (TW:621-640, read off the rendered PDF eqs. S2/S3): psi_R^1 = X/sqrt2, psi_R^2 = Y/sqrt2 on register qubit 1;
  psi_L^1 = Z X, psi_L^2 = Z Y (X/Y on register qubit 2); pair j >= 3: psi_L^j = Z^p X, psi_R^j = Z^p Y, i.e. X/Y on
  register qubit p + 1. S2 (FIX) p = 5,6,4,2,7,3 for j = 3..8; S3 (SYM) p = 6,7,3,4,2,5.
- Register qubit k (1..8) is qiskit qubit 8 - k, Q (message) = 8, P (reference) = 9, so a Pauli label written in
  register order 1..8 is a qiskit label on qubits range(8) and Operator(qc) is kron(reg1, ..., reg8).
- Full circuit (TW:286-361): Bell_PQ, |I> (X, H, S + 1 CNOT, TW:1184-1190), VQA (6 Rx Rz layers, 5 linear
  CNOT layers), e^{+i H_L t0}, fermionic SWAP_QL (eq. S6) with its JW string, e^{-i H_L t0}, e^{i mu V}, e^{-i H_R t1}.
  Readout: measure register qubit 1 and P (T replaced, TW:296-300): 1q only, not built.
- Trotter: Lie-1 for FIX, Suzuki-2 for SYM (TW:320-323). Backward = same term order, negated time.
"""
from __future__ import annotations

import warnings
from functools import lru_cache
from typing import Literal

import numpy as np
from pydantic import BaseModel

warnings.filterwarnings("ignore", message=".*not intended to represent typical.*")

from qiskit import QuantumCircuit, qasm3, transpile  # noqa: E402
from qiskit.circuit.library import PauliEvolutionGate  # noqa: E402
from qiskit.synthesis import LieTrotter, SuzukiTrotter  # noqa: E402
from qiskit.quantum_info import Pauli, SparsePauliOp  # noqa: E402

import os, sys  # noqa: E401,E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sykg  # noqa: E402  (sykg owns the JW maps; imported, never copied)

JW_P = {"S2": tuple(sykg.JW_EXP["fix"]), "S3": tuple(sykg.JW_EXP["sym"])}   # p of psi^j = Z^p X, j = 3..8 (TW S2, S3)
assert JW_P == {"S2": (5, 6, 4, 2, 7, 3), "S3": (6, 7, 3, 4, 2, 5)}, "sykg JW maps changed: re-read TW S2/S3"
Q, PQ = 8, 9                                                   # qiskit indices of the message and reference qubits
N_VQA_LAYERS, N_VQA_CHAINS = 6, 5                              # TW:1184-1190
Q_SYK, N_SYK = 4, 8
EQ11_TERMS = [(1, 2, 4, 7), (1, 2, 6, 8), (1, 3, 4, 5), (1, 3, 7, 8), (1, 4, 5, 6),
              (1, 5, 6, 7), (2, 3, 4, 8), (2, 3, 5, 7), (2, 4, 6, 7), (3, 5, 6, 8)]   # TW:236-246
EQ11_SIGNS = [-1, 1, -1, 1, -1, 1, -1, 1, -1, 1]


def rq(k: int) -> int:
    """register qubit k (1..8) -> qiskit index."""
    return 8 - k


class Candidate(BaseModel):
    id: str
    terms: list[tuple[int, int, int, int]]
    J: list[float]                                 # c_T per term (Eq. 11: sign / sqrt 5)
    jw: str | list[int] = "S2"                     # "S2", "S3" or the host register qubit of pair j = 3..8
    protocol: Literal["FIX", "SYM"] = "FIX"
    trotter: Literal["lie1", "suzuki2"] | None = None   # None: lie1 for FIX, suzuki2 for SYM
    t0: float = 1.8
    t1: float = 2.3                                # FIX readout time; SYM uses t1 = t0
    mu: float = -12.0
    order: list[int] | None = None                 # Trotter term order (indices into terms); None = sykg default
    backward: Literal["same_order", "adjoint"] = "same_order"
    prep: Literal["I", "H"] = "I"                  # |I> start (TW) or an |+>^8 start (variant)
    vqa_seed: int = 0                              # generic VQA angles (counts do not depend on them)
    vqa_order: str | None = None                   # VQA rotation order / CNOT chain; None = vqa.ORDER / vqa.CHAIN
    vqa_chain: str | None = None                   # (e.g. "zx", "down": an alternative block, variant)
    theta: list[float] | None = None               # fitted VQA angles in vqa.py layout (6, 8, 2): [..., 0] thx, [..., 1] thz
    synth: Literal["term", "whole", "whole_nopo", "tket"] = "term"   # per factor, per block (rustiq/LieTrotter),
                                                   # or "tket": per-block GreedyPauliSimp (the shipped setting)
                                                   # (whole_nopo lets qiskit reorder terms: a different Trotter product)

    def trot(self) -> str:
        return self.trotter or ("lie1" if self.protocol == "FIX" else "suzuki2")

    def tread(self) -> float:
        return self.t0 if self.protocol == "SYM" else self.t1

    def seq(self) -> list[int]:
        """Trotter term order (first applied first): explicit order, else sykg.default_term_order (Eq. 11 ->
        EQ11_PDF_ORDER, any other set -> listed order), so circuits and sykg build the same product by default."""
        return list(self.order) if self.order is not None else list(sykg.default_term_order(self.terms))


def eq11(**kw) -> Candidate:
    return Candidate(id=kw.pop("id", "eq11"), terms=EQ11_TERMS, J=[s / np.sqrt(5) for s in EQ11_SIGNS], **kw)


# ---------------------------------------------------------------- Majorana -> Pauli
def hosts(jw: str | list[int]) -> tuple[int, ...]:
    """host register qubit of pair j = 3..8."""
    h = tuple(p + 1 for p in JW_P[jw]) if isinstance(jw, str) else tuple(int(k) for k in jw)
    if sorted(h) != [3, 4, 5, 6, 7, 8]:
        raise ValueError(f"JW hosts must permute 3..8: {h}")
    return h


def majorana_labels(jw: str | list[int]) -> dict[str, dict[int, str]]:
    """Pauli label (register order 1..8) of sqrt2 * psi_{L,R}^j. jw = "S4": one-sided map (TW:645-654) on 4 qubits, L only
    (Table S1 step-count hypothesis; not a two-sided map)."""
    s4 = {j: "Z" * ((j - 1) // 2) + ("X" if j % 2 else "Y") + "I" * (3 - (j - 1) // 2) for j in range(1, 9)}
    if jw == "S4":
        return {"L": s4}
    if jw == "split":      # variant (not GIST's map): R on register 1-4 as S4, L = Z^4 x S4 on register 5-8
        return {"R": {j: s4[j] + "IIII" for j in s4}, "L": {j: "ZZZZ" + s4[j] for j in s4}}
    lab = {"L": {1: "ZX" + "I" * 6, 2: "ZY" + "I" * 6}, "R": {1: "X" + "I" * 7, 2: "Y" + "I" * 7}}
    for j, k in zip(range(3, 9), hosts(jw)):
        lab["L"][j] = "Z" * (k - 1) + "X" + "I" * (8 - k)
        lab["R"][j] = "Z" * (k - 1) + "Y" + "I" * (8 - k)
    return lab


def _product(labels: list[str], coeff: complex) -> tuple[str, complex]:
    op = SparsePauliOp(labels[0], coeffs=[coeff])
    for la in labels[1:]:
        op = op.dot(SparsePauliOp(la))
    op = op.simplify()
    assert len(op.paulis) == 1
    return op.paulis[0].to_label(), complex(op.coeffs[0])


def term_pauli(lab: dict, side: str, term, c: float) -> tuple[str, float]:
    """c psi_a psi_b psi_c psi_d -> (label, real coeff); each Majorana carries 1/sqrt2."""
    la, co = _product([lab[side][j] for j in term], c / 4.0)
    if abs(co.imag) > 1e-12:
        raise RuntimeError("quartic term not Hermitian: convention bug")
    return la, co.real


def h_paulis(c: Candidate, side: str) -> list[tuple[str, float]]:
    lab = majorana_labels(c.jw)
    return [term_pauli(lab, side, t, J) for t, J in zip(c.terms, c.J)]


def v_paulis(jw) -> list[tuple[str, float]]:
    """V = H_int / (q N), H_int = i sum_j psi_L^j psi_R^j (Eq. 3, TW:138-143) -> [(label, coeff)]."""
    lab = majorana_labels(jw)
    out = []
    for j in range(1, 9):
        la, co = _product([lab["L"][j], lab["R"][j]], 1j / 2.0 / (Q_SYK * N_SYK))
        if abs(co.imag) > 1e-12:
            raise RuntimeError("coupling not Hermitian: convention bug")
        out.append((la, co.real))
    return out


# ---------------------------------------------------------------- circuit blocks
def _evo(qc: QuantumCircuit, label: str, angle: float) -> None:
    """append exp(-i angle P) on the register (qubits 0..len(label)-1)."""
    qc.append(PauliEvolutionGate(Pauli(label), time=angle), range(len(label)))


def trotter_seq(c: Candidate, t: float) -> list[tuple[int, float]]:
    """one step of e^{-i H t}: [(term index, time)] in circuit order (first applied first)."""
    o = c.seq()
    if c.trot() == "lie1":
        return [(k, t) for k in o]
    return [(k, t / 2) for k in o[:-1]] + [(o[-1], t)] + [(k, t / 2) for k in reversed(o[:-1])]   # Suzuki-2, merged middle


def evolve(qc: QuantumCircuit, c: Candidate, paulis: list[tuple[str, float]], t: float, backward: bool = False) -> None:
    """e^{-i H t} (backward: e^{+i H t}) as one Trotter step on the register."""
    n = len(paulis[0][0])
    if c.synth == "tket":
        qc.compose(tket_block(tuple(block_rotations(c, paulis, t, backward))), qubits=range(n), inplace=True)
        return
    if c.synth != "term":
        o = c.seq()
        if backward and c.backward == "adjoint":
            o = o[::-1]
        op = SparsePauliOp([paulis[k][0] for k in o], coeffs=[paulis[k][1] for k in o])
        po = c.synth == "whole"
        sy = LieTrotter(reps=1, preserve_order=po) if c.trot() == "lie1" else SuzukiTrotter(order=2, reps=1, preserve_order=po)
        qc.append(PauliEvolutionGate(op, time=-t if backward else t, synthesis=sy), range(n))
        return
    sq = trotter_seq(c, t)
    if backward:
        sq = [(k, -tt) for k, tt in (sq if c.backward == "same_order" else reversed(sq))]
    for k, tt in sq:
        la, co = paulis[k]
        _evo(qc, la, co * tt)


TKET_TOL = 1e-8
TKET_ERRS: list[float] = []                      # max |U - e^{i phi} U_ref| of every block synthesised in this process


@lru_cache(maxsize=4096)
def tket_block(rots: tuple) -> QuantumCircuit:
    """tket block synthesis (pytket 2.18.5): PauliExpBoxes in time order -> GreedyPauliSimp(depth_weight 0, trials 32,
    seed 1, thread_timeout 1e6 so no time-dependent result) -> DecomposeBoxes -> RebaseTket ->
    FullPeepholeOptimise(allow_swaps False). Checked EXACT against the per-term product (<= TKET_TOL) or raises."""
    from pytket import Circuit
    from pytket.circuit import PauliExpBox
    from pytket.passes import DecomposeBoxes, FullPeepholeOptimise, GreedyPauliSimp, RebaseTket
    from pytket.pauli import Pauli as TP
    from pytket.qasm import circuit_to_qasm_str
    from qiskit import qasm2
    from qiskit.quantum_info import Operator
    n = len(rots[0][0])
    pm = {"X": TP.X, "Y": TP.Y, "Z": TP.Z}
    tk = Circuit(n)
    for la, a in rots:                           # PauliExpBox(P, t) = exp(-i pi/2 t P); label char i -> qubit n-1-i
        tk.add_pauliexpbox(PauliExpBox([pm[ch] for ch in la if ch != "I"], 2 * a / np.pi),
                           [n - 1 - i for i, ch in enumerate(la) if ch != "I"])
    for p in (GreedyPauliSimp(depth_weight=0.0, trials=32, seed=1, thread_timeout=1_000_000), DecomposeBoxes(),
              RebaseTket(), FullPeepholeOptimise(allow_swaps=False)):
        p.apply(tk)
    tk.replace_implicit_wire_swaps()
    out = qasm2.loads(circuit_to_qasm_str(tk, header="qelib1"), custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    ref = QuantumCircuit(n)
    for la, a in rots:
        _evo(ref, la, a)
    A, B = Operator(out).data, Operator(ref).data
    k = int(np.argmax(np.abs(B))); ph = A.flat[k] / B.flat[k]
    err = float(np.max(np.abs(A - ph / abs(ph) * B)))
    TKET_ERRS.append(err)
    if err > TKET_TOL:
        raise RuntimeError(f"tket block not exact: {err:.2e}")
    return out


def prep_I(qc: QuantumCircuit, jw="S2") -> None:
    """|I>: (psi_L^j + i psi_R^j)|I> = 0. S2/S3: pairs j >= 3 sit in |0>, j = 1, 2 is (|01> + i|10>)/sqrt2 on reg 1, 2
    (X, H, S + 1 CNOT, TW:1184-1190). Other maps: the stabilizer state of S_j = 2i psi_R^j psi_L^j = +1, synthesised."""
    if jw in ("S2", "S3") or not isinstance(jw, str):
        qc.x(rq(2)); qc.h(rq(1)); qc.s(rq(1)); qc.cx(rq(1), rq(2))
        return
    from qiskit.quantum_info import StabilizerState
    lab = majorana_labels(jw)
    stabs = []
    for j in range(1, 9):
        la, co = _product([lab["R"][j], lab["L"][j]], 1j / 2.0 * 2)
        assert abs(co.imag) < 1e-12 and abs(abs(co.real) - 1) < 1e-12
        stabs.append(("+" if co.real > 0 else "-") + la)
    qc.compose(StabilizerState.from_stabilizer_list(stabs).clifford.to_circuit(), qubits=range(8), inplace=True)


def vqa(qc: QuantumCircuit, seed: int, prep: str = "I", jw="S2", theta=None, order=None, chain=None) -> None:
    """TW:1184-1190 hardware-efficient ansatz on the register in vqa.py's convention (vqa.ORDER, vqa.CHAIN; theta
    reshaped (6, 8, 2), [..., 0] = thx, [..., 1] = thz). theta None = generic angles rng(seed) (counts only)."""
    import vqa as V
    order, chain = order or V.ORDER, chain or V.CHAIN
    if prep == "I":
        prep_I(qc, jw)
    else:
        for k in range(1, 9):
            qc.h(rq(k))
    if theta is None:
        th = np.random.default_rng(seed).uniform(-np.pi, np.pi, (N_VQA_LAYERS, 8, 2))[..., ::-1]
    else:
        th = np.asarray(theta, float).reshape(N_VQA_LAYERS, 8, 2)
    for layer in range(N_VQA_LAYERS):
        for k in range(1, 9):
            z, x = float(th[layer, k - 1, 1]), float(th[layer, k - 1, 0])
            if order == "zx":
                qc.rz(z, rq(k)); qc.rx(x, rq(k))
            else:
                qc.rx(x, rq(k)); qc.rz(z, rq(k))
        if layer < N_VQA_CHAINS:
            for ctl, tgt in V.CHAINS[chain]:                 # 0-based register indices, time order
                qc.cx(rq(ctl + 1), rq(tgt + 1))


def fswap(qc: QuantumCircuit, jw="S2") -> None:
    """fermionic SWAP_QL (eq. S6) between Q and the L injection mode (X/Y of psi_L^1,2), with its JW Z-string
    (S2/S3: mode on register qubit 2, string on qubit 1)."""
    la = majorana_labels(jw)["L"][1]
    k = la.index("X") + 1
    zs = [i + 1 for i, ch in enumerate(la) if ch == "Z"]
    qc.cx(Q, rq(k))
    for z in zs:
        qc.cz(rq(k), rq(z))
    qc.cx(Q, rq(k)); qc.swap(Q, rq(k))


def tfd_block(prep: str = "I", seed: int = 0) -> QuantumCircuit:
    qc = QuantumCircuit(8)
    vqa(qc, seed, prep)
    return qc


def step_circuit(c: Candidate) -> QuantumCircuit:
    """Table S1 style: one Trotter step of e^{-i H_L t0} (TW:1083), 8 qubits."""
    pL = h_paulis(c, "L")
    qc = QuantumCircuit(len(pL[0][0]))
    evolve(qc, c, pL, c.t0)
    return qc


def block_rotations(c: Candidate, paulis: list[tuple[str, float]], t: float, backward: bool = False) -> list[tuple[str, float]]:
    """one Trotter block as [(label, angle)] in time order, each factor exp(-i angle P) (same product as evolve)."""
    sq = trotter_seq(c, t)
    if backward:
        sq = [(k, -tt) for k, tt in (sq if c.backward == "same_order" else reversed(sq))]
    return [(paulis[k][0], paulis[k][1] * tt) for k, tt in sq]


def full_circuit(c: Candidate, with_prep: bool = True, with_evo: bool = True, block_synth=None) -> QuantumCircuit:
    """Full TW circuit, 10 qubits (P, Q, 8 register). block_synth(rotations) -> QuantumCircuit on the register
    replaces each evolution block (secondary synthesis study); None = Candidate.synth."""
    qc = QuantumCircuit(10)
    pL, pR = h_paulis(c, "L"), h_paulis(c, "R")

    def blk(p, t, back=False):
        if block_synth is None:
            evolve(qc, c, p, t, backward=back)
        else:
            qc.compose(block_synth(block_rotations(c, p, t, back)), qubits=range(8), inplace=True)
    if with_prep:
        vqa(qc, c.vqa_seed, c.prep, c.jw, c.theta, c.vqa_order, c.vqa_chain)
        qc.h(PQ); qc.cx(PQ, Q)
    if with_evo:
        blk(pL, c.t0, True)
    fswap(qc, c.jw)
    if with_evo:
        blk(pL, c.t0)
    for la, co in v_paulis(c.jw):
        _evo(qc, la, -c.mu * co)                                  # e^{i mu V}
    if with_evo:
        blk(pR, c.tread())
    return qc



# ---------------------------------------------------------------- routing
def backend(name: str):
    from qiskit_ibm_runtime import fake_provider as fp
    return {"nighthawk": fp.FakeNighthawk, "marrakesh": fp.FakeMarrakesh}[name]()


def greedy_order(terms) -> list[int]:
    """terms grouped into mutually commuting blocks (|T cap T'| even), first fit (prior verifier _greedy_order)."""
    groups: list[list[int]] = []
    for k, t in enumerate(terms):
        for g in groups:
            if all(len(set(terms[i]) & set(t)) % 2 == 0 for i in g):
                g.append(k)
                break
        else:
            groups.append([k])
    return [k for g in groups for k in g]


def count2q(tqc: QuantumCircuit) -> int:
    return sum(1 for i in tqc.data if i.operation.num_qubits == 2 and i.operation.name not in ("barrier", "measure"))


def route(qc: QuantumCircuit, be, seed: int, opt_level: int = 3, layout_method: str | None = None,
          routing_method: str | None = None, hls: str = "default") -> tuple[QuantumCircuit, int, int]:
    from qiskit.transpiler.passes import HLSConfig
    cfg = None if hls == "default" else HLSConfig(PauliEvolution=[hls])
    tq = transpile(qc, backend=be, optimization_level=opt_level, seed_transpiler=seed, approximation_degree=1.0,
                   layout_method=layout_method, routing_method=routing_method, hls_config=cfg)
    return tq, count2q(tq), tq.depth()


def to_qasm3(tqc: QuantumCircuit) -> str:
    return qasm3.dumps(tqc)


# ---------------------------------------------------------------- self checks (dense, seconds)
_M = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]], complex), "Y": np.array([[0, -1j], [1j, 0]]),
      "Z": np.diag([1.0, -1.0]).astype(complex)}


def dense(label: str) -> np.ndarray:
    out = np.array([[1.0 + 0j]])
    for ch in label:
        out = np.kron(out, _M[ch])
    return out


def majoranas(jw) -> dict[str, dict[int, np.ndarray]]:
    lab = majorana_labels(jw)
    return {s: {j: dense(lab[s][j]) / np.sqrt(2) for j in range(1, 9)} for s in "LR"}


def self_check(jw) -> None:
    """Clifford algebra; |I> annihilated; SWAP gate = eq. S6; Pauli terms = dense Majorana products."""
    from qiskit.quantum_info import Operator, Statevector
    m = majoranas(jw)
    al = [m["L"][j] for j in range(1, 9)] + [m["R"][j] for j in range(1, 9)]
    for a in range(16):
        for b in range(a, 16):
            ac = al[a] @ al[b] + al[b] @ al[a]
            assert np.allclose(ac, np.eye(256) if a == b else 0, atol=1e-12), "Clifford algebra broken"
    qc = QuantumCircuit(8); prep_I(qc, jw); sv = Statevector(qc).data
    for j in range(1, 9):
        assert np.allclose((m["L"][j] + 1j * m["R"][j]) @ sv, 0, atol=1e-12), f"|I> not annihilated at j={j}"
    lm = (m["L"][1] + 1j * m["L"][2]) / np.sqrt(2); lp = (m["L"][1] - 1j * m["L"][2]) / np.sqrt(2)   # eq. S5
    e = lambda i, j: np.outer(np.eye(2)[i], np.eye(2)[j])
    S6 = np.kron(e(0, 0), lm @ lp) + np.kron(e(0, 1), lp) + np.kron(e(1, 0), lm) + np.kron(e(1, 1), lp @ lm)
    q3 = QuantumCircuit(9); fswap(q3, jw)                                            # Q = 8, kron (Q, reg1..reg8)
    U = Operator(q3).data
    k = int(np.argmax(np.abs(S6))); ph = S6.flat[k] / U.flat[k]
    assert abs(abs(ph) - 1) < 1e-12 and np.allclose(S6, ph * U, atol=1e-12), "fswap circuit != eq. S6"
    c = eq11(jw=jw)
    for side in "LR":
        for (la, co), t, J in zip(h_paulis(c, side), c.terms, c.J):
            ref = J * m[side][t[0]] @ m[side][t[1]] @ m[side][t[2]] @ m[side][t[3]]
            assert np.allclose(co * dense(la), ref, atol=1e-12), "term Pauli != Majorana product"


def end_to_end_check(c: Candidate, hls: str = "rustiq") -> tuple[float, float]:
    """I_PT of full_circuit(c) without the VQA block, started from sykg's exact |Bell>_PQ |TFD> in the same map, vs
    sykg's Trotter I_PT (one time point): the circuit and the physics library must agree to ~1e-10."""
    import sykg
    from qiskit.quantum_info import Statevector, entropy, partial_trace
    proto = "fix" if c.protocol == "FIX" else "sym"          # physics is map independent (V3)
    m = sykg.TwoSided(c.terms, c.J, proto, order=c.seq())
    mj = majoranas(c.jw)
    H = sum(J * (mj["L"][t[0]] @ mj["L"][t[1]] @ mj["L"][t[2]] @ mj["L"][t[3]]
                 + mj["R"][t[0]] @ mj["R"][t[1]] @ mj["R"][t[2]] @ mj["R"][t[3]]) for t, J in zip(c.terms, c.J))
    qi = QuantumCircuit(8); prep_I(qi, c.jw); iv = Statevector(qi).data
    w, v = np.linalg.eigh(H); tfd = v @ (np.exp(-sykg.BETA * w / 4) * (v.conj().T @ iv)); tfd /= np.linalg.norm(tfd)
    init = np.kron(np.array([1, 0, 0, 1]) / np.sqrt(2), tfd)                       # |Bell>_PQ (x) |TFD>, own map
    # synthesise first: an unsynthesised whole-block PauliEvolutionGate simulates as the EXACT exponential, not the
    # Trotter product the transpiler ships (checked 2026-10-22: 0.3357 exact-like vs 0.2825 Trotter on a K=24 draw)
    from qiskit.transpiler.passes import HLSConfig
    qc = transpile(full_circuit(c, with_prep=False), basis_gates=["cx", "cz", "swap", "rz", "sx", "x"], optimization_level=0,
                   hls_config=None if hls == "default" else HLSConfig(PauliEvolution=[hls]))
    sv = Statevector(init).evolve(qc)
    rho = partial_trace(sv, [k for k in range(10) if k not in (rq(1), PQ)])          # keep P and register qubit 1
    mi = entropy(partial_trace(rho, [0]), base=2) + entropy(partial_trace(rho, [1]), base=2) - entropy(rho, base=2)
    tr = {"lie1": "lie1", "suzuki2": "suz2"}[c.trot()]
    ref = float(m.ipt(c.mu, [c.tread()], tr, c.t0)[0])
    return float(mi), ref
