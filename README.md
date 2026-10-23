# syk-gaussian-n8

Sparse Gaussian SYK Hamiltonians (N = 8 Majoranas, K = 24..25 terms) selected for a large traversable-wormhole signal
(Delta I_PT) at a low two-qubit cost after routing on IBM Nighthawk, with ready-to-run OpenQASM 3 circuits.
An instance is one term set plus one Gaussian coupling draw (the code and file names call it a "listing").
This repository holds the model, the two-sided protocol, the verifier, routing and the FakeNighthawk noise model, plus
every instance, so each reported noiseless number can be recomputed on a laptop. Ohara Labs, October 2026.

## Result

- One N = 8 Gaussian sparse SYK instance with K = 24 terms (`kit/p325_tb_23_24_44300310_1429_19277_201`, backup 155_129 on the same term set) routes on FakeNighthawk to 216 (FIX) / 225 (SYM) two-qubit gates: 23% / 5% more than the binary Eq. 11 of arXiv:2604.10090 in the same pipeline (175 / 215), 30% / 45% fewer than a typical K = 24 draw (309 / 410).
- In a FakeNighthawk noise simulation the FIX circuit keeps 25-35% more teleportation signal than Eq. 11 (18% with term order optimized on both sides); SYM keeps 25-29% more in canonical order and ties with optimized order.
- Simulation only; nothing has run on hardware. The note is `paper/gaussian-sparse-syk-n8-nighthawk-teleportation.pdf`; `uv run verify` reproduces its noiseless numbers in about a minute.

Conventions follow arXiv:2604.10090 (beta = 3, |mu| = 12, {psi_a, psi_b} = delta_ab, base-2 entropies, JW map S2 of its
Eq. S2 for both protocols, FIX t0 = 1.8 and SYM t0 = t1).

## Install

Requires [uv](https://docs.astral.sh/uv/) (it fetches Python 3.12 if needed).

```bash
git clone <this repo> && cd syk-gaussian-n8
uv sync            # exact pins from uv.lock: numpy 2.5.2, scipy 1.18.1, qiskit 2.5.2, qiskit-ibm-runtime 0.49.0,
                   # qiskit-aer 0.18.0, qiskit-qasm3-import 0.6.0, pytket 2.18.5, pydantic 2.13.4
```

## Reproduce

| command | what it does | time |
|---|---|---|
| `uv run verify` | Eq. 11 known answers (r_odd 0.537163, r_even 0.478812, model self test, exact Delta I_PT peaks of Eq. 11 vs Fig. S1 of arXiv:2604.10090, SYM 0.212 and FIX 0.276, to 0.003 bits) and, for every instance in `candidates/` and `kit/`, the recomputed noiseless numbers vs the reported ones: r_odd, r_even, SFF / DENSE, the 8 Delta I_PT peaks and 2 revivals, VQA fidelity F and Trotter error eta of the listed angles, OTOC E3, the shipped circuit's noiseless dI at its readout time, and the two-qubit count of every shipped QASM file. Writes `out/verify.json`, exits 0 iff all pass. `--only <substring>` checks a subset of instances. | about 1 min on 8 cores |
| `uv run route [instance]` | rebuilds the instance's fitted circuit (tket block synthesis), routes it on FakeNighthawk at the shipped transpiler seed (optimization level 3, approximation degree 1) and compares the 2q count and the QASM 3 body with the shipped files, mu = -12 and +12, FIX and SYM. Default: the kit's recommended instance `p325_tb_23_24_44300310_1429_19277_201`. | about 1 min |
| `uv run noise [instance] --protocol FIX` | FakeNighthawk noisy dI = I_PT(mu -12) - I_PT(mu +12) of the shipped circuit pair (density matrix on the active qubits, gate errors incl. T1 / T2 on gate durations, readout as a confusion map on two-qubit Pauli tomography). Checks the ideal dI against the model (1e-5) and the noisy dI (gate and readout noise) against the kit's `predictions.json` (1e-6). | 1 to a few min |

`[instance]` is a folder path, a folder name or a unique substring of one (for example `19277_201`). The three kit
instances also sit in `candidates/`; a name or substring that matches both resolves to `kit/<id>` (same circuits and
`hamiltonian.json`, plus the `predictions.json` that `noise` checks against). Pass `candidates/<id>` to pick that copy.

Expected `uv run verify` result: `known_ok 1  checks 441  failed 0`; instance values agree with the reported ones to
<= 4e-14, except the reported shipped-circuit dI0, which is a simulation of the routed circuits and agrees to <= 6e-7
(tolerance 1e-5). Expected `uv run noise` (default instance, FIX): noisy dI 0.1199.

Routing is platform sensitive: the shipped circuits were routed on x86_64 Linux, where `uv run route` reproduces both
their two-qubit counts and their QASM 3 text exactly. On another architecture the same transpiler seed can route
differently: in one local check on an arm64 Mac, one FIX circuit routed to 199 instead of 200 two-qubit gates, so its
QASM text differed as well. `uv run verify` counts the gates of the shipped files, so its two-qubit checks and all
physics numbers do not depend on the platform.

## Contents

| path | content |
|---|---|
| `src/sykn8/sykg.py` | the N = 8 model: Hamiltonians, spectra and parity sectors, r, SFF, two-sided protocol (exact and Trotter I_PT), JW maps, ensemble draws |
| `src/sykn8/circuits.py` | full 10-qubit circuit (|I>, VQA block, Bell pair, SWAP_QL, e^(i mu V), Trotter blocks), tket block synthesis, routing on FakeNighthawk |
| `src/sykn8/vqa.py`, `diag.py` | VQA TFD (6 Rx Rz layers, 5 linear CNOT layers, 96 angles); OTOC, two-point correlator, size winding |
| `src/sykn8/verify.py` | the frozen verifier that judged every instance, as a library (job entry points removed). Its checks: chaos (level-spacing ratio r in both parity sectors, spectral form factor, no early revival), signal (the 8 noiseless Delta I_PT peaks at least Eq. 11's), the variational TFD fit (fidelity F >= 0.96) and Trotter error, OTOC and size-winding chaos diagnostics, and routed two-qubit cost (at most 0.70x the typical median). Its sha256 pins match this repository's files (circuits.py, diag.py, sykg.py, vqa.py and v2_bars/l1.json, l5.json, m2e.json re-pinned after comment and label-string edits; the run's original hashes are noted next to the pins) |
| `src/sykn8/reproduce.py`, `route.py`, `noise.py` | the three commands above |
| `src/sykn8/v2_bars/` | frozen thresholds the verifier reads (chaos bands, typical baseline, cost bars and routing setting, Eq. 11 VQA fit) |
| `candidates/<id>/` | 16 instances with their verifier verdicts: 14 certified and 2 that failed only the SYM routed-cost check on seed list B (13800064_1823 and 16600062_1987, `status` FAILED in `row.json` and `comparison.csv`). 7 of the 14 certified ones are among the 23 that also passed the robustness gates (see below). The borderline 1718_714's row in `comparison.csv` is a placeholder: its reported numbers are in its `verdict.json`. Each folder: `hamiltonian.json` (terms, couplings, protocols, times, VQA angles, ship seed and layout), `row.json` (reported numbers), shipped QASM (`FIX.qasm`, `SYM.qasm` at mu = -12, `*_mu+12.qasm`), run card; `comparison.csv` |
| `kit/<id>/` | the hardware kit: three instances (19277_201, recommended; 155_129; and the borderline 1718_714, which beats Eq. 11 under noise on routing-seed list A only), Eq. 11 (positive control) and two typical K = 24 draws (negative controls), each with `hamiltonian.json`, `predictions.json` and the four QASM files; report-only `t1_scan.json` (readout-time scan) in the two robust instances and Eq. 11, and `kit/chaos_curves.json` (SFF, two-point function, size winding); S3-map variants under `S3/`; `kit/README.md` says what to run and how to read it |
| `paper/gaussian-sparse-syk-n8-nighthawk-teleportation.pdf` | the paper |

How many instances: the frozen verifier certified 32 instances on 16 term sets. 23 of them, on 14 term sets, also
passed the paper's robustness gates (the cost bar on three routing-seed lists, eight noiseless Delta I_PT peaks at least
1.02x Eq. 11's, stability under 1% coupling noise, r >= 0.386 in both parity sectors). Under FakeNighthawk noise one of
those 14 sets (44300310_1429: instances 19277_201 and 155_129) beats Eq. 11 on all three lists, and 1718_714 (set
16500233_1414) is borderline. `candidates/` packages 14 of the 32 certified instances (7 of the 23) and the 2 rejected
ones; `kit/` carries the three noise-tested instances and the controls.

Running the shipped circuits: submit the ISA circuit as is (no re-transpile), or pass initial_layout and
optimization_level=0. Measure P and register qubit 1 (named in each QASM header). Load with
`qiskit.qasm3.loads(open(path).read())`.

## Citation

Model and protocol: arXiv:2604.10090. This repository: Gyalpo Dongo and Sebastian Alberdi, Ohara Labs, syk-gaussian-n8 (2026).

## License

MIT, see `LICENSE`.
