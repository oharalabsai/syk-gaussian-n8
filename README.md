# syk-gaussian-n8

Gaussian sparse SYK Hamiltonians (N = 8 Majoranas, K = 24 terms) for the traversable-wormhole teleportation protocol of
[arXiv:2604.10090](https://arxiv.org/abs/2604.10090), with routed circuits ready to run on IBM Nighthawk.
Everything here is simulation; nothing has run on hardware yet.

## Result

- **Recommended instance:** `kit/p325_tb_23_24_44300310_1429_19277_201` ("19277_201"). The backup, 155_129, uses the same
  24 terms with different couplings.
- **Cost:** 216 (FIX) / 225 (SYM) two-qubit gates after routing on FakeNighthawk. That is more than the binary Eq. 11 of
  arXiv:2604.10090 compiled the same way (175 / 215), but 30% / 45% fewer than a typical K = 24 draw (309 / 410).
- **Signal under noise:** in a FakeNighthawk noise simulation, the FIX circuit keeps 25-35% more teleportation signal
  than Eq. 11, or 18% with the term order optimized on both sides. SYM keeps 25-29% more in the default order and ties
  once the order is optimized.
- **Eq. 11 in its authors' own term orders** (private communication, 23 Oct 2026): in map S2 our margins are 23-50%
  (FIX, the backward convention worse for us) and 51-67% (SYM). In map S3, where the paper ran SYM, SYM only ties, so
  the SYM claim is for S2 only.

The note is `paper/gaussian-sparse-syk-n8-nighthawk-teleportation.pdf`.

## Install

Requires [uv](https://docs.astral.sh/uv/), which fetches Python 3.12 if needed.

```bash
git clone https://github.com/oharalabsai/syk-gaussian-n8 && cd syk-gaussian-n8
uv sync    # exact versions from uv.lock (qiskit 2.5.2, qiskit-ibm-runtime 0.49.0, qiskit-aer 0.18.0, pytket 2.18.5)
```

## Check the numbers

```bash
uv run verify                    # about 1 min: every noiseless number in the note
uv run route 19277_201           # about 1 min: rebuild and route the circuits, compare with the shipped files
uv run noise 19277_201 --protocol FIX   # a few min: the FakeNighthawk noisy signal
```

- **`verify`** first checks the model against Eq. 11's known answers in arXiv:2604.10090: level-spacing ratios
  0.537 / 0.479 and Fig. S1's signal peaks 0.212 (SYM) / 0.276 (FIX). It also pins our model's Eq. 11 peaks in its
  authors' term orders (regression values; the paper prints none).
  It then recomputes every instance's chaos, signal, TFD-fit and gate-count numbers and compares them with the
  reported ones. The expected last line is `known_ok 1  checks 444  failed 0`. It writes `out/verify.json`.
- **`route`** rebuilds an instance's circuits and routes them on FakeNighthawk at the shipped seed. Routing depends on the
  platform: the shipped files were routed on x86_64 Linux, and an arm64 Mac can differ by a gate. `verify` counts gates in
  the shipped files, so it does not depend on the platform.
- **`noise`** simulates the shipped circuit pair with FakeNighthawk gate and readout noise. The expected result for
  19277_201 FIX is a noisy signal of 0.1199 bits.

An instance argument can be a folder path, a folder name, or any unique part of one (a name found in both `kit/` and
`candidates/` resolves to `kit/`).

## Running on hardware

Start with `kit/README.md`: what to run, what to measure, shot counts and predictions. In short:
- Submit the QASM files as they are. Do not transpile them again; if your stack requires it, pass the layout from the
  header with `optimization_level=0`.
- Measure P and register qubit 1. Each file's header names the physical qubits.
- Load a file with `qiskit.qasm3.loads(open(path).read())`.
- Run Eq. 11 (`kit/eq11_greedy`) in the same session as the control. It is Eq. 11 in our greedy term order, which in our
  compile is a harder bar than the authors' own orders.

## What is here

| path | content |
|---|---|
| `kit/` | the hardware kit: the three noise-tested instances (19277_201, 155_129, and the borderline 1718_714), Eq. 11 and two typical draws as controls |
| `candidates/` | 16 instances with their verifier verdicts: 14 that passed every check and 2 that failed only a gate-count bar |
| `src/sykn8/` | the model and protocol (`sykg.py`), circuits and routing (`circuits.py`), TFD preparation (`vqa.py`), chaos diagnostics (`diag.py`), the verifier that judged every instance (`verify.py`), and the three commands |
| `paper/` | the note |

Conventions follow arXiv:2604.10090: beta = 3, mu = -/+12, FIX uses t0 = 1.8, SYM uses t0 = t1, base-2 entropies.
We use its Jordan-Wigner map S2 for both protocols; the paper ran SYM in S3.

## Citation

Model and protocol: M. Byun, H. Lee and K.-Y. Kim, arXiv:2604.10090. This repository: Gyalpo Dongo and Sebastian
Alberdi, Ohara Labs, syk-gaussian-n8 (2026).

## License

MIT, see `LICENSE`.
