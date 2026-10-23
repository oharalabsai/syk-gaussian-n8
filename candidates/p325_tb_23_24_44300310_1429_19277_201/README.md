# p325_tb_23_24_44300310_1429_19277_201  [CERTIFIED]

K = 24; term set d6a76272f3fb45eb; source: frozen verifier verdict (inside the frozen noise-comparison verdict)

Files: `hamiltonian.json` (terms, J, order, times, mu, VQA angles, layout, seeds, conventions), `FIX.qasm` / `SYM.qasm` (shipped routed circuits, mu = -12, operating time), `FIX_mu+12.qasm` / `SYM_mu+12.qasm` (same seed and time, mu = +12, when built), `panels.pdf`, `row.json`.

| protocol | t | VQA F | eta | 2q | 1q | depth | transpiler seed |
|---|---|---|---|---|---|---|---|
| FIX | 2.6 | 0.9654 | 0.0936 | 216 | 925 | 471 | 1019 |
| SYM | 2.35 | 0.9654 | 0.0017 | 225 | 1069 | 509 | 1018 |

Rebuild and check the shipped circuits: `uv run route candidates/p325_tb_23_24_44300310_1429_19277_201`.

Conventions: see `hamiltonian.json` -> conventions (JW map S2 for both protocols; arXiv:2604.10090 ran SYM in S3).

Submit the shipped QASM as is (no re-transpile; or initial_layout + optimization_level=0). J_T multiplies the ordered Majorana product with i^(q/2) absorbed (Eq. 11 as printed); Eq. 6's extra sign drops F to about 0.45. FakeNighthawk's coupling map is assumed to match the Nighthawk device. Measure P and register qubit 1 (physical qubits in each QASM header).
