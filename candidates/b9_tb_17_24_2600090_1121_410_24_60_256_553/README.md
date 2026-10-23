# b9_tb_17_24_2600090_1121_410_24_60_256_553  [CERTIFIED]

K = 24; term set d122addbfb06f85b; source: frozen verifier verdict (j-d1809ee5be4bd840)

Files: `hamiltonian.json` (terms, J, order, times, mu, VQA angles, layout, seeds, conventions), `FIX.qasm` / `SYM.qasm` (shipped routed circuits, mu = -12, operating time), `FIX_mu+12.qasm` / `SYM_mu+12.qasm` (same seed and time, mu = +12, when built), `panels.pdf`, `row.json`.

| protocol | t | VQA F | eta | 2q | 1q | depth | transpiler seed |
|---|---|---|---|---|---|---|---|
| FIX | 2.55 | 0.9625 | 0.1159 | 215 | 910 | 463 | 1009 |
| SYM | 2.35 | 0.9625 | 0.0043 | 276 | 1188 | 685 | 1011 |

Rebuild and check the shipped circuits: `uv run route candidates/b9_tb_17_24_2600090_1121_410_24_60_256_553`.

Conventions: see `hamiltonian.json` -> conventions (JW map S2 for both protocols; arXiv:2604.10090 ran SYM in S3).

Submit the shipped QASM as is (no re-transpile; or initial_layout + optimization_level=0). J_T multiplies the ordered Majorana product with i^(q/2) absorbed (Eq. 11 as printed); Eq. 6's extra sign drops F to about 0.45. FakeNighthawk's coupling map is assumed to match the Nighthawk device. Measure P and register qubit 1 (physical qubits in each QASM header).
