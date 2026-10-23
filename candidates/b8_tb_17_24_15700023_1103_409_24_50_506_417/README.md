# b8_tb_17_24_15700023_1103_409_24_50_506_417  [CERTIFIED]

K = 24; term set 618253002363307d; source: frozen verifier verdict (j-a741f0f731f0afb2)

Files: `hamiltonian.json` (terms, J, order, times, mu, VQA angles, layout, seeds, conventions), `FIX.qasm` / `SYM.qasm` (shipped routed circuits, mu = -12, operating time), `FIX_mu+12.qasm` / `SYM_mu+12.qasm` (same seed and time, mu = +12, when built), `panels.pdf`, `row.json`.

| protocol | t | VQA F | eta | 2q | 1q | depth | transpiler seed |
|---|---|---|---|---|---|---|---|
| FIX | 2.5 | 0.9808 | 0.1061 | 214 | 888 | 409 | 1019 |
| SYM | 2.4 | 0.9808 | 0.0049 | 272 | 1155 | 689 | 1002 |

Rebuild and check the shipped circuits: `uv run route candidates/b8_tb_17_24_15700023_1103_409_24_50_506_417`.

Conventions: see `hamiltonian.json` -> conventions (JW map S2 for both protocols; arXiv:2604.10090 ran SYM in S3).

Submit the shipped QASM as is (no re-transpile; or initial_layout + optimization_level=0). J_T multiplies the ordered Majorana product with i^(q/2) absorbed (Eq. 11 as printed); Eq. 6's extra sign drops F to about 0.45. FakeNighthawk's coupling map is assumed to match the Nighthawk device. Measure P and register qubit 1 (physical qubits in each QASM header).
