# b12_tb_17_25_13000048_2333_412_25_29_983_124  [CERTIFIED]

**FLAG: certified under the frozen 200-draw bar; fails the 1000-draw bar** (SYM routing-seed list B 303 vs 299.2, CI 297.5-301.7). Bar-noise sensitive.

**FLAG: set 13000048_2333 fails the noise-comparison cost check** (routed 2q above its cost bar on routing-seed lists B and C; frozen noise-comparison verdict). Certified by the frozen verifier, not counted in the noise comparison.

K = 25; term set 01b43c156d1e6382; source: frozen verifier verdict (j-62486d7315243062)

Files: `hamiltonian.json` (terms, J, order, times, mu, VQA angles, layout, seeds, conventions), `FIX.qasm` / `SYM.qasm` (shipped routed circuits, mu = -12, operating time), `FIX_mu+12.qasm` / `SYM_mu+12.qasm` (same seed and time, mu = +12, when built), `panels.pdf`, `row.json`.

| protocol | t | VQA F | eta | 2q | 1q | depth | transpiler seed |
|---|---|---|---|---|---|---|---|
| FIX | 2.4 | 0.9773 | 0.1352 | 218 | 904 | 480 | 1005 |
| SYM | 2.1 | 0.9773 | 0.0077 | 299 | 1249 | 678 | 1015 |

Rebuild and check the shipped circuits: `uv run route candidates/b12_tb_17_25_13000048_2333_412_25_29_983_124`.

Conventions: see `hamiltonian.json` -> conventions (JW map S2 for both protocols; arXiv:2604.10090 ran SYM in S3).

Submit the shipped QASM as is (no re-transpile; or initial_layout + optimization_level=0). J_T multiplies the ordered Majorana product with i^(q/2) absorbed (Eq. 11 as printed); Eq. 6's extra sign drops F to about 0.45. FakeNighthawk's coupling map is assumed to match the Nighthawk device. Measure P and register qubit 1 (physical qubits in each QASM header).
