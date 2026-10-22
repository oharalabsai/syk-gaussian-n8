# p325_tb_23_24_44300310_1429_19277_201  [CERTIFIED]

K = 24; term set d6a76272f3fb45eb; source: frozen verifier verdict (inside the frozen noise-comparison verdict)

Files: `hamiltonian.json` (terms, J, order, times, mu, VQA angles, layout, seeds, conventions), `FIX.qasm` / `SYM.qasm` (shipped routed circuits, mu = -12, operating time), `FIX_mu+12.qasm` / `SYM_mu+12.qasm` (same seed and time, mu = +12, when built), `panels.pdf`, `row.json`.

| protocol | t | VQA F | eta | 2q | 1q | depth | transpiler seed |
|---|---|---|---|---|---|---|---|
| FIX | 2.6 | 0.9654 | 0.0936 | 216 | 925 | 471 | 1019 |
| SYM | 2.35 | 0.9654 | 0.0017 | 225 | 1069 | 509 | 1018 |

Regenerate for another time or mu sign: `python m3/make_circuits.py candidates/p325_tb_23_24_44300310_1429_19277_201/hamiltonian.json --protocol FIX --t 2.0 --mu 12 --route`.

Conventions: see `hamiltonian.json` -> conventions (JW map S2 for both protocols; your SYM runs used S3).

Submit the shipped QASM as is (no re-transpile; or initial_layout + optimization_level=0). J_T multiplies the ordered Majorana product with i^(q/2) absorbed (Eq. 11 as printed); Eq. 6's extra sign drops F to about 0.45. FakeNighthawk's map is assumed identical to ibm_phoenix. Measure P and register qubit 1 (physical qubits in each QASM header).
