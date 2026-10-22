# b8_tb_17_25_10200118_1444_409_25_5_443_610  [CERTIFIED]

K = 25; term set 4bb5b414c16fff2e; source: frozen verifier verdict (j-b4db6e0ae575b105)

Files: `hamiltonian.json` (terms, J, order, times, mu, VQA angles, layout, seeds, conventions), `FIX.qasm` / `SYM.qasm` (shipped routed circuits, mu = -12, operating time), `FIX_mu+12.qasm` / `SYM_mu+12.qasm` (same seed and time, mu = +12, when built), `panels.pdf`, `row.json`.

| protocol | t | VQA F | eta | 2q | 1q | depth | transpiler seed |
|---|---|---|---|---|---|---|---|
| FIX | 2.4 | 0.9667 | 0.0896 | 209 | 907 | 489 | 1016 |
| SYM | 2.35 | 0.9667 | 0.0067 | 286 | 1185 | 688 | 1017 |

Regenerate for another time or mu sign: `python m3/make_circuits.py candidates/b8_tb_17_25_10200118_1444_409_25_5_443_610/hamiltonian.json --protocol FIX --t 2.0 --mu 12 --route`.

Conventions: see `hamiltonian.json` -> conventions (JW map S2 for both protocols; your SYM runs used S3).

Submit the shipped QASM as is (no re-transpile; or initial_layout + optimization_level=0). J_T multiplies the ordered Majorana product with i^(q/2) absorbed (Eq. 11 as printed); Eq. 6's extra sign drops F to about 0.45. FakeNighthawk's map is assumed identical to ibm_phoenix. Measure P and register qubit 1 (physical qubits in each QASM header).
