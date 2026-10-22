# b8_tb_17_24_7900071_455_409_24_2_548_259  [CERTIFIED]

K = 24; term set dd9fbc0b8598f171; source: frozen verifier verdict (j-eba08eb4157100e8)

Files: `hamiltonian.json` (terms, J, order, times, mu, VQA angles, layout, seeds, conventions), `FIX.qasm` / `SYM.qasm` (shipped routed circuits, mu = -12, operating time), `FIX_mu+12.qasm` / `SYM_mu+12.qasm` (same seed and time, mu = +12, when built), `panels.pdf`, `row.json`.

| protocol | t | VQA F | eta | 2q | 1q | depth | transpiler seed |
|---|---|---|---|---|---|---|---|
| FIX | 2.5 | 0.9727 | 0.1072 | 200 | 867 | 409 | 1008 |
| SYM | 2.2 | 0.9727 | 0.0022 | 261 | 1125 | 557 | 1004 |

Regenerate for another time or mu sign: `python m3/make_circuits.py candidates/b8_tb_17_24_7900071_455_409_24_2_548_259/hamiltonian.json --protocol FIX --t 2.0 --mu 12 --route`.

Conventions: see `hamiltonian.json` -> conventions (JW map S2 for both protocols; your SYM runs used S3).

Submit the shipped QASM as is (no re-transpile; or initial_layout + optimization_level=0). J_T multiplies the ordered Majorana product with i^(q/2) absorbed (Eq. 11 as printed); Eq. 6's extra sign drops F to about 0.45. FakeNighthawk's map is assumed identical to ibm_phoenix. Measure P and register qubit 1 (physical qubits in each QASM header).
