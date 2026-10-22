# typ_K24_1_24_1_1 (control-)

negative control: a typical uniform K24 Gaussian draw (no selection), same pipeline

Files: `hamiltonian.json`, `FIX_mu-12.qasm`, `FIX_mu+12.qasm`, `SYM_mu-12.qasm`, `SYM_mu+12.qasm` (routed for FakeNighthawk / ibm_phoenix map, physical-qubit syntax; header = layout, seed, counts), `predictions.json`.

| protocol | t0 / t1 | seed | 2q (mu -12 / +12) | 1q | depth | noiseless dI (circuit) | retention pred | noisy lists A / B / C |
|---|---|---|---|---|---|---|---|---|
| FIX | 1.8 / 2.5 | 1004 | 307/307 | 1152/1153 | 619/621 | 0.064 | 0.013 | 0.017 |
| SYM | 2.5 / 2.5 | 1000 | 471/471 | 1721/1717 | 1040/1041 | 0.064 | 0.0055 | 0.0069 |

Sources: typical draw typ_K24_1_24_1_1 (routed typical-draw table); 720-start S2 fit; noisy dI on routing-seed list A, central five seeds, job j-8d825c8aff59355b; shipped QASM from kit job j-5f38416434f460bf. Noisy cites per seed set: see predictions.json.
