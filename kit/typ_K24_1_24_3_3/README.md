# typ_K24_1_24_3_3 (control-)

negative control: a typical uniform K24 Gaussian draw (no selection), same pipeline

Files: `hamiltonian.json`, `FIX_mu-12.qasm`, `FIX_mu+12.qasm`, `SYM_mu-12.qasm`, `SYM_mu+12.qasm` (routed for FakeNighthawk / ibm_phoenix map, physical-qubit syntax; header = layout, seed, counts), `predictions.json`.

| protocol | t0 / t1 | seed | 2q (mu -12 / +12) | 1q | depth | noiseless dI (circuit) | retention pred | noisy lists A / B / C |
|---|---|---|---|---|---|---|---|---|
| FIX | 1.8 / 2.6 | 1019 | 288/288 | 1054/1057 | 567/566 | 0.065 | 0.015 | 0.015 |
| SYM | 3.7 / 3.7 | 1001 | 406/406 | 1504/1513 | 907/911 | 0.12 | 0.015 | 0.013 |

Sources: typical draw typ_K24_1_24_3_3 (routed typical-draw table); 720-start S2 fit; noisy dI on routing-seed list A, central five seeds, job j-8d825c8aff59355b; shipped QASM from kit job j-5f38416434f460bf. Noisy cites per seed set: see predictions.json.
