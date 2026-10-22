# eq11_greedy (control+)

positive control: the published sparse model, same pipeline as the listings

Files: `hamiltonian.json`, `FIX_mu-12.qasm`, `FIX_mu+12.qasm`, `SYM_mu-12.qasm`, `SYM_mu+12.qasm` (routed for FakeNighthawk / ibm_phoenix map, physical-qubit syntax; header = layout, seed, counts), `predictions.json`.

| protocol | t0 / t1 | seed | 2q (mu -12 / +12) | 1q | depth | noiseless dI (circuit) | retention pred | noisy lists A / B / C |
|---|---|---|---|---|---|---|---|---|
| FIX | 1.8 / 2.35 | 1007 | 175/175 | 697/696 | 345/345 | 0.21 | 0.090 | 0.089 / 0.089 / 0.090 |
| SYM | 2.05 / 2.05 | 1004 | 215/215 | 819/818 | 397/396 | 0.20 | 0.073 | 0.074 / 0.075 / 0.077 |

Stats: FIX_shipped_seed_gate_ro 0.08288, FIX_shipped_seed_ideal 0.2118, SYM_shipped_seed_gate_ro 0.07415, SYM_shipped_seed_ideal 0.2003

Sources: Eq. 11 of arXiv:2604.10090 (GIST), greedy term order, pinned S2 fit src/sykn8/v2_bars/vqa_eq11.json (F 0.9825); noisy-dI reference jobs j-70ce5bd2ed4bcd60 (FIX) / j-c0831ef1138f1f77 (SYM); shipped QASM from kit job j-5f38416434f460bf. Noisy cites per seed set: see predictions.json.
