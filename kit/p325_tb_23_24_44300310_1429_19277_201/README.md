# p325_tb_23_24_44300310_1429_19277_201 (listing)

Listing certified by the frozen verifier that also passed the frozen noise-comparison check (cost, robustness, noisy dI vs Eq. 11).

Files: `hamiltonian.json`, `FIX_mu-12.qasm`, `FIX_mu+12.qasm`, `SYM_mu-12.qasm`, `SYM_mu+12.qasm` (routed for FakeNighthawk / ibm_phoenix map, physical-qubit syntax; header = layout, seed, counts), `predictions.json`.

| protocol | t0 / t1 | seed | 2q (mu -12 / +12) | 1q | depth | noiseless dI (circuit) | retention pred | noisy lists A / B / C |
|---|---|---|---|---|---|---|---|---|
| FIX | 1.8 / 2.6 | 1019 | 216/216 | 925/927 | 471/472 | 0.33 | 0.11 | 0.12 / 0.11 / 0.12 |
| SYM | 2.35 / 2.35 | 1018 | 225/225 | 1069/1071 | 509/509 | 0.28 | 0.095 | 0.096 / 0.095 / 0.096 |

Stats: r_odd 0.4437, r_even 0.4507, sff_rel_DENSE 1.266, otoc_E3 2.94, eta_FIX 0.09363, eta_SYM 0.001687, vqa_F 0.9654, selection_one_in 1.518e+13, term_set_id d6a76272f3fb45eb, FIX_shipped_seed_gate_ro 0.1199, SYM_shipped_seed_gate_ro 0.09385

Sources: frozen noise-comparison verdict (check version 370dcd31; jobs core:j-52bd7902571bc219,eq11:j-70ce5bd2ed4bcd60+j-c0831ef1138f1f77,m4c:j-9f103849485d9e5a+j-f24ce419a53b46e7,noise:j-e10e45a66c692933+j-266bf0011bbfca30); candidates/p325_tb_23_24_44300310_1429_19277_201/row.json. Noisy cites per seed set: see predictions.json.
