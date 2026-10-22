# p324c_tb_23_24_44300310_1429_155_129 (listing)

Listing certified by the frozen verifier that also passed the frozen noise-comparison check (cost, robustness, noisy dI vs Eq. 11).

Files: `hamiltonian.json`, `FIX_mu-12.qasm`, `FIX_mu+12.qasm`, `SYM_mu-12.qasm`, `SYM_mu+12.qasm` (routed for FakeNighthawk / ibm_phoenix map, physical-qubit syntax; header = layout, seed, counts), `predictions.json`.

| protocol | t0 / t1 | seed | 2q (mu -12 / +12) | 1q | depth | noiseless dI (circuit) | retention pred | noisy lists A / B / C |
|---|---|---|---|---|---|---|---|---|
| FIX | 1.8 / 2.45 | 1019 | 216/216 | 922/925 | 481/482 | 0.27 | 0.093 | 0.097 / 0.094 / 0.10 |
| SYM | 2.15 / 2.15 | 1018 | 225/225 | 1049/1052 | 501/501 | 0.25 | 0.085 | 0.087 / 0.086 / 0.088 |

Stats: r_odd 0.4957, r_even 0.3961, sff_rel_DENSE 1.577, otoc_E3 2.56, eta_FIX 0.1469, eta_SYM 0.008961, vqa_F 0.9722, selection_one_in 1.802e+13, term_set_id d6a76272f3fb45eb, FIX_shipped_seed_gate_ro 0.1021, SYM_shipped_seed_gate_ro 0.08545

Sources: frozen noise-comparison verdict (check version 370dcd31; jobs core:j-be96ba8a65030818,eq11:j-70ce5bd2ed4bcd60+j-c0831ef1138f1f77,m4c:j-7b057abe7268073c+j-3120b27978351caf,noise:j-ddf53639327eab62+j-af00b3a565523fc9); candidates/p324c_tb_23_24_44300310_1429_155_129/row.json. Noisy cites per seed set: see predictions.json.
