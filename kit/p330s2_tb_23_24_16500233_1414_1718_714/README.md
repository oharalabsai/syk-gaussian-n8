# p330s2_tb_23_24_16500233_1414_1718_714 (listing)

BORDERLINE option (second term set 16500233_1414): passes the pre-registered rule (frozen noise-comparison check, routing-seed list A: FIX +2.8%, SYM +0.8% vs Eq. 11), but the SYM margin is within seed noise and it falls below Eq. 11 on the held-out seed lists (list B FIX -0.05% / SYM -1.4%; list C -3.6% / -3.1%). 6 listings reached the verifier, 1 passed. Prefer set 44300310_1429 (listing 19277_201 first).

Files: `hamiltonian.json`, `FIX_mu-12.qasm`, `FIX_mu+12.qasm`, `SYM_mu-12.qasm`, `SYM_mu+12.qasm` (routed for FakeNighthawk / ibm_phoenix map, physical-qubit syntax; header = layout, seed, counts), `predictions.json`.

| protocol | t0 / t1 | seed | 2q (mu -12 / +12) | 1q | depth | noiseless dI (circuit) | retention pred | noisy lists A / B / C |
|---|---|---|---|---|---|---|---|---|
| FIX | 1.8 / 2.55 | 1015 | 208/208 | 927/930 | 423/425 | 0.25 | 0.089 | 0.092 / 0.088 / 0.087 |
| SYM | 2.5 / 2.5 | 1010 | 250/250 | 1128/1127 | 577/579 | 0.26 | 0.078 | 0.075 / 0.074 / 0.074 |

Stats: r_odd 0.5635, r_even 0.647, sff_rel_DENSE 0.8893, otoc_E3 2.84, eta_FIX 0.09645, eta_SYM 0.005471, vqa_F 0.9705, selection_one_in 1.175e+13, term_set_id 2f2d8f4f11f259b2, FIX_shipped_seed_gate_ro 0.09352, SYM_shipped_seed_gate_ro 0.07478

Sources: frozen noise-comparison verdict (check version 370dcd31; jobs core:j-a2c2b10d8136a413,eq11:j-70ce5bd2ed4bcd60+j-c0831ef1138f1f77,m4c:j-2eb00e82d036beb4+j-0d7e817b462788df,noise:j-552809eff16dab92+j-9d5acb4db6e413c6); candidates/p330s2_tb_23_24_16500233_1414_1718_714/row.json. Noisy cites per seed set: see predictions.json.
