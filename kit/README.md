# Hardware kit: sparse Gaussian SYK (N = 8) circuits for the traversable-wormhole signal

Ready-to-run OpenQASM 3 circuits, already routed to physical qubits on the Nighthawk coupling map, for three selected sparse Gaussian SYK instances plus controls built with the same pipeline. An instance is one set of 24 quartet terms plus one Gaussian draw of their couplings (the code calls it a "listing").

| instance | folder | role |
|---|---|---|
| 19277_201 | `p325_tb_23_24_44300310_1429_19277_201` | recommended |
| 155_129 | `p324c_tb_23_24_44300310_1429_155_129` | backup, same term set |
| 1718_714 | `p330s2_tb_23_24_16500233_1414_1718_714` | borderline |
| Eq. 11 | `eq11_greedy` | positive control |
| typical 1_1, 3_3 | `typ_K24_1_24_1_1`, `typ_K24_1_24_3_3` | negative controls (typical K = 24 draws) |

Conventions follow arXiv:2604.10090 (beta = 3, |mu| = 12, FIX t0 = 1.8, SYM t0 = t1, base-2 entropies). All predictions come from a FakeNighthawk noise simulation; nothing has run on hardware.

## What to run

- For each folder and protocol P (FIX, SYM), run `P_mu-12.qasm` and `P_mu+12.qasm` as shipped: no re-transpiling (or initial_layout with optimization_level 0). Each header lists the physical qubits, routing seed and gate counts.
- Measure two-qubit Pauli tomography (9 settings, XX to ZZ) on P and register qubit 1 at the readout time t1. Then I_PT = S_P + S_T - S_PT and dI = I_PT(mu = -12) - I_PT(mu = +12).
- Shots: target standard error 0.005 bits on dI per circuit pair, which puts the 0.027-bit FIX gap between 19277_201 and Eq. 11 at about 3.8 sigma. For 19277_201 FIX: 14,900 shots per setting, so 2 signs x 9 settings x 14,900 = 268,200 shots per protocol.
- How we sized the shots: we simulated naive tomography of a two-qubit Werner state with the same noisy I_PT, found SE 0.0136 bits at 2,000 shots per setting, and scaled as 1/sqrt(N). A direct check at 14,900 shots gave SE 0.0051 bits and bias -0.0003 bits.

## Predictions

2q is the CZ count (same for mu = -12 and +12). The shipped-circuit column is the one circuit pair in this folder, which is what the hardware runs; a single routing seed moves by up to about 8% from the median. Noisy dI is the median over routing-seed list A (a fixed list of transpiler seeds) ± half the spread across independent lists A, B and C.

| instance | 2q FIX / SYM | noiseless dI FIX / SYM | noisy dI FIX / SYM | shipped circuit noisy dI FIX / SYM | shots per setting FIX / SYM |
|---|---|---|---|---|---|
| 19277_201 | 216 / 225 | 0.33 / 0.28 | 0.12 ± 0.005 / 0.096 ± 0.0006 | 0.12 / 0.094 | 14900 / 18000 |
| 155_129 | 216 / 225 | 0.27 / 0.25 | 0.097 ± 0.004 / 0.087 ± 0.0007 | 0.10 / 0.085 | 14000 / 17700 |
| 1718_714 (borderline) | 208 / 250 | 0.25 / 0.26 | 0.092 ± 0.002 / 0.075 ± 0.0003 | 0.094 / 0.075 | 15100 / 13200 |
| Eq. 11 | 175 / 215 | 0.21 / 0.20 | 0.089 ± 0.0009 / 0.074 ± 0.001 | 0.083 / 0.074 | 15600 / 17100 |
| typical 1_1 | 307 / 471 | 0.064 / 0.064 | 0.017 / 0.0069 | - | 11800 / 7500 |
| typical 3_3 | 288 / 406 | 0.065 / 0.12 | 0.015 / 0.013 | - | 15200 / 9900 |

Physical qubits (measure P and register qubit 1; register qubit 2 for reference):

| instance | FIX: P / reg 1 / reg 2 | SYM: P / reg 1 / reg 2 |
|---|---|---|
| 19277_201 | $63 / $43 / $33 | $22 / $52 / $41 |
| 155_129 | $63 / $43 / $33 | $22 / $52 / $41 |
| 1718_714 | $50 / $44 / $42 | $29 / $38 / $39 |
| Eq. 11 | $56 / $53 / $55 | $23 / $33 / $34 |
| typical 1_1 | $56 / $54 / $34 | $22 / $54 / $52 |
| typical 3_3 | $61 / $53 / $42 | $22 / $33 / $44 |

- Predicted ordering: 19277_201 > 155_129 > Eq. 11 on FIX and SYM, on all three seed lists. 1718_714 beats Eq. 11 on list A only, hence borderline. Typical draws sit 5 to 10x below.
- The shipped circuit is a single seed (lower median of list A); one seed can move noisy dI by up to about 8%.
- The options share few terms (19277_201 and 155_129 share all 24 and differ only in couplings), so a problem tied to one term structure should not hit every option.

## Noise-scale stress

All gate and readout errors scaled by S (T1, T2 divided by S); noisy dI for 19277_201 / 155_129 / Eq. 11:

- FIX: x0.5 0.19 / 0.16 / 0.14; x1 0.12 / 0.097 / 0.089; x2 0.043 / 0.036 / 0.038.
- SYM: x0.5 0.16 / 0.15 / 0.12; x1 0.096 / 0.087 / 0.074; x2 0.034 / 0.031 / 0.028.
- 155_129 FIX falls to Eq. 11 at about x1.55. 19277_201 beats Eq. 11 at every scale on both protocols, so run it first on a noisier device.

## Jordan-Wigner map, and an S3 correction

- Everything ships in the S2 map of the GIST group (Byun, Lee and Kim, arXiv:2604.10090; Eq. S2 there) for both protocols. In S2, R's psi^1, psi^2 sit on register qubit 1 (logical 7, replacing T). Gate count changes with the map; I_PT does not.
- Correction (2026-10-23): an earlier version said our instances also beat Eq. 11 in S3. That held only against Eq. 11 in our greedy term order. Against Eq. 11 in the term orders GIST actually used (different for H_L and H_R), SYM in S3 is a tie: 19277_201 +8 / +5 / -2% on seed lists A / B / C, 155_129 -2 / -5 / -9%, with Eq. 11 at 241 to 245 2q and noisy dI 0.0577 / 0.0580 / 0.0593.
- In S2, against Eq. 11 in GIST's orders, 19277_201 keeps, under noise, 23 to 60% more dI on FIX (23 to 50% if the backward step reuses the term order, 42 to 60% if it reverses it; the paper does not say which) and 51 to 67% more on SYM (Eq. 11 noisy dI 0.081 to 0.090 with the same order, 0.076 to 0.081 reversed; SYM 0.057 to 0.063). Eq. 11 then needs 195 to 206 (FIX) and 254 to 257 (SYM) two-qubit gates, more than in our greedy order. Against it the shot plan's FIX gap is 2.9 to 4.3 sigma.
- So run the S2 circuits. The `S3/` subfolders (19277_201, 155_129, Eq. 11) ship as report only; S3 SYM carries no claim. Never compare across maps.

## Retention rule to test

noisy dI = dI0 x a x exp(-b x 2q), with dI0 the noiseless dI and 2q the pair's mean CZ count; FIX a = 1.088, b = 0.00540; SYM a = 1.209, b = 0.00561 (fit to 125 simulated rows per protocol, rms 0.003 bits). A device that scales all ratios equally changes a, not the ordering. Known miss: it over-predicts SYM on one term set outside this kit by 14%, so 2q count alone does not set retention.

## Extra files (report only)

- `t1_scan.json` (19277_201, 155_129, eq11_greedy): the pipeline rerun at shipped t1, ±0.1 and ±0.2, with 2q and noiseless and noisy dI per point; no QASM off-centre. Both instances beat Eq. 11 at every offset; FIX dI rises about 2% at +0.2, SYM peaks within 0.1 of shipped t1.
- `chaos_curves.json`: spectral form factor against the typical K = 24 band, two-point function G(t) on [0, 20] and size winding for the same three Hamiltonians; its `keys` field explains the layout.

## Caveats

- FakeNighthawk (qiskit-ibm-runtime 0.49.0 snapshot, density matrix, gate and readout noise) is a proxy for the Nighthawk device, not device data.
- With idle noise added, 155_129 FIX ties Eq. 11 (0.0642 vs 0.0644); 19277_201 holds.
- The wins use our greedy term order on both sides. With SYM order optimized on both sides it is a tie at best (Eq. 11 routes to 134 to 136 2q); 19277_201 FIX stays +18%.
- These are selected, not typical, Hamiltonians: about 1 in 1.5e13 of the screened bank for 19277_201.
- The Werner stand-in matches I_PT only: treat shot counts as an order-of-magnitude budget.

Each folder's `predictions.json` cites the source of every number.
