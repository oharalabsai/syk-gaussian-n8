# Hardware kit: selected sparse Gaussian SYK (N = 8) circuits for the traversable-wormhole signal

What it is: ready-to-run routed circuits (OpenQASM 3, physical qubits on the Nighthawk map) for the listings that passed the frozen noise-comparison check (3 so far, on 2 distinct term set(s)), plus controls built with the same pipeline: Eq. 11 (positive, the published sparse model) and two typical K = 24 draws (negative). Each listing beat Eq. 11's noisy dI in the FakeNighthawk simulation on both protocols on routing-seed list A (the frozen gate), in canonical term order (order-optimised SYM: tie at best); the listings of set 44300310_1429 also hold on held-out lists B and C, the one marked borderline does not (see its card). Primary listing: 19277_201 (largest worst-protocol margin over Eq. 11).

## What to run
- For each directory and protocol (FIX, SYM): run `<P>_mu-12.qasm` and `<P>_mu+12.qasm` as is (no re-transpile; or initial_layout + optimization_level 0). The header names the physical qubits, the routing seed and the counts.
- Measurement: two-qubit Pauli tomography (9 settings XX..ZZ) of P and register qubit 1 (physical qubits in the header) at the readout time t1. I_PT = S_P + S_T - S_PT (base 2), dI = I_PT(mu = -12) - I_PT(mu = +12).
- Shots: SE(dI)^2 = SE(I_-)^2 + SE(I_+)^2 with SE(I) ~ c / sqrt(N) at N shots per setting; c from a parametric bootstrap of the naive tomography (PSD clipped) on a two-qubit Werner state whose I_PT equals the predicted noisy I_PT (FakeNighthawk routing-seed list B rows). Target SE(dI) = 0.005 bits per circuit pair (the 19277_201 vs Eq. 11 FIX gap 0.027 bits is then about 3.8 sigma). Worked number (primary listing 19277_201, FIX): predicted noisy I_PT(mu -12) and I_PT(mu +12) give Werner weights p = 0.275 / 0.128; at N0 = 2000 shots per setting the bootstrap SE of dI is 0.0136 bits, so N = 2000 x (0.0136 / 0.005)^2 = 14900 shots per setting (checked: SE 0.0051, bias -0.0003 bits). Total per protocol: 2 signs x 9 settings x 14900 = 268,200 shots.

| listing | 2q FIX (mu -/+) | 2q SYM | noiseless dI FIX / SYM | FakeNighthawk noisy dI FIX / SYM (routing-seed list A median ± half the lists A-C range) | retention pred FIX / SYM | shots per setting FIX / SYM |
|---|---|---|---|---|---|---|
| 155_129 (listing) | 216/216 | 225/225 | 0.27 / 0.25 | 0.097 ± 0.004 / 0.087 ± 0.0007 | 0.093 / 0.085 | 14000 / 17700 |
| 19277_201 (listing) | 216/216 | 225/225 | 0.33 / 0.28 | 0.12 ± 0.005 / 0.096 ± 0.0006 | 0.11 / 0.095 | 14900 / 18000 |
| 1718_714 (borderline) (listing) | 208/208 | 250/250 | 0.25 / 0.26 | 0.092 ± 0.002 / 0.075 ± 0.0003 | 0.089 / 0.078 | 15100 / 13200 |
| Eq. 11 (control+) | 175/175 | 215/215 | 0.21 / 0.20 | 0.089 ± 0.0009 / 0.074 ± 0.001 | 0.090 / 0.073 | 15600 / 17100 |
| typical 1_1 (control-) | 307/307 | 471/471 | 0.064 / 0.064 | 0.017 / 0.0069 | 0.013 / 0.0055 | 11800 / 7500 |
| typical 3_3 (control-) | 288/288 | 406/406 | 0.065 / 0.12 | 0.015 / 0.013 | 0.015 / 0.015 | 15200 / 9900 |

Shipped circuit vs shipped circuit (one routing seed each, the lower-median seed; what the hardware will run): 155_129 FIX 0.10 / SYM 0.085; 19277_201 FIX 0.12 / SYM 0.094; 1718_714 (borderline) FIX 0.094 / SYM 0.075; Eq. 11 FIX 0.083 / SYM 0.074 (FakeNighthawk gate_ro). The routing-seed list A medians above are the gated statistic; a single seed moves by up to ~8%.

## Noise-scale stress (FakeNighthawk errors x0.5 / x1 / x2)
- Every FakeNighthawk gate and readout error scaled x S and T1 / T2 divided by S (routing unchanged); median noisy dI, gate_ro, routing-seed list A central five, listing 19277_201 / listing 155_129 / Eq. 11. S = 1 reproduces this kit to 7e-16.
- FIX: x0.5 0.19 / 0.16 / 0.14; x1 0.12 / 0.097 / 0.089; x2 0.043 / 0.036 / 0.038 (155_129 falls 7% BELOW Eq. 11; 19277_201 stays 1.1x above).
- SYM: x0.5 0.16 / 0.15 / 0.12; x1 0.096 / 0.087 / 0.074; x2 0.034 / 0.031 / 0.028 (both above Eq. 11).
- FIX crossover (XCHK C, jobs r120 family xchkc): 155_129 vs Eq. 11 is +5% at x1.25, +1% (tie) at x1.5, -3% at x1.75; it falls to Eq. 11 at about x1.55. 19277_201 stays 1.2x above Eq. 11 across x1.25-x1.75.
- Reading: if the device is noisier than FakeNighthawk, run listing 19277_201 first; it beats Eq. 11 at every scale on both protocols. Jobs j-bec542824fab004f, j-3fa46a8a7c4104a2, j-cb4f4accddde3a48, j-d8065e07a0b0779e, j-939c4e1f0dfb8f0c, j-3a33cedc8d565913.

## Encoding (JW map)
- All shipped QASM and every hamiltonian.json use GIST's S2 Jordan-Wigner map (Eq. S2 of arXiv:2604.10090) for BOTH protocols, FIX and SYM (one map, chosen on the baselines only, applied identically to Eq. 11, typical draws and candidates). GIST's paper ran SYM in S3.
- I_PT is map invariant: S2 vs S3 agree to 1e-14, so the physics is the same; only the gate count changes.
- In S2, Majoranas psi^1, psi^2 of R sit on register qubit 1 (logical qubit 7; it replaces T) and those of L on register qubit 2 (logical 6). Tomography = P and register qubit 1; physical qubits per circuit in the table below.
- Run the shipped circuits as is; recompiling (another map, layout or seed) changes the 2q count and the prediction.
- S3 costs on our pipeline (report only, routing-seed list A lower median): set 44300310_1429 (listing 155_129) FIX 238 / SYM 301 vs Eq. 11 222 / 274 (same pipeline, j-8da943a7aaca0553). S3 noisy comparison: PASS, report only: S3 on BOTH sides, FakeNighthawk gate_ro, routing-seed list A central-five median (the SYM listing medians use 4 of the 5 central seeds: two shards hung and were cancelled): Eq. 11 FIX 0.070 (222 2q) / SYM 0.049 (274); 155_129 0.083 (238, +19%) / 0.057 (301, +15%); 19277_201 0.091 (238, +30%) / 0.063 (301, +28%); every listing seed beats every Eq. 11 seed. Frozen certification is S2.
- If you prefer S3: both listings of set 44300310_1429 and the Eq. 11 control have an S3/ subfolder (not the borderline listing) with the 4 S3 circuits (lower-median routing-seed list A seed, run as is), the S3 VQA angles and the two physical qubits to measure per circuit (P = logical 9; T is replaced by register qubit 1 = logical 7, the same qubit as in S2 since S3 only reorders pairs 3-8). Compare S3 with S3 and S2 with S2, never across maps.

| option | protocol | P | register qubit 1 (R psi^1,2; measured) | register qubit 2 (L psi^1,2) |
|---|---|---|---|---|
| 155_129 | FIX | $63 | $43 | $33 |
| 155_129 | SYM | $22 | $52 | $41 |
| 19277_201 | FIX | $63 | $43 | $33 |
| 19277_201 | SYM | $22 | $52 | $41 |
| 1718_714 (borderline) | FIX | $50 | $44 | $42 |
| 1718_714 (borderline) | SYM | $29 | $38 | $39 |
| Eq. 11 | FIX | $56 | $53 | $55 |
| Eq. 11 | SYM | $23 | $33 | $34 |
| typical 1_1 | FIX | $56 | $54 | $34 |
| typical 1_1 | SYM | $22 | $54 | $52 |
| typical 3_3 | FIX | $61 | $53 | $42 |
| typical 3_3 | SYM | $22 | $33 | $44 |

## Term overlap between options (independent failure)
Jaccard J = |A n B| / |A u B| on the quartet term sets, as labelled / max over all 8! relabelings of the 8 Majoranas (relabel-aware). Low values = the options share few terms, so a defect tied to one term structure should not hit both. Listings on the same set share every term (J = 1) and differ only in couplings J_T.

| | set d6a76272 (155_129, 19277_201) | set 2f2d8f4f (1718_714 (borderline)) | Eq. 11 | typical 1_1 | typical 3_3 |
|---|---|---|---|---|---|
| set d6a76272 (155_129, 19277_201) | - | 0.20 / 0.55 | 0.21 / 0.36 | 0.23 / 0.45 | 0.20 / 0.45 |
| set 2f2d8f4f (1718_714 (borderline)) | 0.20 / 0.55 | - | 0.06 / 0.36 | 0.20 / 0.50 | 0.23 / 0.50 |
| Eq. 11 | 0.21 / 0.36 | 0.06 / 0.36 | - | 0.13 / 0.42 | 0.13 / 0.36 |
| typical 1_1 | 0.23 / 0.45 | 0.20 / 0.50 | 0.13 / 0.42 | - | 0.20 / 0.55 |
| typical 3_3 | 0.20 / 0.45 | 0.23 / 0.50 | 0.13 / 0.36 | 0.20 / 0.55 | - |

## Predicted ordering (canonical greedy term order, both sides)
- 19277_201 > 155_129 > Eq. 11 on FIX and on SYM, on all three routing-seed sets:
  - routing-seed list A: FIX: 19277_201 0.12 > 155_129 0.097 > 1718_714 (borderline) 0.092 > Eq. 11 0.089; SYM: 19277_201 0.096 > 155_129 0.087 > 1718_714 (borderline) 0.075 > Eq. 11 0.074
  - routing-seed list B: FIX: 19277_201 0.11 > 155_129 0.094 > Eq. 11 0.089 > 1718_714 (borderline) 0.088; SYM: 19277_201 0.095 > 155_129 0.086 > Eq. 11 0.075 > 1718_714 (borderline) 0.074
  - routing-seed list C: FIX: 19277_201 0.12 > 155_129 0.10 > Eq. 11 0.090 > 1718_714 (borderline) 0.087; SYM: 19277_201 0.096 > 155_129 0.088 > Eq. 11 0.077 > 1718_714 (borderline) 0.074
- Typical K24 draws sit 5-10x below (median over 20 draws gate_ro FIX 0.014 / SYM 0.010).

## The hardware test: retention rule
noisy dI = dI0 x a exp(-b x 2q), dI0 = the circuit's noiseless dI, 2q = mean cz count of the mu -/+ pair; FIX a 1.088 b 0.00540, SYM a 1.209 b 0.00561 (least squares over 125 FakeNighthawk rows per protocol, rms 0.003 bits). Test: measured dI / (dI0 x a exp(-b 2q)) for every circuit pair here. The rule predicts the order above; a device that keeps the order but scales all ratios equally changes a, not the claim. Known miss: it over-predicts SYM on set 8400088_2290 (not in this kit) by 14% (retention 0.278 at 239 2q vs Eq. 11 0.370 at 215), so retention is not 2q count alone.

## t1 scan (report only)
What it is: `t1_scan.json` in 19277_201, 155_129 and eq11_greedy re-runs the shipped pipeline at 5 readout times, the shipped t1 and t1 +/- 0.1, +/- 0.2 (FIX: t0 stays 1.8; SYM: t0 = t1 moves with it). Each circuit is rebuilt from `hamiltonian.json` and routed at the shipped seed. No QASM ships for the off-centre points, and no verdict or prediction above changes.
How to read it: one row per t1 with `t0`, `2q`, `dI_noiseless` (model, fitted TFD) and `dI_gate_ro` (FakeNighthawk gate + readout, the `uv run noise` function); `shipped: true` marks the shipped circuit, whose values equal this kit's predictions (gate_ro to 1e-15).
- Both instances beat Eq. 11 at every matched offset, on both protocols (noisy dI ratio at least x1.44 / x1.23 FIX and x1.22 / x1.14 SYM for 19277_201 / 155_129).
- FIX: dI keeps rising slightly past the shipped t1 (the shipped t1 is the mu < 0 I_PT argmax, not the dI argmax): +0.2 gives 0.1199 -> 0.1217 (19277_201), 0.1021 -> 0.1043 (155_129), 0.0829 -> 0.0844 (Eq. 11), about +2% for all three.
- SYM: noisy dI peaks within 0.1 of the shipped t1: 0.0945 at 2.25 vs 0.0938 shipped (19277_201, which routes to 223 instead of 225 2q there), 0.0863 at 2.25 vs 0.0855 (155_129), 0.0750 at 2.15 vs 0.0742 (Eq. 11).
- 2q is flat across the scan (instances FIX 216 / SYM 225, Eq. 11 175 / 215) except that one point; gate_ro / noiseless stays within 0.34-0.40.

`chaos_curves.json` (report only) holds the chaos diagnostics of the same three Hamiltonians: the 18-window time-averaged SFF with the K = 24 typical band and the DENSE mean, the two-point function G(t) on [0, 20] and the size winding at t = 2.21 and at each shipped SYM time; its `keys` field says how to read it.

## Honest caveats
- FakeNighthawk (qiskit-ibm-runtime 0.49.0 snapshot, density matrix, gate + readout) is a proxy, not device data; its map is assumed identical to ibm_phoenix. Idle-noise variant: 155_129 FIX ties Eq. 11 (0.0642 vs 0.0644); 19277_201 robust.
- The win is in the canonical greedy term order on both sides. With order-optimised SYM on both sides it is a tie at best (Eq. 11 routes to 134-136 2q); 19277_201 FIX still +18% there.
- Selection rarity (frozen verifier selection_1_in_M, row.json one_in): 155_129 about 1 in 1.8e+13, 19277_201 about 1 in 1.5e+13, 1718_714 (borderline) about 1 in 1.2e+13 within the screened bank. These are selected, not typical, Hamiltonians.
- 2 distinct term set(s) so far (2f2d8f4f, d6a76272); 3 listings. More sets may be added.
- Shot estimate uses a Werner stand-in matched in I_PT only; the real reduced states have local Bloch vectors, so treat it as an order-of-magnitude budget; at the planned shots the dI estimator's bias is small and negative on every row (-0.0006 to -0.00005 bits, under 0.15 SE; kit.json shots.*.bias).

Sources: each directory's predictions.json carries per-number cites.
