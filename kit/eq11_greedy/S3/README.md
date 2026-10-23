# eq11_greedy: S3 variant (secondary; the frozen certification and the main circuits are S2)

Superseded 2026-10-23: against Eq. 11 in its authors' term orders (241-245 2q; noisy dI 0.0577 / 0.0580 / 0.0593 on routing-seed lists A / B / C), S3 SYM ties for 19277_201 (+8 / +5 / -2%) and loses for 155_129 (-2 / -5 / -9%). The PASS below is against Eq. 11 in our greedy order only; S3 FIX was not rechecked. Run the S2 circuits.

S3 noisy-sim comparison: PASS, report only: S3 on BOTH sides, FakeNighthawk gate_ro, routing-seed list A central-five median (the SYM listing medians use 4 of the 5 central seeds: two shards hung and were cancelled): Eq. 11 FIX 0.070 (222 2q) / SYM 0.049 (274); 155_129 0.083 (238, +19%) / 0.057 (301, +15%); 19277_201 0.091 (238, +30%) / 0.063 (301, +28%); every listing seed beats every Eq. 11 seed. Frozen certification is S2.

Same Hamiltonian, times and protocol as the parent folder; JW map S3 on both protocols, VQA refit in S3 (angles.json, same register-index convention as hamiltonian.json conventions.tfd). Run as is, no re-transpile.

| protocol | mu | 2q (cz) | depth | measure P | measure R psi^1,2 host |
|---|---|---|---|---|---|
| FIX | -12 | 222 | 414 | $35 (logical 9) | $53 (logical 7, register qubit 1) |
| FIX | +12 | 222 | 413 | $35 (logical 9) | $53 (logical 7, register qubit 1) |
| SYM | -12 | 274 | 564 | $24 (logical 9) | $42 (logical 7, register qubit 1) |
| SYM | +12 | 274 | 563 | $24 (logical 9) | $42 (logical 7, register qubit 1) |
