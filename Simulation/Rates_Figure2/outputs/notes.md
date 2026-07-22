# Phase Transition Run Notes

## Theory Mapping
- Corollary 1 uses slope `-(a-1)/(3a+1)` with polynomial eigenvalue decay `lambda_k ~ k^{-a}`.
- Here `decayR = 1.00` means score SDs decay like `k^{-decayR}`, so variances decay like `k^{-2*decayR}` and `a = 2.0`.
- The theoretical log-log slope is `-0.1429`.

## Final Settings
- Grid: `n = c(20, 30, 40, 50, 70, 100, 150)`, `N = c(3, 4, 5, 7, 10)`, `m = c(200)`.
- Replications: `200`; seed base: `20260323`; cores: `5`.
- Design: `random`; latent = `TRUE`; methodXi = `CE`.
- Stability knobs: `eps_u = 0.05`, `u_cap = 0.95`, `Z_cap = 8.0`, `t_grid_len = 51`.
- Runtime: 1098.2 seconds (18.30 minutes).
- Modular method file sourced from `/Users/hangzhou/Documents/research/Hans/random object/JRSSB_RR1/new_code/simulation_rate/codes/transport_fpca_modular.R`.

## Tail Slopes From Replication-Level IMSE
- Slopes use `n >= 70`.

| N | tail n values | slope | SE |
| --- | --- | --- | --- |
|  3 | 70,100,150 | -0.2090 | 0.0338 |
|  4 | 70,100,150 | -0.1734 | 0.0357 |
|  5 | 70,100,150 | -0.1189 | 0.0218 |
|  7 | 70,100,150 | -0.1315 | 0.0224 |
| 10 | 70,100,150 | -0.0939 | 0.0148 |

## Mean-Curve Slopes For Plotting
| N | all-n curve slope | tail-only curve slope |
| --- | --- | --- |
|  3 | -0.2589 | -0.2386 |
|  4 | -0.2437 | -0.1874 |
|  5 | -0.2226 | -0.1389 |
|  7 | -0.1684 | -0.1544 |
| 10 | -0.1372 | -0.1081 |

## Diagnostics
- IMSE was evaluated on the fixed grid `t_grid_out = seq(0, 1, length.out = 51)` for every setting.
- Mean generated sample variance across cells ranged from `0.0829` to `0.0840`.
- True transport values stayed in [`0.0000`, `1.0000`]; predicted transport values stayed in [`0.0000`, `1.0000`].
- Mean clamp / saturation diagnostics: `frac_abs_u_raw_gt1` in [`0.0003`, `0.0010`], `frac_u_grid_at_cap` in [`0.0000`, `0.0000`], `frac_pred_inner_boundary` in [`0.0000`, `0.0000`].

## Flattening Diagnosis
- This run writes versioned outputs so the earlier sparse-grid files remain available for comparison.
- The main dense-grid plot focuses on `N <= 10` so that small-N curvature is visible and the earlier `N = 20` flattening does not dominate the figure.

## Failure Check
- No failed fits.
