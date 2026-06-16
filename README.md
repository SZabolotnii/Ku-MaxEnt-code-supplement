# Generating-Element Maximum Entropy — Reproducibility Repository

Verification code for the paper:

> **Generating-Element Maximum Entropy for Non-Gaussian Uncertainty Evaluation**
> Serhii Zabolotnii, 2026.

**Repository:** <https://github.com/SZabolotnii/Ku-MaxEnt-code-supplement>
(`git clone https://github.com/SZabolotnii/Ku-MaxEnt-code-supplement.git`)

**One-line summary.** Moment-constrained maximum-entropy (MaxEnt) density
reconstruction is governed by the *choice of generating element* of the
underlying Kunchenko decomposition space — not by the solver — and matching
that element to the target's tail class (fractional-power / trigonometric /
logarithmic-rational, under one dual Newton solver) decides which densities are
representable and how well-conditioned the dual problem is.

---

## What this repository verifies

A single dual MaxEnt solver (damped Newton with backtracking line search and a
ridge-regularized dual Hessian) is exercised across four experiments:

1. **Heavy tails (standard Cauchy)** — fractional-power PATP restores solver
   feasibility where classical monomial MaxEnt is infeasible, and reconstructs
   the body.
2. **Bimodal Gaussian mixture** — the scan-selected PATP member beats the
   six-moment monomial baseline on accuracy *and* conditioning; trigonometric
   (T-MaxEnt) is best-conditioned.
3. **Matching principle across tail classes** (Cauchy, Student-*t*(2),
   symmetric α-stable(1.5), Gaussian) — the matched logarithmic element
   `log(1+(x/s)^2)` is the only element that recovers an algebraic tail index.
4. **M&V sampling-design optimization** — an analytical product-moment MaxEnt
   fitness makes a genetic-algorithm evaluation exactly deterministic and far
   faster than Monte Carlo.

All code is **base R** (no contributed packages). Every run is deterministic
under a fixed seed.

---

## One-to-one mapping: script → paper figure / table / claim

| Script (`R/`)                | Paper artifact | Reproduces |
|------------------------------|----------------|------------|
| `patp_maxent_simulation.R`   | Exp 1 & 2 core | Canonical `solve_maxent` solver + the all-odd Form-B basis; Cauchy α-sweep and mixture α-sweep console baselines that the ablations extend. |
| `ablation_cauchy_metrics.R`  | **Table 1** (`tab:cauchy`), **Fig 1** (`fig:cauchy`) discussion | Per-α Cauchy metrics: κ_H, dual potential Γ, body MSE on [-10,10], KS on [-45,45], Q75/90/95/99 errors. Confirms α=0 body MSE ≈ 2.46×10⁻³, α=0 KS ≈ 0.092. |
| `ablation_mixture.R`         | **Table 2** (`tab:mixture`), **Table 3** (`tab:ablation`) | Mixture baselines at 6 constraints (T-MaxEnt, PATP 0.9/1.0, monomial, Legendre); T-MaxEnt (p,S) frequency ablation with empirical-CF-amplitude rule; ridge sensitivity. |
| `ablation_ga_noise.R`        | **Table 6** (`tab:mv`) | Experiment 4: MaxEnt-vs-Monte-Carlo GA fitness — per-call speedup, deterministic-fitness assertion (sd = 0), constraint violations (3/10 MC vs 0/10 MaxEnt). |
| `extra_ablations.R`          | **§ Replication over seeds**; **Table 5** (`tab:designmap`) support | Linear-exponent-map control (does the PATP *quadratic* map matter?) + 20-seed replication of mixture and Cauchy headlines: MSE ratio mono/selected **8.5 ± 5.8**, α* median 0.7, Cauchy KS 0.068 ± 0.026, α=0 converges 18/20, four-moment monomial infeasible **19/20** (certificate m̂₄ > L⁴ = 6.25×10⁶; secondary m̂₂ > 2500 witness fires 7/20). |
| `opmm_alpha_cauchy.R`        | **Revision (reviewer Q1/Q6): oPMMα selection, §sec:selection / Table sel** | Variance-optimal generating-element selection V(α)=∇Tᵀ H⁻¹ Σ(α) H⁻¹ ∇T (PMM/delta method). Targeting the reported body-coverage functional picks the most-accurate member **17/18 seeds (body-MSE penalty 1.04×)** vs the dual-potential heuristic **0/18 (2.17×)** — the principled replacement for Γ-minimization on heavy tails. |
| `opmm_alpha_mixture.R`       | **Revision (reviewer Q1/Q6): selection control, §sec:selection / Table sel** | Mixture (light-tailed) control. Held-out CV log-score is the most reliable shape selector (**6/20**, penalty 2.29×) vs Γ (4/20, 2.57×); a deliberately saturated V target P(\|X\|<3)≈1 degenerates (1/20, 14.9×), showing oPMMα needs an *informative* reported functional. |
| `free_exponent_baseline.R`   | Design justification (Exp 2 text / FM-MEM comparison) | FM-MEM-style Nelder-Mead search over free exponents vs the 11-point 1-D α scan — shows the scan is competitive at a fraction of the solver calls. |
| `log_rational_element.R`     | **Table 4** (`tab:genelement`), **Fig 3** (`fig:genelement`) | Experiment 3: 6 elements × 4 targets KS table; Cauchy tail-slope recovery; `fig_genelement_cauchy.pdf`. Confirms LogRat-1 Cauchy KS ≈ 0.015 vs power-PATP ≈ 0.087, λ̂ ≈ −0.94 → slope −1.88. |
| `log_rational_unbounded.R`   | **Revision (reviewer Q3/W3): §sec:genelement** | Closed-form matched fit on ℝ (no truncation): Z=s·B(½,−λ−½), moment ψ(−λ)−ψ(−λ−½). On Cauchy λ̂=−0.989, slope −1.98 (vs truncated −1.89), KS over ℝ **0.0035** (vs 0.0147); removes quantile saturation — truncated q99=23.5 (−26%) vs unbounded 34.5 (+8.5%, true 31.8). |
| `sensitivity_L_N.R`          | **Revision (reviewer Q2/W6): Table sens, §sec:sensitivity** | Truncation L × sample-size N grid (10 seeds/cell) for Cauchy and mixture: convergence, conditioning, body MSE, KS, quantile error, monomial-infeasibility rate. Shows body/shape/α* robust to L,N; only extreme-tail quantiles are L-sensitive (q99 −55%@L=20 → −35%@L=200); monomial infeasibility tracks m̂₄>L⁴. |
| `fm_mem_ga.R`                | **Revision (reviewer Q4/W5): §sec:exp2** | FM-MEM head-to-head: GA joint exponent+multiplier optimization (Zhang 2020 style) vs the 11-point PATP α-scan on the mixture (5 seeds). Realizable Γ-selected: FM-GA MSE 1.0e-4 @ 209 solves ≈ tied with scan 9.7e-5 @ 10 solves (**21× cost, no gain**). Oracle-MSE GA 5.9e-5 vs scan-oracle 7.6e-5 (only 22%, unrealizable). Confirms the scan captures the realizable benefit of free exponents at 1/21 the cost. |
| `mv_qmc_baseline.R`          | **Revision (reviewer Q5/W7): §sec:exp3** | Stronger M&V baseline: randomized quasi-MC (Halton+Cranley-Patterson) fitness + per-seed final designs for MC/RQMC/MaxEnt arms. RQMC cuts probe sd 12.6→4.3 (2.9×) and incurs 0/10 violations (vs plain MC 3/10) — so the analytical evaluator's durable edge is exact determinism (sd=0), true cumulants, and lowest cost sd (6.6 vs 10.2 vs 17.4), NOT a lower violation rate or speed over QMC. ⚠ slow (~few min: true-q05 checks at N=2e5). Env `MV_GA_SEEDS=0` for probe only. |
| `tmaxent_auto_freq.R`        | **Revision (reviewer Q6): §sec:tmaxent + §sec:exp2** | Automated T-MaxEnt frequency rule replacing the fixed 0.05: N-aware 3σ admissibility (\|ψ̂(jp)\|≥3/√N) + held-out log-score selection. Over 10 seeds it lands in the accurate p∈{0.5,0.7} region (mean PDF-MSE penalty 1.7× vs oracle) vs 3.3× for the fixed (0.5,3) default, and rejects noise-floor configs (1.0,4)/(1.0,5) with no hand-set threshold. |
| `parity_matched_patp.R`      | **Fig 1** (`fig:cauchy`), **Fig 2** (`fig:mixture`); rows of **Table 1** & **Table 2** | Parity-matched PATP (PM-PATP) Cauchy + mixture sweeps and the **final paper figures** `fig_cauchy_logpdf.pdf` and `fig_mixture_fits.pdf`. Confirms Cauchy α*=0.3 (Γ=2.395), mixture α*=0.7 (MSE 8.0×10⁻⁵, ~7× over monomial). |

Output figures (written to `outputs/`):
`fig_cauchy_logpdf.pdf`, `fig_mixture_fits.pdf`, `fig_genelement_cauchy.pdf`.

> Run order note: `parity_matched_patp.R` runs **last** so that the
> PM-PATP/T-MaxEnt variants of `fig_cauchy_logpdf.pdf` and `fig_mixture_fits.pdf`
> (the ones used in the manuscript) supersede the intermediate versions written
> by the two ablation scripts.

---

## How to run

Requirements: **R ≥ 4.5** (developed and verified on **R 4.5.3**), base R only —
no package installation needed.

From the repository root:

```sh
Rscript run_all.R
```

This sources every experiment in order in an isolated environment, writes all
figures to `outputs/`, and writes a per-script transcript to
`outputs/logs/<script>.log` plus a timing/status summary to
`outputs/run_all_summary.txt`.

To run a single experiment directly (figures still land in `outputs/`):

```sh
GENELEMENT_OUT_DIR=outputs Rscript R/log_rational_element.R
```

(Figure-writing scripts honor the `GENELEMENT_OUT_DIR` environment variable;
it defaults to `outputs` relative to the working directory.)

### Reproducibility settings

- **Fixed seed `20260612`** for all single-realization headline experiments.
- Seed-replication studies use seeds `1..20` (mixture/Cauchy in
  `extra_ablations.R`) and `1001..1010` (GA runs in `ablation_ga_noise.R`).
- The solver is deterministic; no randomness beyond the seeded data draws.

### Expected runtime

~3 minutes total on a single core of an Apple M-series CPU. Almost all of it is
`ablation_ga_noise.R` (~3 min: nested Monte-Carlo GA runs); every other script
finishes in well under a second.

---

## What to look for (key diagnostics)

| Quantity | Expected | Where |
|---|---|---|
| Cauchy body MSE at α=0 | ≈ 2.46×10⁻³ | `ablation_cauchy_metrics.log`, `parity_matched_patp.log` |
| Cauchy KS at α=0 | ≈ 0.092 (single seed); 0.068 ± 0.026 (20 seeds) | `parity_matched_patp.log`, `extra_ablations.log` |
| Cauchy potential-selected α* | 0.3 (Γ = 2.395) | `parity_matched_patp.log` |
| LogRat-1 Cauchy KS vs power-PATP | **0.015 vs 0.087** | `log_rational_element.log` |
| LogRat-1 fitted λ̂ → tail slope 2λ | λ̂ ≈ −0.94 → slope ≈ −1.88 (true Cauchy −2) | `log_rational_element.log` |
| `E_Cauchy[log(1+X²)]` (population identity) | 2 ln 2 ≈ 1.386 (empirical at N=1000, seed: ≈ 1.42) | LogRat-1 target moment |
| α-stable(1.5) true tail slope | −2.5 | `log_rational_element.log` target table |
| Mixture potential-selected α* | 0.7 (MSE 8.0×10⁻⁵, κ_H 2.4×10⁴) | `parity_matched_patp.log`, `extra_ablations.log` |
| Mixture MSE ratio monomial/selected (20 seeds) | **8.5 ± 5.8** (favorable all 20) | `extra_ablations.log` |
| GA MaxEnt fitness determinism | sd = 0 (assertion passes) | `ablation_ga_noise.log` |
| GA noise-induced constraint violations | 3/10 (MC) vs 0/10 (MaxEnt) | `ablation_ga_noise.log` |

**Note on the α-stable reference density.** `log_rational_element.R` computes
the symmetric α-stable PDF by **trapezoidal** characteristic-function inversion
(`du = 0.005`, `u = seq(0, 40, by = du)`, endpoints half-weighted). This is
intentional and must be preserved: an earlier rectangle-rule version injected a
`du/(2π)` DC offset and was a bug.

---

## Repository layout

```
.
├── R/                          verification scripts (base R)
│   ├── patp_maxent_simulation.R     core solver + Exp 1/2
│   ├── ablation_cauchy_metrics.R    Exp 1 extended metrics
│   ├── ablation_mixture.R           Exp 2 baselines + freq ablation
│   ├── ablation_ga_noise.R          Exp 4 GA (MaxEnt vs MC)
│   ├── extra_ablations.R            linear-map control + 20-seed replication
│   ├── free_exponent_baseline.R     FM-MEM free-exponent baseline
│   ├── log_rational_element.R       Exp 3 generating-element comparison
│   └── parity_matched_patp.R        PM-PATP final figures (run last)
├── run_all.R                   master driver
├── outputs/                    figures + logs (gitignored; created at run time)
├── README.md
├── LICENSE                     MIT (code only)
├── CITATION.cff
└── SESSION_INFO.md             captured sessionInfo()
```

### A note on the duplicated solver

`solve_maxent` is **deliberately copied verbatim** into each script rather than
factored into a shared file: each experiment is meant to be a self-contained,
auditable artifact, and there are small, intentional variant settings
(`ablation_mixture.R` exposes the ridge as an argument for its sensitivity
sub-study; `log_rational_element.R` raises `max_iter` to 200). Keeping them
separate avoids silently changing any experiment's numerical behavior.

## License

Code: MIT (see `LICENSE`). The manuscript text and its published figures are
**not** covered by this license.
