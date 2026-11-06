# CE653: RC beam reliability — four standalone MATLAB scripts

## Run

Open any of the following files in MATLAB and run the **whole script**:

- `beam_mvfosm.m`: mean-value first-order second-moment method.
- `beam_hasofer_lind.m`: Hasofer–Lind/FORM with nonlinear marginal transformations and damped HLRF iterations.
- `beam_monte_carlo.m`: independent crude Monte Carlo.
- `beam_latin_hypercube.m`: independently randomized Latin hypercube replicates.

Each script is self-contained and leaves a `results` structure in the workspace. They clear the workspace at the start, so save a result before running another script, for example `save('mvfosm_results.mat','results')`. They do not automatically write output files. MATLAB R2016b or later supports the local functions used here. No Statistics, Optimization or Symbolic Math Toolbox is required.

Edit the local `beam_parameters()` function near the end of each script to change inputs. Keep its four copies consistent. Simulation sample sizes and seeds are near the top of the relevant script. FORM also prints iteration history, the design point, and squared direction cosines. MVFOSM prints first-order variance contributions.

## Beam and performance function

The model uses the **original supplied 6 m example**, not the separate 5 m project baseline: nominal width 250 mm, overall depth 450 mm, clear cover 30 mm, three 20 mm tension bars, 8 mm stirrups, M25 concrete and Fe500 steel. Span is fixed at 6 m.

With lengths in mm (except span L in m), strengths in MPa, line loads in kN/m and RC unit weight gamma in kN/m^3:

```
As    = 3*pi*phi_l^2/4                         [mm^2]
d     = D - c - phi_st - phi_l/2               [mm]
wSelf = gamma_RC*b*D*1e-6                     [kN/m]
R     = theta_R*F*As*(d - 0.77*F*As/(b*C))*1e-6 [kN m]
S     = (G + Q + wSelf)*L^2/8                  [kN m]
g     = R - S                                 [kN m]
Failure: g <= 0
```

`G` is superimposed dead load, excluding beam self-weight. `gamma_RC` is effective reinforced-concrete unit weight, already including reinforcement; steel weight is not added again. The final 450 mm depth is used for self-weight, rather than the 600 mm preliminary trial depth in the design slides. Nominal self-weight is therefore 2.8125 kN/m.

All 11 basic variables are assumed independent. One common random longitudinal diameter represents the three tension bars; this is a common-batch simplification, not three independent diameter realizations. The same diameter sample determines both As and d. The same b and D samples enter both resistance and self-weight, so those derived quantities are dependent automatically. Cover c is clear cover to the outside of the stirrup.

There are no additional random variables for As, d or self-weight. There are no design partial factors such as 1.5 on loads or 0.87 on steel strength in this physical resistance-demand model.

## Input distributions

Mean and SD below are arithmetic moments in the stated units; SD is not variance. The code computes the unrounded values where a formula is given.

| Variable | Meaning / unit | Distribution | Mean | SD | Basis |
|---|---|---|---:|---:|---|
| C | Concrete cube strength, MPa | Lognormal | 32.309393 | 4.846409 | Lognormal/COV 0.15 assumption; mean derived from 5th percentile 25 MPa, with [5] providing a transferable distribution precedent |
| F | Steel yield strength, MPa | Normal | 560 | 30 | [1], selected high-standard-production model, mean 500+2(30) |
| b | Width, mm | Normal | 250.75 | 5.50 | [2], mean 1.003 times nominal, SD 4+0.006 times nominal |
| D | Overall depth, mm | Normal | 451.35 | 6.70 | [2], same dimensional model |
| c | Clear cover, mm | Normal | 30 | 3 | [3], adopted COV 0.10 |
| phi_l | Longitudinal diameter, mm | Normal | 20 | 0.20 | [4], adopted COV 0.01 |
| phi_st | Stirrup diameter, mm | Normal | 8 | 0.08 | [4], same COV transferred to nominal 8 mm |
| gamma_RC | RC unit weight, kN/m^3 | Normal | 25 | 1.25 | [5], COV 0.05 |
| G | Superimposed dead line load, kN/m | Normal | 5.25 | 0.525 | [6], mean 1.05 times nominal 5; COV 0.10 |
| Q | Live line load, kN/m | Gumbel maximum | 9.195402 | 2.298851 | [6], mean 10/(1+0.35*0.25); COV 0.25 |
| theta_R | Resistance model factor | Normal | 1.01 | 0.0465 | [7], supplied excerpt, Eq. 7.40 |

These are adopted research models and transfer assumptions; they are not a fitted joint distribution from measurements of this particular beam. In particular, the load model in [6] is transferred from its original design context, not prescribed by IS 456. The live-load reference period must be matched to any subsequently chosen target beta; these scripts do not independently establish an annual or 50-year reliability claim.

For concrete, using ln(C) ~ Normal(lambda,zeta^2), with numerical C in MPa:

```
zeta  = sqrt(log(1+0.15^2)) = 0.1491663800
lambda = log(25) + 1.64485362695*zeta = 3.4642326861
mean(C) = exp(lambda + zeta^2/2)
```

The mean 32.309393 MPa is a derived value, not a mean quoted directly from [5]. Do not interpret M25 as the mean, or substitute a cylinder-strength distribution without conversion.

For Gumbel maxima, `F_Q(q)=exp(-exp(-(q-a)/s))`, with scale `s=SD*sqrt(6)/pi=1.792406440` and location `a=mean-0.5772156649*s=8.160797224`. The inverse CDF is `a-s*log(-log(p))`. This sign matters: MATLAB's extreme-value conventions must not be substituted without checking whether they describe minima or maxima.

## How beta is calculated

- MVFOSM: `beta = g(meanX)/sqrt(grad(g)*Cov(X)*grad(g)')`. Its reported normal-tail Pf is an approximation. It uses input moments, not the full nonnormal tail shapes.
- Hasofer–Lind/FORM: minimize `norm(u)` subject to `g(T(u))=0`, with independent standard-normal u and `T_i(u_i)=F_i^(-1)(Phi(u_i))`. Uses exact lognormal and Gumbel marginal transformations. Five starting points, a merit-function line search, and convergence/gradient checks are included. Multiple starts are not a proof of a global minimum. `alpha` points toward increasing safety, so for this safe origin `u_star ~= -beta*alpha`.
- Monte Carlo and LHS: estimate `Pf = number of failures / number of samples`, then `beta=-Phi^(-1)(Pf)`. They do not calculate beta as the simulated mean of g divided by its SD. The latter would be a different approximation.

Monte Carlo uses 20 million samples in batches of 100,000 and an exact binomial 95% interval. LHS uses 20 independent replicates of 1 million samples. It estimates uncertainty between replicate Pf estimates, with an approximate Student-t interval; a binomial interval would be inappropriate for dependent points within an LHS replicate. LHS requires more working memory (several hundred MB) than batched Monte Carlo.

Pf is only a few per million. Small trial runs may produce no failures: an infinite point estimate of beta is not evidence of infinite safety. Increase N or the LHS replicate sizes/count for greater precision. Twenty million samples still produce only roughly 55 failures on average, so neither simulation should be reported with excessive decimal precision. LHS is not guaranteed to reduce rare-event uncertainty in every run.

## Numerical checks and expected results

The four `.m` files passed a MISS_HIT MATLAB-language syntax parse. MATLAB/Octave was unavailable: **these files have not been executed in MATLAB**. Independent Python numerical checks verified the model gradient, marginal transformations, LHS stratification, and consistency of shared model/input functions. A separate constrained optimizer reproduced the FORM distance.

| Method | Beta in numerical cross-check | Details |
|---|---:|---|
| MVFOSM | 6.346975 | g(mean)=110.535379 kN m; linearized SD(g)=17.415442 kN m |
| FORM | 4.566737 | Independent optimization residual approximately zero |
| Monte Carlo | 4.556553 | 52 failures / 20 million; Pf=2.60e-6 |
| Replicated LHS | 4.529942 | 59 failures / 20 million; Pf=2.95e-6 |

Simulation values above are **Python validation runs**, not claimed MATLAB outputs. MATLAB's different random streams will give different counts. MC 95% Pf interval: [1.9418e-6,3.4096e-6]. LHS approximate 95% Pf interval: [2.3723e-6,3.5277e-6]. An additional conditional integration check, integrating out the live load analytically, gave Pf approximately 2.7757e-6 (standard error 8.84e-9), corresponding to beta approximately 4.5428. Both simulation intervals included that reference. FORM and simulation need not agree exactly because FORM approximates the failure surface by a tangent plane.

## Resistance model scope

The scripts preserve the agreed Ranganathan yielded-steel flexural equation. They diagnose samples violating a simple steel-yield strain compatibility condition, using xu=F*As/(0.546*C*b), eps_cu=0.0035 and Es=200000 MPa. This diagnostic does not replace the model, reject samples, or introduce another failure event. In the Python MC check, 0.01412% of samples exceeded that domain; none of its 52 failures did. This does not establish validity for every design or every tail realization. Substantial redesign, especially increasing reinforcement alone, requires reviewing the resistance branch and ductility assumptions. These beta values refer to this flexural model, not a complete assessment of shear, anchorage, serviceability or deterioration.

Unbounded normal/Gumbel distributions are retained as adopted. Nonphysical realizations are counted and warned about, not silently resampled or clipped. Only floating-point inverse-CDF endpoints are protected.

## References and links

[1] JCSS, Probabilistic Model Code, Part 3: Reinforcing Steel, Section 3.2.2. https://www.jcss-lc.org/publications/jcsspmc/rebar.pdf

[2] JCSS, Probabilistic Model Code, Part 3: Dimensions. https://jcss-pmc.github.io/PMC/part-03/dimensions.html

[3] Dan, W.; Yue, X.; Yu, M.; Li, T.; Zhang, J. (2023). Prediction and Global Sensitivity Analysis of Long-Term Deflections in Reinforced Concrete Flexural Structures Using Surrogate Models. Materials 16(13), 4671, Table 5. https://doi.org/10.3390/ma16134671

[4] Zhang, D.; Zeng, Y.; Fang, M.; Jin, W. (2019). Service Life Prediction of Precast Concrete Structures Exposed to Chloride Environment. Advances in Civil Engineering, Article 3216328, Table 2. The stirrup nominal diameter is adapted to the present beam. https://doi.org/10.1155/2019/3216328

[5] Miceli, E.; Gino, D.; Castaldo, P. (2025). Reliability assessment of robustness for reinforced concrete moment resisting frames. Developments in the Built Environment 21, 100639, Table 3. https://doi.org/10.1016/j.dibe.2025.100639

[6] Correia, R. S.; Bono, G. F. F.; Paliga, C. M. (2022). Evaluation of reliability of optimized reinforced concrete beams. Revista IBRACON de Estruturas e Materiais 15(4), e15409, Table 2. https://doi.org/10.1590/S1983-41952022000400009

[7] Ranganathan, R. Structural Reliability Analysis and Design, Example 7.9, Eqs. 7.39–7.40; resistance model and model-factor statistics taken from the user's supplied pages. Book record: https://books.google.com/books?id=7SPRPAAACAAJ

MATLAB documentation: https://www.mathworks.com/help/matlab/matlab_prog/local-functions-in-scripts.html ; https://www.mathworks.com/help/matlab/ref/erfcinv.html ; https://www.mathworks.com/help/matlab/ref/betaincinv.html
