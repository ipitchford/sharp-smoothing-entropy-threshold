# A sharp smoothing inequality in Gauss space and the exact threshold for local entropy concavity

Candidate proofs of the following statements. Not refereed by a journal; produced by an AI system
(Claude, Anthropic) on 10 November 2026 and revised on 11 November 2026 after an internal AI-assisted
review (not peer review). Appendix B of the paper says how the work was produced, checked and revised,
and Appendix C is the priority report.

In plain terms: mix two independent copies of a symmetric random quantity with a slider t and add
Gaussian noise. When the inputs are very close to Gaussian, the entropy of the result is a concave
function of t exactly when the inputs carry at most four sevenths of the output variance. With three or
more copies it is not concave on the whole range of weights at any signal fraction. The effects are
tiny (the proved neighbourhood is 1e-8 relative to the Gaussian density), so the interest is structural.

- Theorem A: for bounded g on the real line and Ornstein-Uhlenbeck levels r < s < 1, the sharp exponent
  E(r,s) in ||T_sqrt(s) g||^2 <= C ||g||_inf^(2-E) ||T_sqrt(r) g||^E, including the critical level
  s = b(r), where the L^2(dx) form is Luo's Theorem F.1 and the supremum-norm form was known only below it.
  Remark 1.1 compares E with the exponent 2 log s / log r that holds for all square-integrable g.
- Theorem B: a weighted form at the critical level.
- Theorem C: two copies, local entropy concavity at and below the signal fraction 4/7, with one radius.
  New at 4/7; the scheme of the proof is Schubr's. Proposition 6.6: the radius 1e-8 is admissible.
- Theorem D: a form of Luo's theorem (arXiv:2609.11418v2, Theorem 1.2) that two-copy concavity fails near
  the Gaussian for every signal fraction in (4/7, 1]. The proof (Appendix D) is Luo's with two steps changed.
- Corollary: the local-concavity set of the Evidence Press release "The critical signal fraction for local
  entropy concavity" (DOI 10.5281/zenodo.23910629) is exactly (0, 4/7], in the relative supremum norm
  and in every relative C^m_b norm (Remark 6.9). This also answers Luo's Problem 10.2 negatively for
  every q <= 4/7.
- Theorem E: three or more copies. For each signal fraction in (0, 1) there are strongly log-concave
  symmetric densities arbitrarily close to the Gaussian whose entropy is not concave on the simplex of
  weights; the failure lies near a face of the simplex. The obstruction is that of Madiman, Nayar and
  Tkocz (2020); new are three summands, noise, and the limit criterion at a face. New in content for
  q <= 4/7, where two copies are locally concave.
- Corollary E': for three or more copies no signal fraction in (0, 1] gives local concavity on the whole
  simplex (Theorem E, and Theorem D on an edge of the simplex).
- Proposition F and Remark 5.3: the best constant in Theorem A at the critical level is not one; for
  small r it exceeds 1.48 (computed with +-1 sign-polynomial test functions).

## Contents

- `paper.pdf`, `paper.md`        the paper
- `figures/`                     the three figures of the paper (PNG)
- `code/ingredients.py`          closed-form ingredients at (r,s) = (2/7,4/7); finite maxima over even k <= 400; level sums
- `code/radius_bound.py`         Proposition 6.6: closed-form bounds for all k, every displayed rounded value, and the radius 1e-8
- `code/holder.py`               the exponent E(r,s): algebra, optimality along h_n, the inequality on test functions
- `code/baseline.py`             Remark 1.1: E(r,s) against the L^2 baseline 2 log s / log r
- `code/best_constant.py`        Remark 5.3: exact +-1 functions sign(P_n) and the critical ratio as r -> 0
- `code/adversarial.py`          Theorems A and B at (2/7,4/7) on 759 bounded functions, with gradient ascent
- `code/bilinear.py`             the quadratic term M_t: coefficient formula and the bounds of Lemma 6.4 for all t
- `code/entropy_curvature.py`    Theorem C end to end: F''(t) at q = 4/7 for admissible perturbations
- `code/positive_control.py`     midpoint curvature at q = 4/7, 0.7, 0.8 on h_N (N <= 80) from exact 90-digit input coefficients: detects the failure above 4/7
- `code/above_threshold.py`      Theorem D: every formula of Appendix D at 60 digits, and the entropy itself at q = 1, 0.9
- `code/three_copies.py`         Theorem E: the boundary criterion for three copies (`--fine` for finer numerical parameters)
- `code/compare_runs.py`         convergence table: standard against fine parameters
- `code/make_figures.py`         draws the figures from the outputs above
- `run_all.sh`                   runs everything into a fresh directory and compares with `results/`
- `results/*.json`, `results/*.log`  outputs of one complete run of `run_all.sh` (11 November 2026, macOS arm64)
- `ENVIRONMENT.txt`              software versions used, and the record of the replays
- `refs.bib`                     the bibliography, copied from Crossref, DataCite and arXiv records
- `AI_INDEX.md`                  machine-oriented map of claims, evidence and limits
- `REVIEW_SUMMARY.md`            objections raised in the internal AI-assisted review and how each was handled
- `RESEARCH_GATES.json`          contribution map, sources read, limitation triage
- `REPLAY_RECEIPT.md`            the clean replay made before this release was frozen
- `CITATION.cff`, `.zenodo.json` citation metadata; `LICENSE.md` licences; `.github/workflows/replay.yml` CI
- `MANIFEST.sha256`              checksums of every file

## Reproduce

Python 3 with numpy, scipy, mpmath and matplotlib (tested versions in `ENVIRONMENT.txt`). With uv:

    PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.13 --with numpy==2.5.3 --with scipy==1.18.1 \
        --with mpmath==1.4.1 --with matplotlib==3.11.2 sh run_all.sh fresh_output

This takes about five minutes on one core of the machine listed in `ENVIRONMENT.txt`. Each script prints
`ALL OK` and exits with status 0 only if every check passes; the checks are explicit tests, not `assert`
statements, so they also run under `python -O`. `run_all.sh` ends with `ALL REPRODUCED` when every
regenerated JSON file is identical to the one in `results/`; a second run in the environment of
`ENVIRONMENT.txt` did so. On other platforms or library versions last digits may differ. On 11 November
2026 the scripts of the first version were also replayed on that macOS machine; every check passed, and
the JSON files agreed with the Linux originals to rounding level. To rebuild the PDF:

    pandoc paper.md -o paper.pdf --pdf-engine=tectonic

The computations are diagnostics. Series are truncated, integrals are sums on finite grids, and maxima
are taken over finite sets; nothing is interval arithmetic and no computation is part of a proof.

## What is not established here

- The indicative radius 1.5e-3 in Section 8 is not certified, and no computation measures the true
  radius. What is proved is the radius 1e-8. The numerical range of entropy_curvature.py is a quadrature
  on [-9,9] and was not recomputed with a wider window; the positive control tests levels N <= 80 only.
- Three copies: concavity on a fixed interior region of the simplex, and for inputs that match the
  Gaussian fourth moment, is open.
- Remark 5.4 (higher dimensions) is a sketch, apart from the inclusion of the d-dimensional set in (0, 4/7];
  Remark 6.10 is a heuristic.
- Priority for Theorems A (supremum norm at the critical level; general exponent), B, C and for the new
  parts of Theorem E rests on the bounded search recorded in Appendix C. The search did not cover
  MathSciNet, zbMATH or citation indexes.
- No journal refereeing and no formal verification. The analytic steps that no computation can test
  are listed in Appendix B of the paper.

## How to cite

Anonymous (2026). *A sharp smoothing inequality in Gauss space and the exact threshold for local entropy concavity*
(version 0.1.0-candidate, unrefereed). Zenodo. https://doi.org/10.5281/zenodo.24144687 (this version; all versions:
https://doi.org/10.5281/zenodo.24144686). Source: https://github.com/ipitchford/sharp-smoothing-entropy-threshold
(tag v0.1.0-candidate). Release page: https://evidencepress.org/releases/sharp-smoothing-entropy-threshold/.
Machine-readable citation: `CITATION.cff`.

This release extends the Evidence Press release *The critical signal fraction for local entropy concavity*
(https://doi.org/10.5281/zenodo.23910629), which proved local concavity below 4/7 and left 4/7 open.

## Status and licences

Unrefereed candidate. The proofs are written analytic arguments; they were checked only by internal AI-assisted review
(see `REVIEW_SUMMARY.md`), which is not peer review, independent reproduction or formal verification. Text, figures and
data are CC0-1.0; code is MIT (`LICENSE.md`). Creator: Anonymous.
