# AI index: a sharp smoothing inequality in Gauss space and the exact threshold for local entropy concavity (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) and the
paper [paper.pdf](paper.pdf) (source [paper.md](paper.md); its "In plain terms" section and results table come first).
DOI 10.5281/zenodo.24144687 (concept 10.5281/zenodo.24144686). Status: unrefereed candidate.

## Normalisation

X, Y (and X_1..X_k) i.i.d. with even density f, mean 0, variance 1; G standard Gaussian; signal fraction q in (0,1];
two copies: F_{q,f}(t) = h(sqrt(qt) X + sqrt(q(1-t)) Y + sqrt(1-q) G), natural logarithms; k copies: weights t in the
simplex, F(t) = h(sum_i sqrt(q t_i) X_i + sqrt(1-q) G). f = phi(1+g) with g admissible: real, even, bounded,
int g phi = int x^2 g phi = 0. T_rho g(x) = E g(rho x + sqrt(1-rho^2) G) (Ornstein-Uhlenbeck), N He_n = n He_n,
norms in L^2(phi). For 0<r<1: a_r = (1-r^2)^(-1/2), lambda_r = r/(1+sqrt(1-r^2)), b(r) = sqrt(2r(1+r)) - r,
tau_r(s) = (1-a_r(1-s))/(1+a_r(1-s)). C = { q : some delta>0 makes F_{q,phi(1+g)} concave on [0,1] for every
admissible g with ||g||_inf <= delta }; C_k the same for concavity on the whole k-simplex (C_2 = C).
Total SNR = q/(1-q); q = 4/7 is SNR 4/3.

## Claims and where they are proved

| ID | Claim | Location in [paper.md](paper.md) | Status / credit |
|---|---|---|---|
| A | Theorem A: for bounded g and 0<r<s<1, ‖T_sqrt(s) g‖^2 <= C ‖g‖_inf^(2-E) ‖T_sqrt(r) g‖^E with E(r,s) = 2 log tau_r(s) / log lambda_r; E decreases from 2 to 0 and E >= 1 iff s <= b(r); E optimal; at s = b(r) a logarithmic gain (1 + log(‖g‖_inf/‖T_sqrt(r) g‖))^(-1/4) | Section 1.3; proofs Sections 3-5 | new in bounded search (supremum norm at the critical level; general exponent). Critical level b(r) and the L^2(dx) linear case are Luo's Theorem F.1 |
| B | Theorem B: ‖N T_sqrt(s) g‖ ‖T_sqrt(s) g‖ <= C_r ‖g‖_inf ‖(N+1) T_sqrt(r) g‖ at s = b(r) | Section 1.3; Section 5 | new in bounded search |
| C | Theorem C: one delta_* > 0 such that for every q in (0,4/7], admissible ‖g‖_inf <= delta_*, t in [0,1]: F'' <= -(1/16) ‖N(T_sqrt(qt) g + T_sqrt(q(1-t)) g)‖^2; delta_* = 1e-8 admissible (Proposition 6.6) | Section 1.3; Section 6 | new at q = 4/7 and the uniform radius; proof scheme is Schubr's |
| D | Theorem D: for q in (4/7,1], arbitrarily small admissible smooth strongly log-concave g with F(t_0) > F(1/2) | Section 1.3; Appendix D | Luo's theorem (arXiv:2609.11418v2, Thm 1.2) and Luo's proof with two steps changed |
| Cor | C = (0,4/7]; also in every relative C^m_b norm (Remark 6.9); negative answer to Luo's Problem 10.2 for all q <= 4/7 | Section 1.3; Sections 6.7, 9.1 | follows from C and D |
| E | Theorem E: k >= 3, q in (0,1): arbitrarily small admissible smooth g, strongly log-concave, with F strictly convex on the segment (eps+tau, eps-tau, 1-2eps, 0, ...) for small eps; limit criterion F'' -> -4 hat h_4 K_4(h) at a face | Section 1.3; Section 7 | obstruction of Madiman-Nayar-Tkocz (2020); criterion is Abu-Ajamieh's sign condition read at a face; three summands, noise and the face limit new in bounded search |
| E' | Corollary E': C_k is empty for every k >= 3 | Section 1.3 | from E and D (edge of the simplex) |
| F | Proposition F / Remark 5.3: best constant in Theorem A at s = b(r) exceeds 1 for small r (computed limits up to 1.487) | Sections 1.3, 5 | proof of the limit formula; values computed |
| R1.1 | Remark 1.1: E >= E_2 = 2 log s / log r (square-integrable baseline) | Section 1.3 | proved; strictness computed on a grid only |
| R5.4 | Remark 5.4: d-dimensional local-concavity set contained in (0,4/7] | Section 5 | proved; rest of the remark a sketch |

Not claimed: concavity on a fixed interior region of the simplex for k >= 3; anything for moment-matched three-copy
inputs; the indicative radius 1.5e-3 (not certified); any global (non-local) result; asymmetric inputs (they fail for
every q, Luo Thm 1.1); Remark 6.10 (heuristic); priority beyond the bounded search of Appendix C.

## Evidence files

- Check scripts in [code/](code/), outputs of one complete run in [results/](results/), figures in [figures/](figures/).
  One command, [run_all.sh](run_all.sh), reruns everything into a fresh directory and compares each JSON byte for byte
  with [results/](results/). Software versions: [ENVIRONMENT.txt](ENVIRONMENT.txt). What each script tests and what it
  cannot show: Section 8 and Appendix A of the paper.
- Per claim: A, B -> [code/holder.py](code/holder.py), [code/adversarial.py](code/adversarial.py),
  [code/ingredients.py](code/ingredients.py); C -> [code/radius_bound.py](code/radius_bound.py),
  [code/bilinear.py](code/bilinear.py), [code/entropy_curvature.py](code/entropy_curvature.py),
  [code/positive_control.py](code/positive_control.py); D -> [code/above_threshold.py](code/above_threshold.py);
  E -> [code/three_copies.py](code/three_copies.py), [code/compare_runs.py](code/compare_runs.py);
  F -> [code/best_constant.py](code/best_constant.py); R1.1 -> [code/baseline.py](code/baseline.py).
- [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md): clean replay before freezing (normal and PYTHONOPTIMIZE=1, ALL REPRODUCED).
- [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md): internal AI-assisted review objections and dispositions.
- [RESEARCH_GATES.json](RESEARCH_GATES.json): contribution map, sources read, limitation triage.
- [refs.bib](refs.bib): bibliography from registry records. [MANIFEST.sha256](MANIFEST.sha256): SHA-256 of every file.
  [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json): citation metadata. [LICENSE.md](LICENSE.md): CC0-1.0 text
  and data, MIT code. CI: [.github/workflows/replay.yml](.github/workflows/replay.yml).

## Limits

No human peer review, independent human reproduction or formal verification. The computations are finite-precision
diagnostics (truncated series, grids, finite maxima), not interval certificates, and no computation is part of a proof.
The effects are minute (proved radius 1e-8 relative; largest computed failure signal about 6e-8 nats). The four 2026
sources the paper builds on (Luo, Schubr, Abu-Ajamieh, the earlier Evidence Press release) are unrefereed.
Predecessor release: https://doi.org/10.5281/zenodo.23910629 (proved (0,4/7) in C; left 4/7 open).
