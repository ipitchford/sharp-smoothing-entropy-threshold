# Review summary (public, rating-free)

Internal AI-assisted review of the first version (10 November 2026) and of its revision (11 November 2026). This is
not human peer review, independent reproduction or formal verification.

**Who reviewed.** First round: five reviewer roles run as separate AI agents (editor, methodology, domain specialist,
cross-disciplinary perspective, devil's advocate), one adversarial review by a second AI model, and a clean replay of
every shipped script on a second machine (macOS arm64). Second round, on the revision: a re-review of everything that
changed, a separate proof audit of the new and rewritten proofs, and a second adversarial audit by the second model.

## Outcome

No reviewer in either round found a mathematical error in any theorem, and none found a counterexample. The first-round
replay passed every check; its JSON outputs agreed with the original Linux outputs to rounding level. The requested
changes concerned scope, framing for readers, attribution, the evidential status of some computations, and several
rounded constants. One round of revision, with extensions, addressed all of them; the second round found no error in
any theorem and led to the further repairs listed at the end.

## Objections and how each was handled

| Objection | Disposition |
|---|---|
| For three or more copies and signal fraction above 4/7, failure already follows from Luo's two-copy theorem by restricting to an edge of the simplex; the new content of Theorem E was not stated plainly, and q = 1 was not covered | New Corollary E' with proof: for every k >= 3 no signal fraction in (0, 1] gives local concavity on the whole simplex (Theorem E for q < 1, Theorem D on an edge for q in (4/7, 1]). The text now says that the new content of Theorem E is q <= 4/7, where two copies are locally concave and three are not |
| The three-copy failures sit near codimension-two faces, in a window that shrinks with the size of the perturbation; "no threshold for three copies" overstated this | Subtitle, abstract, Section 1.3, Theorem E and Section 9.1 now say "on the whole simplex, including near its faces"; concavity on a fixed interior region is stated as open (Section 9.3) |
| The Theorem E counterexamples have nonzero fourth cumulant, so they say nothing about inputs that match Gaussian moments | Remark 7.4(9) says so; the moment-matched three-copy question is listed as open |
| The threshold might depend on how "close to the Gaussian" is measured | Remark 6.9 with proof: the local-concavity set is (0, 4/7] in every relative C^m_b norm as well. Remark 6.10 (labelled a heuristic) explains why neighbourhoods that damp high Hermite levels exponentially could move it, consistent with Abu-Ajamieh's fixed direction being concave up to 5/8. Theorem C is noted to answer Luo's Problem 10.2 negatively for every q <= 4/7 |
| The end-to-end numerical test could not fail, so it could not support the suggestion that the true radius is much larger than 1e-8 | The suggestion was withdrawn. A new positive control (`positive_control.py`) shows that the curvature computation detects failure above 4/7 on Luo's test functions and gives clearly negative values at 4/7; its limits (levels up to N = 80, one amplitude) are stated. The proved radius remains 1e-8 |
| Relation to Gaussian noise stability and to the elementary square-integrable baseline was missing | Remark 1.1: noise-stability form, two-line baseline exponent 2 log s / log r, a comparison table (`baseline.py`) and the contrast with the discrete cube |
| Whether the best constant in Theorem A equals one | Remark 5.3 and Proposition F: it does not for small r; sign-polynomial test functions give ratios up to about 1.49 as r -> 0 (`best_constant.py`). Whether the supremum over r is finite is listed as open |
| Higher dimensions | Remark 5.4 proves that the d-dimensional local-concavity set is contained in (0, 4/7]; the rest of that remark is labelled a sketch |
| Rounded values in the explicit-radius proof (Proposition 6.6) were not all rounded in the safe direction | Every rounded value that enters an upper bound is now rounded up, and `radius_bound.py` checks each displayed value |
| Several small slips and imprecisions: two index slips in Lemma 4.2; an integrability statement missing before a duality step in Lemma 4.3; a missing factor in the constant of Theorem A(d); the size of the gap in Appendix D at q = 1; the status of two numerical constants (11.85 and 7.99) | All corrected in the text; the two constants are now described as truncated estimates, the second for admissible functions only |
| Attribution: the results of Abu-Ajamieh (2026) were not engaged in detail; further results of Schubr (2026) and the full strength of Luo's theorem were not described | Abu-Ajamieh's endpoint identity, sign condition and other theorems are now stated and related to Theorem E (his quantity B equals sqrt(24) times K_4); Schubr's Theorems 1.3-1.6 are described and his equation credited for Lemma 6.1; Theorem D is labelled "a form of Luo's Theorem 1.2", with his stronger statement described |
| The paper was hard to read for non-specialists | New "In plain terms" section (the question, a communication-channel reading, the answer, the tool, the size of the effects, the role of symmetry), a results-at-a-glance table and a shorter title |
| Descriptions of several tests in Appendix A and Section 8 were imprecise; one script used `assert` | Ledger rows and Section 8 corrected to say exactly what was computed; all scripts use explicit checks and also run under `python -O` |
| References | Three works added from registry records (O'Donnell 2014; Bobkov, Chistyakov and Götze 2013; Nourdin, Peccati and Swan 2014), each cited where used |

## Second round (re-review, proof audit, second adversarial audit)

No error was found in any theorem. Further changes:

- The positive control was rewritten to use exact 90-digit input coefficients and a coefficient pairing for the
  curvature, because double-precision quadrature of the inputs had set a floor that invalidated some of its first
  values. The earlier values were withdrawn.
- The window of the end-to-end harness (`entropy_curvature.py`) is stated as a limitation.
- Remark 5.3 no longer asserts what the extremal functions are; the strictness of E > E_2 in Remark 1.1 is labelled as
  computed, not proved.
- Abu-Ajamieh's radius and Schubr's Theorem 1.6 are stated exactly; Section 8 lists each radius constant as an upper
  bound together with its evaluated value.
- Wording of the plain-language summary and the results table was tightened; remarks were renumbered.

## What the review did not do

No human specialist read the proofs, no proof-assistant formalisation was made, and the numerical work consists of
finite-precision diagnostics, not interval certificates. The priority search was bounded (Appendix C of the paper).
The analytic steps that no computation can test are listed in Appendix B of the paper.
