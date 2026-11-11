#!/usr/bin/env python3
"""Direct numerical tests of the two endpoint inequalities at (r, q) = (2/7, 4/7).

For a bounded function g sampled on a grid, compute its Hermite coefficients in L^2(phi),
G_n = <g, e_n>_phi (n < NMAX), and the ratios

  rho0(g) = ||T_sqrt(q) g||^2 / ( ||g||_inf ||T_sqrt(r) g|| )                      (plain endpoint lemma)
  rho1(g) = ||N T_sqrt(q) g|| ||T_sqrt(q) g|| / ( ||g||_inf ||(N+1) T_sqrt(r) g|| )  (weighted version)

The proof gives rho0 <= J_end and rho1 <= C_E1 (constants evaluated in ingredients.py).
The script also tests the per-level coefficient bounds
  |c_k| <= S_k(r)^(-1/2) ||T_sqrt(r) g||            and
  |c_k| <= c_L(k) (k+1)^(-1) S_k(r)^(-1/2) ||(N+1) T_sqrt(r) g||
on every test function.  Double precision; diagnostics, not proofs.
"""
import json, math, sys
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(20261110)
r, q = 2 / 7, 4 / 7
a = 1 / math.sqrt(1 - r * r)
NMAX = 700
x = np.linspace(-16, 16, 32001)
dx = x[1] - x[0]

def hermite_functions(kmax, y):
    psi = np.zeros((kmax + 1, y.size))
    psi[0] = np.exp(-y * y / 4) / (2 * math.pi) ** 0.25
    if kmax >= 1:
        psi[1] = y * psi[0]
    for k in range(1, kmax):
        psi[k + 1] = (y * psi[k] - math.sqrt(k) * psi[k - 1]) / math.sqrt(k + 1)
    return psi

sqrtphi = np.exp(-x * x / 4) / (2 * math.pi) ** 0.25
PSI = hermite_functions(NMAX - 1, x) * sqrtphi * dx        # G = PSI @ g
n = np.arange(NMAX)
wq, wr = q ** n, r ** n
KB = 120                                                    # oscillator levels used as test basis
HB = math.sqrt(2) * (2 * math.pi) ** 0.25 * hermite_functions(KB, math.sqrt(2 * a) * x)  # h_k(x)
L2sq = 2 * math.sqrt(math.pi / a)

ing = json.load(open(sys.argv[1]))
perk = ing["per_k"]

def norms(g):
    G = PSI @ g
    A = float(np.sum(n ** 2 * wq * G * G)); B = float(np.sum(wq * G * G))
    C = float(np.sum((n + 1) ** 2 * wr * G * G)); D0 = float(np.sum(wr * G * G))
    return G, A, B, C, D0

def ratios(g):
    G, A, B, C, D0 = norms(g)
    sup = float(np.max(np.abs(g)))
    return B / (sup * math.sqrt(D0)), math.sqrt(A * B) / (sup * math.sqrt(C))

def level_bounds(g):
    """max over even k <= KB of |c_k| / (its two proven upper bounds); must be <= 1."""
    G, A, B, C, D0 = norms(g)
    worst0 = worst1 = 0.0
    for k in range(0, KB + 1, 2):
        ck = abs(float(np.sum(g * HB[k]) * dx)) / L2sq
        Sr = float(perk[str(k)]["S_r"])
        if Sr < 1e-200 or ck < 1e-13:
            continue
        b0 = math.sqrt(D0 / Sr)
        b1 = perk[str(k)]["cL"] / (k + 1) * math.sqrt(C / Sr)
        worst0 = max(worst0, ck / b0); worst1 = max(worst1, ck / b1)
    return worst0, worst1

res = {"families": {}, "ok": True}
def record(name, gs):
    best0 = best1 = 0.0; lb0 = lb1 = 0.0
    for g in gs:
        r0, r1 = ratios(g); best0 = max(best0, r0); best1 = max(best1, r1)
        w0, w1 = level_bounds(g); lb0 = max(lb0, w0); lb1 = max(lb1, w1)
    res["families"][name] = {"count": len(gs), "max_rho0": best0, "max_rho1": best1,
                             "max_level_bound_plain": lb0, "max_level_bound_weighted": lb1}
    if best0 > ing["J_end"] or best1 > ing["C_E1"] or lb0 > 1 + 1e-6 or lb1 > 1 + 1e-6:
        res["ok"] = False
    print(f"{name:34s} n={len(gs):4d} rho0<={best0:.4f} rho1<={best1:.4f} level<={lb0:.4f},{lb1:.4f}")

# ---- Gram matrix at the level r: <T_sqrt(r) h_j, T_sqrt(r) h_k> = delta_jk S_k(r)   (fact R2)
GH = (PSI @ HB[:61].T)                                   # Hermite coefficients of h_0..h_60
Gram = (GH * wr[:, None]).T @ GH
dg = np.sqrt(np.diag(Gram))
# off-diagonal entries relative to the larger diagonal entry: double-precision roundoff in the
# quadrature is about 1e-15 in absolute terms, so compare |G_jk| with max(G_jj, G_kk)
off = np.abs(Gram - np.diag(np.diag(Gram))) / np.maximum.outer(np.diag(Gram), np.diag(Gram))
near = max(abs(Gram[j, k]) / (dg[j] * dg[k]) for j in range(31) for k in range(31) if 0 < abs(j - k) <= 4 and max(j, k) <= 16)
diag_err = max(abs(Gram[k, k] / float(perk[str(k)]["S_r"]) - 1) for k in range(0, 31, 2))
res["gram_offdiag_over_larger_diagonal_k<=30"] = float(off[:31, :31].max())
res["gram_near_diagonal_normalised_k<=16"] = float(near)
res["gram_diag_rel_err_k<=30"] = float(diag_err)
print("Gram at level r: max |G_jk|/max(G_jj,G_kk) (k<=30) %.2e; near-diagonal normalised (k<=16) %.2e; diagonal rel err %.2e"
      % (off[:31, :31].max(), near, diag_err))
if off[:31, :31].max() > 1e-6 or near > 1e-9 or diag_err > 1e-6:
    res["ok"] = False

def hka(k, aa):
    return hermite_functions(k, math.sqrt(2 * aa) * x)[k]

record("h_k, a=a_r, k even <=120", [HB[k] for k in range(0, KB + 1, 2)])
record("sign(h_k), a=a_r", [np.sign(HB[k]) for k in range(0, KB + 1, 2)])
record("clip(c h_k), a=a_r", [np.clip(c * HB[k], -1, 1) for k in range(0, KB + 1, 4) for c in (2, 5, 20)])
record("h_k other a", [hka(k, aa) for aa in (0.3, 0.5, 0.8, 1.0, 1.2, 1.6, 2.5) for k in range(0, 81, 4)])
record("cos(w x)", [np.cos(w * x) for w in np.linspace(0.2, 14, 70)])
record("chirps cos(b x^2 + w x)", [np.cos(b * x * x + w * x) for b in (0.1, 0.3, 1, 3) for w in (0, 2, 5, 9)])
record("spikes", [(np.abs(x - x0) < w).astype(float) for x0 in (0, 1, 2.5) for w in (1e-2, 0.1, 0.5, 2)])
record("gaussian packets", [np.exp(-al * (x - x0) ** 2) * np.cos(w * x)
                            for al in (0.05, 0.3, 1, 4) for x0 in (0, 1.5, 4) for w in (0, 3, 6, 10)])
record("windowed h_k (bulk, amplitude 1)",
       [np.clip(k ** 0.25 * HB[k] * np.exp(-(x / (c * math.sqrt(k))) ** 4), -1, 1)
        for k in range(8, KB + 1, 8) for c in (0.3, 0.6, 1.0)])
rand = []
for _ in range(200):
    K = int(rng.integers(4, KB)); m = int(rng.integers(1, 8))
    idx = rng.integers(0, K + 1, size=m); co = rng.normal(size=m) * rng.uniform(0.2, 30)
    g = np.tanh(sum(c * HB[i] for c, i in zip(co, idx)))
    rand.append(g)
record("tanh(random sums of h_k)", rand)

# ---- local maximisation of rho1 and rho0 over g = tanh(sum_k theta_k h_k), k even <= KO
def optimise(which, KO, starts):
    ks = np.arange(0, KO + 1, 2); H = HB[ks]
    def f(theta):
        g = np.tanh(theta @ H)
        G, A, B, C, D0 = norms(g)
        if which == 1:
            val = 0.5 * math.log(A) + 0.5 * math.log(B) - 0.5 * math.log(C)
            dG = n ** 2 * wq * G / A + wq * G / B - (n + 1) ** 2 * wr * G / C
        else:
            val = math.log(B) - 0.5 * math.log(D0)
            dG = 2 * wq * G / B - wr * G / D0
        grad = H @ ((1 - g * g) * (PSI.T @ dG))
        return -val, -grad
    best = (0.0, None)
    for s in range(starts):
        th0 = rng.normal(size=ks.size) * rng.uniform(0.3, 5)
        o = minimize(f, th0, jac=True, method="L-BFGS-B", options={"maxiter": 400})
        g = np.tanh(o.x @ H)
        val = ratios(g)[which]
        if val > best[0]:
            best = (val, g)
    return best
for which in (0, 1):
    for KO in (12, 40, 80):
        val, g = optimise(which, KO, 12)
        record(f"optimised rho{which}, levels<={KO}", [g])

json.dump(res, open(sys.argv[2], "w"), indent=1)
print("constants: J_end = %.4f, C_E1 = %.4f" % (ing["J_end"], ing["C_E1"]))
print("ALL OK" if res["ok"] else "FAILED")
sys.exit(0 if res["ok"] else 1)
