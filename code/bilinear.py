#!/usr/bin/env python3
"""Checks of the bilinear term M_t at unit smoothing (rho = 1).

 (1) the coefficient formula  M_n(t) = sum_{j+k=n} h_j h_k sqrt(binom(n,j) t^j (1-t)^k)
     against the conditional-expectation definition, by Gauss-Hermite quadrature;
 (2) the second-derivative bracket
       d^2/dt^2 sqrt(b_n(j;t)) = [ (j-nt)^2/(4 t^2 (1-t)^2) - (j/t^2 + k/(1-t)^2)/2 ] sqrt(b_n(j;t))
     against central finite differences;
 (3) the three norm bounds, for all 0 <= t <= 1 (h with even Hermite degrees >= 4):
       ||M|| <= ||h||^2,  ||N M|| <= 2 ||N h|| ||h||,  ||N^{-1} M_tt|| <= 2.97 ||N h|| ||h||.
Double precision diagnostics.  Exit code is non-zero on failure.
"""
import json, math, sys
import numpy as np
from numpy.polynomial.hermite_e import hermegauss, hermeval
from scipy.special import gammaln

rng = np.random.default_rng(7)
out = {"ok": True}

def sqrt_b(n, j, t):
    return math.exp(0.5 * (gammaln(n + 1) - gammaln(j + 1) - gammaln(n - j + 1))) * t ** (j / 2) * (1 - t) ** ((n - j) / 2)

def M_coeffs(h, t):
    J = len(h) - 1
    M = np.zeros(2 * J + 1)
    for j in range(J + 1):
        if h[j] == 0: continue
        for k in range(J + 1):
            if h[k] == 0: continue
            M[j + k] += h[j] * h[k] * sqrt_b(j + k, j, t)
    return M

def Mtt_coeffs(h, t):
    J = len(h) - 1
    M = np.zeros(2 * J + 1)
    for j in range(J + 1):
        if h[j] == 0: continue
        for k in range(J + 1):
            if h[k] == 0: continue
            n = j + k
            al, be = j // 2, k // 2          # d^2/dt^2 of t^al (1-t)^be, valid at t = 0 and t = 1 too
            d2 = (al * (al - 1) * t ** (al - 2) * (1 - t) ** be - 2 * al * be * t ** (al - 1) * (1 - t) ** (be - 1)
                  + be * (be - 1) * t ** al * (1 - t) ** (be - 2))
            M[n] += h[j] * h[k] * math.exp(0.5 * (gammaln(n + 1) - gammaln(j + 1) - gammaln(k + 1))) * d2
    return M

def e(n, x):
    c = np.zeros(n + 1); c[n] = 1
    return hermeval(x, c) / math.sqrt(math.factorial(n))

# (1) definition check
xg, wg = hermegauss(80); wg = wg / math.sqrt(2 * math.pi)
h = np.zeros(9); h[4], h[6], h[8] = 0.7, -0.4, 0.25
worst = 0.0
for t in (0.25, 0.4, 0.5, 0.7):
    c, s = math.sqrt(t), math.sqrt(1 - t)
    Y = xg[:, None]; Yp = xg[None, :]                 # Y = c G1 + s G2, Yp = -s G1 + c G2
    G1 = c * Y - s * Yp; G2 = s * Y + c * Yp
    hv = lambda z: sum(h[n] * e(n, z) for n in (4, 6, 8))
    Mx = np.sum(wg[None, :] * hv(G1) * hv(G2), axis=1)   # conditional expectation as a function of Y
    Mc = M_coeffs(h, t)
    for n in range(0, 17, 2):
        proj = float(np.sum(wg * Mx * e(n, xg)))
        worst = max(worst, abs(proj - Mc[n]))
out["coefficient_formula_max_abs_err"] = worst
if worst > 1e-9: out["ok"] = False

# (2) derivative bracket vs finite differences, (3) norm bounds
r1 = r2 = r3 = fd = 0.0
for trial in range(400):
    J = int(rng.choice([8, 12, 20, 40, 60]))
    h = np.zeros(J + 1)
    idx = np.arange(4, J + 1, 2)
    style = trial % 4
    if style == 0: h[idx] = rng.normal(size=idx.size)
    elif style == 1: h[idx] = rng.normal(size=idx.size) * 0.6 ** idx
    elif style == 2: h[rng.choice(idx)] = 1.0
    else:
        m = rng.choice(idx); h[idx] = np.exp(-(idx - m) ** 2 / 18.0) * rng.choice([-1, 1], size=idx.size)
    t = [0.0, 1.0, 1e-3, 0.02, 0.5][trial % 5] if trial % 3 == 0 else rng.uniform(0.0, 1.0)
    M = M_coeffs(h, t); Mtt = Mtt_coeffs(h, t)
    if trial < 90 and 0.05 < t < 0.95:
        d = 1e-4
        num = (M_coeffs(h, t + d) - 2 * M + M_coeffs(h, t - d)) / d ** 2
        fd = max(fd, float(np.max(np.abs(num - Mtt)) / (1e-12 + np.max(np.abs(Mtt)))))
    n = np.arange(M.size)
    hh = float(np.sum(h * h)); Nh = math.sqrt(float(np.sum((np.arange(J + 1) * h) ** 2)))
    r1 = max(r1, math.sqrt(float(np.sum(M * M))) / hh)
    r2 = max(r2, math.sqrt(float(np.sum((n * M) ** 2))) / (2 * Nh * math.sqrt(hh)))
    r3 = max(r3, math.sqrt(float(np.sum((Mtt[1:] / n[1:]) ** 2))) / (Nh * math.sqrt(hh)))
out.update({"finite_difference_rel_err": fd, "max ||M||/||h||^2": r1,
            "max ||NM||/(2||Nh|| ||h||)": r2, "max ||N^-1 M_tt||/(||Nh|| ||h||)": r3})
if fd > 1e-4 or r1 > 1 + 1e-12 or r2 > 1 + 1e-12 or r3 > 2.97:
    out["ok"] = False
json.dump(out, open(sys.argv[1], "w"), indent=1)
print(json.dumps(out, indent=1)); print("ALL OK" if out["ok"] else "FAILED")
sys.exit(0 if out["ok"] else 1)
