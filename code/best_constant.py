#!/usr/bin/env python3
"""Lower bounds for the best constant of the linear inequality at the critical level (Remark 5.3).

For a bounded g whose Hermite coefficients vanish below degree n, and s = b(r),

    ||T_sqrt(s) g||^2 / (||g||_inf ||T_sqrt(r) g||)  ->  (sqrt 2)^n |g_n| / ||g||_inf     as r -> 0.

This script builds exact +-1-valued functions g_n = sign(P_n), with P_n a polynomial of degree n whose
positive roots (the switching points of g_n) are solved at 50 digits so that every Hermite coefficient
of g_n below degree n vanishes.  For n = 1 this is g = sign(x) (Proposition F); for n = 2 it is
sign(x^2 - c^2) with c the median of |X| (the even example of Remark 5.3).  For each n in {1,2,4,6,8}
it computes the Hermite coefficients of g_n up to degree n + 120 by exact formulas
(integral of He_k phi over an interval is a difference of He_{k-1} phi), and evaluates the ratio at
r = 1e-2, 1e-5, 1e-8 together with the limit (sqrt 2)^n |g_n|.

Checks (prints ALL OK and exits 0 only if all hold):
  - every coefficient of g_n below degree n vanishes to 1e-30;
  - at n = 1 and r = 1e-2 the series value agrees with Sheppard's formula (2/pi) arcsin(b(r)) / sqrt((2/pi) arcsin r)
    to 1e-20, and the limit equals 2/sqrt(pi) to 1e-30;
  - at r = 1e-8 the ratio agrees with the limit to relative 1e-3, for every n;
  - the limits increase with n and exceed one.
mpmath at 50 digits; the series truncation at degree n + 120 is far beyond what r <= 1e-2 needs.

Usage: python3 code/best_constant.py OUT.json
"""
import json, sys
import numpy as np
from numpy.polynomial import hermite_e as HE
import mpmath as mp

mp.mp.dps = 50

def He(k, x):
    a, b = mp.mpf(1), x
    if k == 0: return a
    for j in range(1, k): a, b = b, x * b - j * a
    return b

phi = lambda x: mp.exp(-x ** 2 / 2) / mp.sqrt(2 * mp.pi)

def coef(k, pts):
    """k-th normalised Hermite coefficient <g, He_k/sqrt(k!)> of the +-1 function with switching points pts on (0,inf),
    value +1 on (max pts, inf), alternating inward, and extended evenly (k even) or oddly (k odd)."""
    m = len(pts); edges = [mp.mpf(0)] + list(pts) + [mp.inf]
    tot = mp.mpf(0)
    for i in range(len(edges) - 1):
        sgn = (-1) ** (m - i)
        a, b = edges[i], edges[i + 1]
        if k == 0:
            I = mp.ncdf(b) - mp.ncdf(a)
        else:
            Fa = He(k - 1, a) * phi(a); Fb = 0 if b == mp.inf else He(k - 1, b) * phi(b)
            I = Fa - Fb
        tot += sgn * I
    return 2 * tot / mp.sqrt(mp.factorial(k))

def switching_points(n):
    m = n // 2
    if m == 0: return []
    zs = sorted(mp.mpf(float(t)) for t in np.real(HE.hermeroots([0] * n + [1])) if t > 1e-9)
    par = n % 2
    eqs = lambda *p: [coef(2 * j + par, p) for j in range(m)]
    sol = mp.findroot(eqs, zs)
    return [sol] if m == 1 and not isinstance(sol, mp.matrix) else [sol[i] for i in range(m)]

def main():
    if len(sys.argv) != 2:
        print("usage: best_constant.py OUT.json"); return 2
    ok = True
    rs = [mp.mpf("1e-2"), mp.mpf("1e-5"), mp.mpf("1e-8")]
    out = {"dps": mp.mp.dps, "rows": [], "ok": True}
    limits = []
    for n in (1, 2, 4, 6, 8):
        pts = switching_points(n)
        low = max([abs(coef(j, pts)) for j in range(n % 2, n, 2)], default=mp.mpf(0))
        if not low < mp.mpf("1e-30"): ok = False
        gh = {k: coef(k, pts) for k in range(n, n + 121, 2)}
        lim = mp.sqrt(2) ** n * abs(gh[n])
        ratios = {}
        for r in rs:
            s = mp.sqrt(2 * r * (1 + r)) - r
            Ss = sum(v ** 2 * s ** k for k, v in gh.items()); Sr = sum(v ** 2 * r ** k for k, v in gh.items())
            ratios[mp.nstr(r, 3)] = Ss / mp.sqrt(Sr)
        if not abs(ratios["1.0e-8"] / lim - 1) < mp.mpf("1e-3"): ok = False
        row = {"n": n, "switching points": [mp.nstr(p, 20) for p in pts], "largest lower coefficient": float(low),
               "g_n": mp.nstr(gh[n], 20), "limit (sqrt2)^n |g_n|": float(lim),
               "ratio at r": {k: float(v) for k, v in ratios.items()}}
        if n == 1:
            r = rs[0]; s = mp.sqrt(2 * r * (1 + r)) - r
            shep = (2 / mp.pi) * mp.asin(s) / mp.sqrt((2 / mp.pi) * mp.asin(r))
            row["Sheppard value at r=1e-2"] = float(shep)
            row["series minus Sheppard at r=1e-2"] = float(abs(ratios["0.01"] - shep))
            if not abs(ratios["0.01"] - shep) < mp.mpf("1e-20"): ok = False
            if not abs(lim - 2 / mp.sqrt(mp.pi)) < mp.mpf("1e-30"): ok = False
        out["rows"].append(row); limits.append(lim)
        print(f"n={n}  limit={mp.nstr(lim, 6)}  ratio r=1e-2: {mp.nstr(ratios['0.01'], 6)}  "
              f"r=1e-5: {mp.nstr(ratios['1.0e-5'], 6)}  r=1e-8: {mp.nstr(ratios['1.0e-8'], 6)}  low={mp.nstr(low, 3)}")
    if not all(limits[i] < limits[i + 1] for i in range(len(limits) - 1)): ok = False
    if not limits[0] > 1: ok = False
    out["ok"] = ok
    with open(sys.argv[1], "w") as fh:
        json.dump(out, fh, indent=1)
    print("ALL OK" if ok else "FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
