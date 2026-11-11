#!/usr/bin/env python3
"""End-to-end check of the concavity statement at q = 4/7.

For admissible perturbations g (even, bounded, with zero Gaussian mean and second moment) scaled to
sup-norm delta, the script builds the exact Hermite coefficients of the output density ratio
1 + R_t = 1 + L_t + M_t and of its first two t-derivatives, evaluates

    F''(t) = -<R_tt, log(1+R)> - int R_t^2/(1+R) phi

on a grid, and compares it with -D^2/16, D = ||N L_t||.  It also compares F'' with a second
difference of F(t) = h(phi) - int (1+R) log(1+R) phi computed independently from the same series.
Double precision; a diagnostic of the statement, not a proof.
"""
import json, math, sys
import numpy as np
from scipy.special import gammaln

q = 4 / 7
NG, NOUT = 400, 240                     # input degrees kept, output degrees kept (q^(n/2) decay)
x = np.linspace(-9, 9, 7201); dx = x[1] - x[0]
phi = np.exp(-x * x / 2) / math.sqrt(2 * math.pi)
def hf(kmax, y):
    psi = np.zeros((kmax + 1, y.size)); psi[0] = np.exp(-y * y / 4) / (2 * math.pi) ** 0.25
    psi[1] = y * psi[0]
    for k in range(1, kmax): psi[k + 1] = (y * psi[k] - math.sqrt(k) * psi[k - 1]) / math.sqrt(k + 1)
    return psi
xg = np.linspace(-16, 16, 32001); dxg = xg[1] - xg[0]
PSIG = hf(NG - 1, xg) * (np.exp(-xg * xg / 4) / (2 * math.pi) ** 0.25) * dxg     # coefficients of g
PSIX = hf(NOUT - 1, x) / np.sqrt(phi)                                             # e_n(x) on the grid
lgam = gammaln(np.arange(NOUT + 2) + 1.0)

def admissible(base):
    """Subtract a combination of two bounded even functions so that degrees 0 and 2 vanish."""
    c1, c2 = np.exp(-xg * xg / 2), (1 - xg * xg) * np.exp(-xg * xg / 2)
    phig = np.exp(-xg * xg / 2) / math.sqrt(2 * math.pi)
    mom = lambda f: np.array([np.sum(f * phig) * dxg, np.sum(xg * xg * f * phig) * dxg])
    A = np.column_stack([mom(c1), mom(c2)])
    co = np.linalg.solve(A, mom(base))
    g = base - co[0] * c1 - co[1] * c2
    return g / np.max(np.abs(g))

def series(gh, t):
    """Hermite coefficients of R, R_t, R_tt at weight t (gh = q^(n/2) ghat_n, even n >= 4)."""
    R = np.zeros(NOUT); R1 = np.zeros(NOUT); R2 = np.zeros(NOUT); L = np.zeros(NOUT)
    ev = [n for n in range(4, NOUT, 2) if gh[n] != 0.0]
    for n in ev:
        m = n // 2
        l0 = t ** m + (1 - t) ** m
        l1 = m * (t ** (m - 1) - (1 - t) ** (m - 1))
        l2 = m * (m - 1) * (t ** (m - 2) + (1 - t) ** (m - 2))
        L[n] = gh[n] * l0; R[n] += gh[n] * l0; R1[n] += gh[n] * l1; R2[n] += gh[n] * l2
    for j in ev:
        for k in ev:
            n = j + k
            if n >= NOUT: break
            w = gh[j] * gh[k] * math.exp(0.5 * (lgam[n] - lgam[j] - lgam[k]))
            al, be = j // 2, k // 2
            R[n] += w * t ** al * (1 - t) ** be
            R1[n] += w * (al * t ** (al - 1) * (1 - t) ** be - be * t ** al * (1 - t) ** (be - 1))
            R2[n] += w * (al * (al - 1) * t ** (al - 2) * (1 - t) ** be - 2 * al * be * t ** (al - 1) * (1 - t) ** (be - 1)
                          + be * (be - 1) * t ** al * (1 - t) ** (be - 2))
    return R, R1, R2, L

def Fpp(gh, t):
    R, R1, R2, L = series(gh, t)
    Rx, R1x, R2x = R @ PSIX, R1 @ PSIX, R2 @ PSIX
    val = -np.sum(R2x * np.log1p(Rx) * phi) * dx - np.sum(R1x ** 2 / (1 + Rx) * phi) * dx
    D2 = float(np.sum((np.arange(NOUT) * L) ** 2))
    return float(val), D2
def F(gh, t):
    R = series(gh, t)[0] @ PSIX
    return -float(np.sum((1 + R) * np.log1p(R) * phi) * dx)

H = hf(60, math.sqrt(2 * 7 / (3 * math.sqrt(5))) * xg)
bases = {"cos(2x)": np.cos(2 * xg), "cos(4x)": np.cos(4 * xg), "cos(7x)": np.cos(7 * xg),
         "sign(cos(3x))": np.sign(np.cos(3 * xg)), "h_8": H[8], "h_20": H[20], "h_40": H[40],
         "sign(h_12)": np.sign(H[12]), "sign(h_30)": np.sign(H[30]),
         "indicator(|x|<1)": (np.abs(xg) < 1).astype(float), "chirp cos(x^2)": np.cos(xg * xg),
         "tanh(5 h_16 + 3 h_24)": np.tanh(5 * H[16] + 3 * H[24])}
delta_star = json.load(open(sys.argv[1]))["delta_star_indicative"]
out = {"delta_star": delta_star, "rows": [], "ok": True}
worst = {}
qn = q ** (np.arange(NG) / 2)
for name, base in bases.items():
    g = admissible(base)
    gc = PSIG @ g
    resid = max(abs(gc[0]), abs(gc[2]))
    for delta in (delta_star, 0.01, 0.05, 1 / 12):
        gh = np.zeros(NOUT); gh[:NOUT] = (delta * gc * qn)[:NOUT]; gh[:4] = 0.0; gh[1::2] = 0.0
        for t in (0.0, 0.1, 0.3, 0.5):
            v, D2 = Fpp(gh, t)
            ratio = v / D2
            row = {"g": name, "delta": delta, "t": t, "F''/D^2": ratio}
            if t in (0.3, 0.5) and delta >= 0.05 and D2 > 1e-6:   # below this the second difference is rounding noise
                d = 2e-3
                fd = (F(gh, t + d) - 2 * F(gh, t) + F(gh, t - d)) / d ** 2
                row["second_difference/D^2"] = fd / D2
                if abs(fd - v) > 2e-3 * abs(v): out["ok"] = False; row["FAIL_fd"] = True
            out["rows"].append(row)
            worst[delta] = max(worst.get(delta, -1e9), ratio)
    if resid > 1e-9: out["ok"] = False
out["max F''/D^2 by delta"] = {str(k): v for k, v in worst.items()}
if worst[delta_star] > -1 / 16: out["ok"] = False
json.dump(out, open(sys.argv[2], "w"), indent=1)
print(json.dumps(out["max F''/D^2 by delta"], indent=1)); print("threshold -1/16 =", -1 / 16)
print("ALL OK" if out["ok"] else "FAILED"); sys.exit(0 if out["ok"] else 1)
