#!/usr/bin/env python3
"""Three copies: numerical confirmation that concavity on the simplex fails near the Gaussian
at every tested signal fraction q, exactly as predicted by the boundary criterion

    F''(0)  ->  -4 h_4 K_4(h)     as eps -> 0,
    K_4(h) = < e_4 (*) (1+h), log(1+h) > = E_phi[(1+h) d^4/dx^4 log(1+h)] / sqrt(24),

along the segment t(tau) = (eps+tau, eps-tau, 1-2 eps), |tau| <= eps, where h = T_sqrt(q) g and
(*) is the Wick product  (u (*) v)_n = sum_{j+k=n} sqrt(C(n,j)) u_j v_k.

For each q the script builds an admissible g = d (g0 + theta g1) with sup-norm at most d, where
g0 has no Hermite component of degree below 8 and g1 has none of degree 0, 2, 6 but a nonzero
component of degree 4 (theta is chosen so that h_4 is about minus one half of K_4 at theta = 0), and

 (1) evaluates h_4, K_4 (two ways) and the second-order quantity Y_2 (two ways);
 (2) evaluates F'' at tau = 0 and tau = +-eps/2 on the segment from truncated Hermite series for
     eps = 1e-2, ..., 1e-8 and compares with -4 h_4 K_4.  The limit is approached from below, with
     a correction of relative size (eps/eps_0)^2 coming from the degree-8 part of g, so F'' turns
     positive only below an eps_0 that depends on q and on g;
 (3) confirms non-concavity independently of the F'' formula: the first derivative
     F'(eps/2) = -<R_tau, log(1+R)> is positive (F is symmetric in tau, so a concave F would be
     maximal at tau = 0); and, where the effect is large enough to resolve in double precision,
     F(0) < (F(tau) + F(-tau))/2 from directly computed values of F;
 (4) evaluates the two-copy curvature F''/D^2 for the same g and q (it must be negative).
Double precision.  The coefficient identities are exact; their evaluation is not: series are truncated
(NG input degrees, NOUT output degrees), integrals are Riemann sums on finite grids, and sup norms are
maxima over the grid.  The amplitude d = 1/4 exceeds the bound 1/12 assumed in the proposition and is
exploratory.  g0 has leading Hermite degree 8 (the existence proof in the paper approximates degree 6;
any g0 without components of degree 0, 2, 4 and with Y(T_sqrt(q) g0) != 0 serves).  Diagnostics.
"""
import json, math, sys
import numpy as np
from scipy.special import gammaln

# Numerical parameters.  Default: input Hermite degrees < 500, output degrees < 260, input grid on
# [-16,16] with mesh 1e-3, output grid on [-9,9] with mesh 2.5e-3.  With --fine (used only for the
# convergence table): degrees < 600 and < 320, grids on [-18,18] and [-10,10], both meshes halved.
FINE = "--fine" in sys.argv
NG, NOUT = (600, 320) if FINE else (500, 260)
XG_MAX, XG_N, X_MAX, X_N = (18, 72001, 10, 16001) if FINE else (16, 32001, 9, 7201)
xg = np.linspace(-XG_MAX, XG_MAX, XG_N); dxg = xg[1] - xg[0]
x = np.linspace(-X_MAX, X_MAX, X_N); dx = x[1] - x[0]
phi = np.exp(-x * x / 2) / math.sqrt(2 * math.pi)
phig = np.exp(-xg * xg / 2) / math.sqrt(2 * math.pi)
def hf(kmax, y):
    psi = np.zeros((kmax + 1, y.size)); psi[0] = np.exp(-y * y / 4) / (2 * math.pi) ** 0.25
    psi[1] = y * psi[0]
    for k in range(1, kmax): psi[k + 1] = (y * psi[k] - math.sqrt(k) * psi[k - 1]) / math.sqrt(k + 1)
    return psi
PSIG = hf(NG - 1, xg) * np.sqrt(phig) * dxg          # Hermite coefficients of a function on xg
PSIX = hf(NOUT - 1, x) / np.sqrt(phi)                # e_n(x) on x
lg = gammaln(np.arange(NOUT + 1) + 1.0)
nn = np.arange(NOUT)
SQB = np.zeros((NOUT, NOUT))                         # SQB[n, j] = sqrt(C(n, j)), j <= n
for n in range(NOUT):
    j = np.arange(n + 1); SQB[n, :n + 1] = np.exp(0.5 * (lg[n] - lg[j] - lg[n - j]))

def wick(u, v):
    out = np.zeros(NOUT)
    for n in range(NOUT):
        out[n] = np.dot(SQB[n, :n + 1] * u[:n + 1], v[n::-1])
    return out
def on_grid(c): return c @ PSIX
def corrected(base, degrees):
    """Subtract bounded even correctors exp(-x^2/2) x^(2j) so the listed Hermite degrees vanish."""
    cor = [np.exp(-xg * xg / 2) * xg ** (2 * j) for j in range(len(degrees))]
    co = lambda f: np.array([(PSIG @ f)[d] for d in degrees])
    A = np.column_stack([co(c) for c in cor])
    g = base - sum(b * c for b, c in zip(np.linalg.solve(A, co(base)), cor))
    return 0.5 * g / np.max(np.abs(g))
def coeffs(g, q):
    c = (PSIG @ g)[:NOUT] * q ** (nn / 2); c[1::2] = 0.0; c[:4] = np.where(np.abs(c[:4]) < 1e-10, 0.0, c[:4])
    return c

def Y2_series(h):
    a = sum(SQB[n, 4] * h[n - 4] * h[n] for n in range(8, NOUT))
    hx = on_grid(h); e4 = PSIX[4]
    return a - 0.5 * float(np.sum(e4 * hx * hx * phi) * dx)
def Y2_closed(h):
    n = nn.astype(float)
    d = 3 * np.sum(h * h * n * (n - 1)) + 4 * np.sum(h[:-2] * h[2:] * n[:-2] * np.sqrt((n[:-2] + 1) * (n[:-2] + 2)))
    return -d / math.sqrt(24)
def K4(h):
    one = np.zeros(NOUT); one[0] = 1.0; e4 = np.zeros(NOUT); e4[4] = 1.0
    lhs = on_grid(wick(e4, one + h)); hx = on_grid(h)
    k_series = float(np.sum(lhs * np.log1p(hx) * phi) * dx)
    L = np.log1p(hx)                                  # second route: E[(1+h) d^4 log(1+h)] / sqrt(24)
    d4 = np.gradient(np.gradient(np.gradient(np.gradient(L, dx), dx), dx), dx)
    k_direct = float(np.sum((1 + hx) * d4 * phi) * dx) / math.sqrt(24)
    return k_series, k_direct

def segment(h, eps, tau):
    """Coefficients of P = A (*) B (*) W and its first two tau-derivatives; P = 1 + R."""
    m = nn / 2.0
    def fac(t):
        a = np.zeros(NOUT); a1 = np.zeros(NOUT); a2 = np.zeros(NOUT); a[0] = 1.0
        ev = np.arange(4, NOUT, 2); mm = ev // 2
        a[ev] = h[ev] * t ** mm; a1[ev] = h[ev] * mm * t ** (mm - 1); a2[ev] = h[ev] * mm * (mm - 1) * t ** (mm - 2)
        return a, a1, a2
    A, A1, A2 = fac(eps + tau); B, B1, B2 = fac(eps - tau); W = fac(1 - 2 * eps)[0]
    P = wick(wick(A, B), W)
    P1 = wick(wick(A1, B) - wick(A, B1), W)
    P2 = wick(wick(A2, B) - 2 * wick(A1, B1) + wick(A, B2), W)
    return P, P1, P2
def F_and_curv(h, eps, tau):
    P, P1, P2 = segment(h, eps, tau)
    Px, P1x, P2x = on_grid(P), on_grid(P1), on_grid(P2)
    F = -float(np.sum(Px * np.log(Px) * phi) * dx)
    F2 = -float(np.sum(P2x * np.log(Px) * phi) * dx) - float(np.sum(P1x ** 2 / Px * phi) * dx)
    F1 = -float(np.sum(P1x * np.log(Px) * phi) * dx)
    return F, F2, F1

def two_copy_ratio(h, t):
    m = np.arange(4, NOUT, 2) // 2; ev = np.arange(4, NOUT, 2)
    def fac(tt):
        a = np.zeros(NOUT); a1 = np.zeros(NOUT); a2 = np.zeros(NOUT); a[0] = 1.0
        a[ev] = h[ev] * tt ** m; a1[ev] = h[ev] * m * tt ** (m - 1); a2[ev] = h[ev] * m * (m - 1) * tt ** (m - 2)
        return a, a1, a2
    A, A1, A2 = fac(t); B, B1, B2 = fac(1 - t)
    P = wick(A, B); P1 = wick(A1, B) - wick(A, B1); P2 = wick(A2, B) - 2 * wick(A1, B1) + wick(A, B2)
    Px, P1x, P2x = on_grid(P), on_grid(P1), on_grid(P2)
    F2 = -float(np.sum(P2x * np.log(Px) * phi) * dx) - float(np.sum(P1x ** 2 / Px * phi) * dx)
    Lc = np.zeros(NOUT); Lc[ev] = h[ev] * (t ** m + (1 - t) ** m)
    return F2 / float(np.sum((nn * Lc) ** 2))

out = {"ok": True, "cases": []}
out["numerics"] = {"input degrees": NG, "output degrees": NOUT, "input grid": [-XG_MAX, XG_MAX, XG_N],
                   "output grid": [-X_MAX, X_MAX, X_N], "positions sampled on each segment": ["0", "+eps/2", "-eps/2"],
                   "sup norms": "maximum over the input grid"}
g0t = corrected(np.cos(2.5 * xg), [0, 2, 4, 6])
g1t = corrected(np.cos(1.5 * xg), [0, 2, 6])
EPS = (1e-2, 1e-3, 1e-4, 1e-6, 1e-8)
for q in (0.2, 0.5, 4 / 7, 0.9):
    for d in (0.25, 1 / 12):
        h0 = coeffs(d * g0t, q); h1 = coeffs(d * g1t, q)
        Y0 = Y2_series(h0); Y0c = Y2_closed(h0)
        K0 = K4(h0)[0]                                       # K_4 at theta = 0 (h_4 = 0 there)
        theta = -0.5 * K0 / h1[4]                            # then h_4 = -K0/2
        g = d * (g0t + theta * g1t)
        sup = float(np.max(np.abs(g)))
        if sup > d: g *= d / sup; sup = d
        h = coeffs(g, q)
        gc = PSIG @ g
        low = float(max(abs(gc[0]), abs(gc[2])))
        k_s, k_d = K4(h)
        pred = -4 * h[4] * k_s
        rows = []
        for eps in EPS:
            F0, F2, _ = F_and_curv(h, eps, 0.0)
            Fp, F2p, F1p = F_and_curv(h, eps, eps / 2); Fm, F2m, F1m = F_and_curv(h, eps, -eps / 2)
            gap = 0.5 * (Fp + Fm) - F0
            expected_gap = 0.5 * F2 * (eps / 2) ** 2
            rows.append({"eps": eps, "F''(0)": F2, "F''(0)/prediction": F2 / pred, "F''(+eps/2)": F2p, "F''(-eps/2)": F2m,
                         "F'(+eps/2)": F1p, "F'(-eps/2)": F1m, "F'(+eps/2)/(F''(0) eps/2)": F1p / (F2 * eps / 2),
                         "chord gap": gap, "chord gap/expected": gap / expected_gap,
                         "chord gap resolvable": bool(abs(expected_gap) > 3e-15)})
        pos = [r["eps"] for r in rows if r["F''(0)"] > 0 and r["F''(+eps/2)"] > 0 and r["F''(-eps/2)"] > 0
               and r["F'(+eps/2)"] > 0 and r["F'(-eps/2)"] < 0]
        two = max(two_copy_ratio(h, t) for t in (0.0, 0.05, 0.2, 0.5))
        case = {"q": q, "d": d, "sup norm of g": sup, "theta": theta, "low-degree residual": low,
                "Y2(h0) series": Y0, "Y2(h0) closed form": Y0c, "h_4": h[4], "K_4 (series)": k_s,
                "K_4 (direct)": k_d, "h_4*K_4": h[4] * k_s, "prediction -4 h_4 K_4": pred,
                "segment": rows, "two-copy max F''/D^2": two}
        last = rows[-1]
        ok = (sup <= d * (1 + 1e-9) and low < 1e-9 and abs(Y0 / Y0c - 1) < 1e-6 and abs(k_s / k_d - 1) < 1e-3
              and h[4] * k_s < 0 and set(pos) >= {1e-6, 1e-8}
              and abs(last["F''(0)/prediction"] - 1) < 1e-3 and two < -1 / 16
              and all(r["chord gap"] > 0 and abs(r["chord gap/expected"] - 1) < 0.2 for r in rows
                      if r["chord gap resolvable"] and r["F''(0)"] > 0))
        case["largest tested eps with convexity confirmed"] = max(pos) if pos else None
        case["ok"] = bool(ok); out["ok"] = out["ok"] and bool(ok); out["cases"].append(case)
        res = [r["eps"] for r in rows if r["chord gap resolvable"] and r["F''(0)"] > 0]
        print(f"q={q:.3f} d={d:.4f} sup={sup:.3f} h4={h[4]:+.2e} K4={k_s:+.2e} limit={pred:+.2e} F''(0)/limit: "
              + " ".join(f"{r['F'+chr(39)*2+'(0)/prediction']:+.3f}" for r in rows)
              + f" | convex for eps<={case['largest tested eps with convexity confirmed']} | chord ok at {res} | two-copy {two:+.3f} | {'ok' if ok else 'FAIL'}")
json.dump(out, open([a for a in sys.argv[1:] if not a.startswith("--")][0], "w"), indent=1)
print("ALL OK" if out["ok"] else "FAILED"); sys.exit(0 if out["ok"] else 1)
