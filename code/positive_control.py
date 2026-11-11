#!/usr/bin/env python3
"""Positive control for the end-to-end computation of F'' (Theorems C and D).

entropy_curvature.py finds F''(t)/D^2 < 0 at q = 4/7 for every test function.  That is informative only if
the same kind of computation can see a failure where one is known to exist.  This script evaluates the
midpoint curvature F''(1/2) for the perturbations used in Theorem D,

    g = delta * h~_N,   h~_N = (h_N - c_0 u_0 - c_2 u_2) / S,   u_0 = e^{-x^2/2},  u_2 = (1 - x^2) e^{-x^2/2},

with a = a_{q/2}, delta = 1/12, N in {30, 40, 50, 80}, at q = 4/7 (Theorem C: concave) and at q = 0.7, 0.8
(Theorem D: concavity fails near the Gaussian).  S is the maximum of |h_N - c_0 u_0 - c_2 u_2| on the grid
[-16,16] (32001 points), so that ||g||_inf = delta on that grid.

Input coefficients.  The Hermite coefficients <h_N, e_n>, <u_0, e_n>, <u_2, e_n> are computed exactly from
Gaussian generating functions in mpmath (90 digits), and c_0, c_2 are solved exactly so that the degree-0 and
degree-2 coefficients of h~_N vanish.  There is therefore no quadrature floor in the input coefficients.  As an
accuracy check the whole computation is repeated at 60 digits, and the coefficients of h_N are compared with
double-precision quadrature for N = 30.

The curvature.  At t = 1/2 the first derivative R_t vanishes by symmetry (checked), so
F''(1/2) = -<R_tt, log(1+R)>.  The script evaluates
    pairing  = -sum_n R_tt,n R_n            over all output degrees (inputs < 400, outputs < 800: no truncation of M),
    nonlinear correction = <R_tt, R - log1p(R)>  by quadrature on [-9,9] (7201 points),
and reports F''(1/2) = pairing + correction together with D^2 = ||N L_1/2||^2 over all degrees.
For comparison it also reports the value of the original harness of entropy_curvature.py (outputs < 240,
quadrature of the synthesized series on [-9,9]) with the same exact inputs, which shows the bias of that window.

"resolved" (per row) means: the 60- and 90-digit input coefficients agree to 1e-40 relative to the largest
coefficient, and D^2 and F''(1/2) computed from them agree to 1e-12; the input tail beyond degree 400 carries
less than 1e-20 of sum_n q^n g_n^2; and the nonlinear correction (the only part evaluated on a window) is below
1e-3 of the pairing.

Checks (prints ALL OK and exits 0 only if all hold):
  - every row is resolved;  R_t(1/2) = 0 to 1e-12 relative (double rounding);  quadrature check at N = 30 to 1e-12;
  - at q = 4/7, F''(1/2)/D^2 < -0.4 for every N;
  - at q = 0.7 and at q = 0.8, F''(1/2) > 0 for at least one N.
Usage: python3 code/positive_control.py OUT.json
"""
import json, math, sys
import numpy as np
import mpmath as mp
from scipy.special import gammaln

NG, NBIG, NOUT = 400, 800, 240          # input degrees, output degrees (complete), harness output degrees
NTAIL = 600                             # input degrees used to measure the tail beyond NG
x = np.linspace(-9, 9, 7201); dx = x[1] - x[0]
phi = np.exp(-x * x / 2) / math.sqrt(2 * math.pi)
def hf(kmax, y):
    psi = np.zeros((kmax + 1, y.size)); psi[0] = np.exp(-y * y / 4) / (2 * math.pi) ** 0.25
    psi[1] = y * psi[0]
    for k in range(1, kmax): psi[k + 1] = (y * psi[k] - math.sqrt(k) * psi[k - 1]) / math.sqrt(k + 1)
    return psi
xg = np.linspace(-16, 16, 32001); dxg = xg[1] - xg[0]
PSIX = hf(NBIG - 1, x) / np.sqrt(phi)                  # e_n(x) on the output grid
lgam = gammaln(np.arange(NBIG + 2) + 1.0)

def exact_coeffs(N, a, dps):
    """Hermite coefficients (n < NTAIL) of h_N, u_0, u_2, at dps digits, from generating functions.
    E_phi[ sum_k h_k z^k/sqrt(k!) * sum_n e_n w^n/sqrt(n!) ] = sqrt(2/(1+a)) exp(al z^2 + be w^2 + ga z w)."""
    with mp.workdps(dps):
        a = mp.mpf(a)
        al = (a - 1) / (2 * (1 + a)); be = -a / (2 * (1 + a)); ga = mp.sqrt(2 * a) / (1 + a)
        pre = mp.sqrt(2 / (1 + a))
        fN = mp.factorial(N)
        H = [mp.mpf(0)] * NTAIL; U0 = [mp.mpf(0)] * NTAIL; U2 = [mp.mpf(0)] * NTAIL
        for n in range(N % 2, NTAIL, 2):
            s = mp.mpf(0)
            for m in range(N % 2, min(N, n) + 1, 2):
                i, j = (N - m) // 2, (n - m) // 2
                s += ga ** m / mp.factorial(m) * al ** i / mp.factorial(i) * be ** j / mp.factorial(j)
            H[n] = pre * mp.sqrt(fN * mp.factorial(n)) * s
        # u_0 = e^{-x^2/2}:   E_phi[u_0 e^{Xw - w^2/2}] = 2^{-1/2} e^{-w^2/4}
        # x^2 u_0:            2^{-1/2} (1/2 + w^2/4) e^{-w^2/4}
        for n in range(0, NTAIL, 2):
            m = n // 2
            c = (mp.mpf(-1) / 4) ** m / mp.factorial(m)                     # [w^n] e^{-w^2/4}
            cm1 = (mp.mpf(-1) / 4) ** (m - 1) / mp.factorial(m - 1) if m >= 1 else mp.mpf(0)
            sq = mp.sqrt(mp.factorial(n)) / mp.sqrt(2)
            U0[n] = sq * c
            X2 = sq * (c / 2 + cm1 / 4)
            U2[n] = U0[n] - X2
        # c_0, c_2 with degree-0 and degree-2 coefficients of h_N - c_0 u_0 - c_2 u_2 equal to zero
        if N % 2 == 0:
            A = mp.matrix([[U0[0], U2[0]], [U0[2], U2[2]]]); rhs = mp.matrix([H[0], H[2]])
            c = mp.lu_solve(A, rhs); c0, c2 = c[0], c[1]
        else:
            c0 = c2 = mp.mpf(0)
        G = [H[n] - c0 * U0[n] - c2 * U2[n] for n in range(NTAIL)]
        return G, H, float(c0), float(c2)

def series(gh, t, nout):
    """Coefficients of R, R_t, R_tt, L below degree nout (same formulas as entropy_curvature.py)."""
    R = np.zeros(nout); R1 = np.zeros(nout); R2 = np.zeros(nout); L = np.zeros(nout)
    ev = [n for n in range(4, min(nout, gh.size), 2) if gh[n] != 0.0]
    for n in ev:
        m = n // 2
        l0 = t ** m + (1 - t) ** m
        l1 = m * (t ** (m - 1) - (1 - t) ** (m - 1))
        l2 = m * (m - 1) * (t ** (m - 2) + (1 - t) ** (m - 2))
        L[n] = gh[n] * l0; R[n] += gh[n] * l0; R1[n] += gh[n] * l1; R2[n] += gh[n] * l2
    for j in ev:
        for k in ev:
            n = j + k
            if n >= nout: break
            w = gh[j] * gh[k] * math.exp(0.5 * (lgam[n] - lgam[j] - lgam[k]))
            al, be = j // 2, k // 2
            R[n] += w * t ** al * (1 - t) ** be
            R1[n] += w * (al * t ** (al - 1) * (1 - t) ** be - be * t ** al * (1 - t) ** (be - 1))
            R2[n] += w * (al * (al - 1) * t ** (al - 2) * (1 - t) ** be - 2 * al * be * t ** (al - 1) * (1 - t) ** (be - 1)
                          + be * (be - 1) * t ** al * (1 - t) ** (be - 2))
    return R, R1, R2, L

def evaluate(Gmp, scale, q, delta):
    """Exact-input evaluation: returns dict with pairing value, correction, harness value, D^2 and diagnostics."""
    G = np.array([float(v) for v in Gmp]) / scale
    qn = q ** (np.arange(NTAIL) / 2)
    full = delta * G * qn; full[:4] = 0.0; full[1::2] = 0.0
    tail = float(np.sum(full[NG:] ** 2) / np.sum(full ** 2))
    gh = full[:NG].copy()
    R, R1, R2, L = series(gh, 0.5, NBIG)
    D2 = float(np.sum((np.arange(NBIG) * L) ** 2))
    pairing = -float(np.dot(R2, R))
    Rx, R2x = R @ PSIX, R2 @ PSIX
    corr = float(np.sum(R2x * (Rx - np.log1p(Rx)) * phi) * dx)
    r1 = float(np.sqrt(np.sum(R1 ** 2)) / max(np.sqrt(np.sum(R ** 2)), 1e-300))
    # the original harness: outputs < NOUT, quadrature of the synthesized series on [-9,9]
    Rh, R1h, R2h, _ = (v[:NOUT] for v in (R, R1, R2, L))
    Rhx, R1hx, R2hx = Rh @ PSIX[:NOUT], R1h @ PSIX[:NOUT], R2h @ PSIX[:NOUT]
    harness = float(-np.sum(R2hx * np.log1p(Rhx) * phi) * dx - np.sum(R1hx ** 2 / (1 + Rhx) * phi) * dx)
    return {"pairing": pairing, "correction": corr, "F''": pairing + corr, "D2": D2, "harness": harness,
            "input tail": tail, "|R_t|/|R|": r1, "max |R| on [-9,9]": float(np.max(np.abs(Rx)))}

def main():
    if len(sys.argv) != 2:
        print("usage: positive_control.py OUT.json"); return 2
    delta = 1 / 12
    qs = [("4/7", 4 / 7), ("0.7", 0.7), ("0.8", 0.8)]
    Ns = [30, 40, 50, 80]
    out = {"delta": delta, "t": 0.5, "digits": [90, 60], "input degrees": NG, "output degrees": NBIG,
           "harness": {"output degrees": NOUT, "window": [-9, 9], "points": 7201}, "rows": [], "checks": {}, "ok": True}
    ok = True
    PSIG = hf(NG - 1, xg) * (np.exp(-xg * xg / 4) / (2 * math.pi) ** 0.25) * dxg      # quadrature (check only)
    for qname, q in qs:
        a = 1 / math.sqrt(1 - (q / 2) ** 2)
        Hg = hf(max(Ns), math.sqrt(2 * a) * xg) * (math.sqrt(2) * (2 * math.pi) ** 0.25)   # h_k(x) on [-16,16]
        for N in Ns:
            G90, H90, c0, c2 = exact_coeffs(N, a, 90)
            G60, _, _, _ = exact_coeffs(N, a, 60)
            gx = Hg[N] - c0 * np.exp(-xg * xg / 2) - c2 * (1 - xg * xg) * np.exp(-xg * xg / 2)
            scale = float(np.max(np.abs(gx)))
            e90 = evaluate(G90, scale, q, delta); e60 = evaluate(G60, scale, q, delta)
            dD = abs(e60["D2"] / e90["D2"] - 1); dF = abs(e60["F''"] / e90["F''"] - 1)
            gmax = max(abs(v) for v in G90)
            dG = float(max(abs(u - v) for u, v in zip(G60, G90)) / gmax)
            low = float(max(abs(G90[0]), abs(G90[2])))
            row = {"q": qname, "N": N, "a": a, "S (sup on grid)": scale,
                   "low-degree residual (90 digits)": float(low),
                   "F''(1/2)": e90["F''"], "D^2": e90["D2"], "F''(1/2)/D^2": e90["F''"] / e90["D2"],
                   "pairing -<R_tt,R>": e90["pairing"], "nonlinear correction on [-9,9]": e90["correction"],
                   "|R_t(1/2)| / |R(1/2)|": e90["|R_t|/|R|"], "max |R| on [-9,9]": e90["max |R| on [-9,9]"],
                   "input tail beyond degree 400": e90["input tail"],
                   "60 vs 90 digits: max coefficient difference / max coefficient": dG,
                   "60 vs 90 digits: rel. change of D^2": dD, "60 vs 90 digits: rel. change of F''": dF,
                   "harness F''(1/2)/D^2 (outputs < 240, [-9,9])": e90["harness"] / e90["D2"]}
            res = (dG < 1e-40 and dD < 1e-12 and dF < 1e-12 and e90["input tail"] < 1e-20
                   and abs(e90["correction"]) < 1e-3 * abs(e90["pairing"]))
            row["resolved"] = bool(res)
            if not res: ok = False
            if not e90["|R_t|/|R|"] < 1e-12: ok = False
            if not low < 1e-60: ok = False
            if N == 30:
                gq = PSIG @ (Hg[N] * (math.sqrt(2) * (2 * math.pi) ** 0.25) ** -1)
                hq = np.array([float(v) for v in H90[:NG]]) / (math.sqrt(2) * (2 * math.pi) ** 0.25)
                qerr = float(np.max(np.abs(gq - hq)))
                row["exact vs double quadrature, max abs difference of coefficients of h_N (normalised psi)"] = qerr
                if not qerr < 1e-12: ok = False
            out["rows"].append(row)
            fpp, rat, har = e90["F''"], row["F''(1/2)/D^2"], row["harness F''(1/2)/D^2 (outputs < 240, [-9,9])"]
            print(f"q={qname:4s} N={N:3d}  F''(1/2)={fpp:+.4e}  F''/D^2={rat:+.4f}  harness={har:+.4f}  "
                  f"D^2={e90['D2']:.4e}  corr/pair={abs(e90['correction'] / e90['pairing']):.1e}  "
                  f"60v90={dG:.1e}  tail={e90['input tail']:.1e}  resolved={res}")
    c47 = [r["F''(1/2)/D^2"] for r in out["rows"] if r["q"] == "4/7"]
    out["checks"]["q=4/7: largest F''(1/2)/D^2"] = max(c47)
    if not max(c47) < -0.4: ok = False
    for qname in ("0.7", "0.8"):
        pos = [r["N"] for r in out["rows"] if r["q"] == qname and r["F''(1/2)"] > 0]
        out["checks"][f"q={qname}: N with F''(1/2) > 0"] = pos
        if not pos: ok = False
    out["checks"]["all rows resolved"] = all(r["resolved"] for r in out["rows"])
    out["ok"] = ok
    with open(sys.argv[1], "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out["checks"], indent=1))
    print("ALL OK" if ok else "FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
