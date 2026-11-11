#!/usr/bin/env python3
"""Theorem A against the L^2 log-convexity baseline (Remark 1.1).

For every g in L^2(phi), Hoelder's inequality on the Hermite series gives
    ||T_sqrt(s) g||^2 <= ||g||_2^(2-E2) ||T_sqrt(r) g||^E2,    E2(r,s) = 2 log s / log r,
and for bounded g one may replace ||g||_2 by ||g||_inf.  This script tabulates the sharp exponent
E(r,s) = 2 log tau_r(s) / log lambda_r of Theorem A next to E2(r,s) for a few pairs, including the
critical level s = b(r) at r = 0.01, 2/7 and 0.8.

Checks (prints ALL OK and exits 0 only if all hold): E(r,b(r)) = 1 to 1e-25; E2(r, sqrt r) = 1 to 1e-25;
E(r,s) > E2(r,s) for every tabulated pair; and E(r,s) > E2(r,s) on a grid of 99 x 50 pairs r < s < 1.
mpmath at 30 digits.

Usage: python3 code/baseline.py OUT.json
"""
import json, sys
import mpmath as mp

mp.mp.dps = 30

def lam(r): return r / (1 + mp.sqrt(1 - r * r))
def b(r): return mp.sqrt(2 * r * (1 + r)) - r
def tau(r, s):
    a = 1 / mp.sqrt(1 - r * r)
    return (1 - a * (1 - s)) / (1 + a * (1 - s))
def E(r, s): return 2 * mp.log(tau(r, s)) / mp.log(lam(r))
def E2(r, s): return 2 * mp.log(s) / mp.log(r)

def main():
    if len(sys.argv) != 2:
        print("usage: baseline.py OUT.json"); return 2
    ok = True
    pairs = [("0.01", "b(r)", mp.mpf("0.01"), b(mp.mpf("0.01"))),
             ("2/7", "3/7", mp.mpf(2) / 7, mp.mpf(3) / 7),
             ("2/7", "b(r)=4/7", mp.mpf(2) / 7, mp.mpf(4) / 7),
             ("2/7", "5/7", mp.mpf(2) / 7, mp.mpf(5) / 7),
             ("0.8", "b(r)", mp.mpf("0.8"), b(mp.mpf("0.8")))]
    rows = []
    for rn, sn, r, s in pairs:
        e, e2 = E(r, s), E2(r, s)
        if not e > e2: ok = False
        rows.append({"r": rn, "s": sn, "s value": float(s), "E (Theorem A)": float(e), "E2 = 2 log s / log r": float(e2),
                     "E - E2": float(e - e2)})
        print(f"r={rn:5s} s={sn:9s} ({mp.nstr(s, 4):>6s})  E={mp.nstr(e, 4):>6s}  E2={mp.nstr(e2, 4):>6s}")
    crit = max(abs(E(r, b(r)) - 1) for r in [mp.mpf(k) / 100 for k in range(1, 100)])
    sq = max(abs(E2(r, mp.sqrt(r)) - 1) for r in [mp.mpf(k) / 100 for k in range(1, 100)])
    if not crit < mp.mpf("1e-25"): ok = False
    if not sq < mp.mpf("1e-25"): ok = False
    worst = None
    for k in range(1, 100):
        r = mp.mpf(k) / 100
        for j in range(1, 51):
            s = r + (1 - r) * mp.mpf(j) / 51
            d = E(r, s) - E2(r, s)
            if worst is None or d < worst: worst = d
    if not worst > 0: ok = False
    out = {"rows": rows, "max |E(r,b(r)) - 1|": float(crit), "max |E2(r,sqrt r) - 1|": float(sq),
           "smallest E - E2 on the grid": float(worst), "ok": ok}
    with open(sys.argv[1], "w") as fh:
        json.dump(out, fh, indent=1)
    print("smallest E - E2 on the grid:", mp.nstr(worst, 3))
    print("ALL OK" if ok else "FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
