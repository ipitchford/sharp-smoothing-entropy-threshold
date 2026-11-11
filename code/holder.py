#!/usr/bin/env python3
"""Checks for the sharp exponent E(r,s) of the bounded-input smoothing inequality

    ||T_sqrt(s) g||^2 <= C ||g||_inf^(2-E) ||T_sqrt(r) g||^E,      0 < r < s < 1,
    E(r,s) = 2 log tau_r(s) / log lambda_r   (between 0 and 2; equal to 1 at s = b(r)).

 (1) algebra: tau_r(b(r)) = sqrt(lambda_r) for a grid of r (so E >= 1 exactly when s <= b(r));
 (2) optimality: along g = h_{n,a_r} the observed exponent log N_n / log D_n tends to E(r,s)
     from above, where N_n = ||T_sqrt(s) h_n||^2, D_n = ||T_sqrt(r) h_n||;
 (3) the inequality itself on about 300 bounded test functions for several (r,s);
 (4) monotonicity of Theta -> log tanh(Theta+d)/log tanh(Theta) (not used in the proof; recorded
     because it explains why no other oscillator parameter does better).
mpmath (50 digits) for (1),(2),(4); double precision for (3).  Diagnostics, not proofs.
"""
import json, math, sys
import numpy as np
import mpmath as mp
mp.mp.dps = 50
out = {"ok": True, "pairs": []}

def lam(r): return r / (1 + mp.sqrt(1 - r * r))
def b(r): return mp.sqrt(2 * r * (1 + r)) - r
def a_of(r): return 1 / mp.sqrt(1 - r * r)
def tau(a, u): return (1 - a + a * u) / (1 + a - a * u)
def nu(a, u): return (a - 1 + a * u) / (1 + a + a * u)
def Pi(a, u): return (1 + a) ** 2 - a * a * u * u
def E(r, s):
    return 2 * mp.log(tau(a_of(r), s)) / mp.log(lam(r))

# (1)
worst = max(abs(tau(a_of(r), b(r)) ** 2 - lam(r)) for r in [mp.mpf(i) / 100 for i in range(1, 100)])
out["max |tau_r(b(r))^2 - lambda_r| over r=0.01..0.99"] = mp.nstr(worst, 5)
if worst > mp.mpf(10) ** -40: out["ok"] = False

# (2)
cm = [mp.mpf(1)]
for m in range(1, 1300): cm.append(cm[-1] * (2 * m - 1) / (2 * m))
def S(k, a, u):
    t, n = tau(a, u), nu(a, u)
    return 2 * mp.fsum(cm[j] * cm[k - j] * t ** j * n ** (k - j) for j in range(k + 1)) / mp.sqrt(Pi(a, u))
pairs = []
for r in (mp.mpf(1) / 10, mp.mpf(2) / 7, mp.mpf(1) / 2, mp.mpf(4) / 5):
    for s in ((3 * r + b(r)) / 4, (r + b(r)) / 2, b(r), (2 * b(r) + 1) / 3, (b(r) + 2) / 3, mp.mpf(99) / 100):
        pairs.append((r, s))
for r, s in pairs:
    a = a_of(r); e = E(r, s); obs = []
    for n in (100, 400, 1200):
        N = S(n, a, s); D = mp.sqrt(S(n, a, r))
        obs.append(float(mp.log(N) / mp.log(D)))
    rec = {"r": float(r), "s": float(s), "b(r)": float(b(r)), "E": float(e), "observed_exponent_n=100,400,1200": obs}
    # the observed exponent must approach E from above as n grows (N_n carries an extra n^(-1/2))
    if not (obs[0] >= obs[1] >= obs[2] >= float(e) - 1e-12 and obs[2] - float(e) < 0.02 * max(float(e), 0.05) + 0.01):
        out["ok"] = False; rec["FAIL"] = True
    out["pairs"].append(rec)

# (4)
f = lambda th, d: mp.log(mp.tanh(th + d)) / mp.log(mp.tanh(th))
mono = all(f(th + mp.mpf('0.01'), d) > f(th, d) for d in (mp.mpf('0.05'), mp.mpf('0.5'), 2)
           for th in [mp.mpf(i) / 50 for i in range(1, 200)])
out["ratio log tanh(Theta+d)/log tanh(Theta) increasing on grid"] = bool(mono)
if not mono: out["ok"] = False

# (3)
rng = np.random.default_rng(11)
x = np.linspace(-16, 16, 32001); dx = x[1] - x[0]; NMAX = 700
def hf(kmax, y):
    psi = np.zeros((kmax + 1, y.size)); psi[0] = np.exp(-y * y / 4) / (2 * math.pi) ** 0.25
    if kmax >= 1: psi[1] = y * psi[0]
    for k in range(1, kmax): psi[k + 1] = (y * psi[k] - math.sqrt(k) * psi[k - 1]) / math.sqrt(k + 1)
    return psi
PSI = hf(NMAX - 1, x) * (np.exp(-x * x / 4) / (2 * math.pi) ** 0.25) * dx
n = np.arange(NMAX)
tests = []
for (r, s) in [(2 / 7, 0.35), (2 / 7, 0.45), (2 / 7, 4 / 7), (2 / 7, 0.7), (2 / 7, 0.9), (0.1, 0.2), (0.1, float(b(mp.mpf(1) / 10))), (0.1, 0.6), (0.5, 0.6), (0.5, 0.85), (0.5, 0.95)]:
    a = float(a_of(mp.mpf(r))); e = float(E(mp.mpf(r), mp.mpf(s)))
    H = hf(80, math.sqrt(2 * a) * x)
    fam = [H[k] for k in range(0, 81, 2)] + [np.sign(H[k]) for k in range(0, 81, 2)]
    fam += [np.cos(w * x) for w in np.linspace(0.3, 12, 40)]
    fam += [(np.abs(x - x0) < w).astype(float) for x0 in (0, 1.5) for w in (0.01, 0.1, 1)]
    for _ in range(150):
        m = int(rng.integers(1, 7)); idx = rng.integers(0, 81, size=m); co = rng.normal(size=m) * rng.uniform(0.3, 30)
        fam.append(np.tanh(sum(c * H[i] for c, i in zip(co, idx))))
    best = 0.0
    for g in fam:
        G = PSI @ g; sup = float(np.max(np.abs(g)))
        N = float(np.sum(s ** n * G * G)); D = math.sqrt(float(np.sum(r ** n * G * G)))
        if D < 1e-140: continue
        best = max(best, N / (sup ** (2 - e) * D ** e))
    tests.append({"r": r, "s": s, "E": e, "functions": len(fam), "max N/(sup^(2-E) D^E)": best})
    if best > 50: out["ok"] = False
out["inequality_tests"] = tests
# (5) the constant in the linear inequality at the critical level is not one:
#     g = sign(x) has ||T_sqrt(u) g||^2 = (2/pi) arcsin(u) (Sheppard), so at s = b(r)
#     ratio(r) = (2/pi) arcsin(b(r)) / sqrt((2/pi) arcsin(r)),  which tends to 2/sqrt(pi) = 1.128 as r -> 0
sgn = []
for rr in (mp.mpf('0.0001'), mp.mpf('0.001'), mp.mpf('0.01'), mp.mpf('0.03'), mp.mpf('0.05'), mp.mpf('0.1'), mp.mpf(2) / 7):
    ratio = (2 / mp.pi) * mp.asin(b(rr)) / mp.sqrt((2 / mp.pi) * mp.asin(rr))
    sgn.append({"r": float(rr), "b(r)": float(b(rr)), "ratio for g = sign(x)": float(ratio)})
out["sign_function_constant"] = sgn; out["limit 2/sqrt(pi)"] = float(2 / mp.sqrt(mp.pi))
if not (sgn[2]["ratio for g = sign(x)"] > 1.05 and abs(sgn[0]["ratio for g = sign(x)"] - float(2 / mp.sqrt(mp.pi))) < 0.02): out["ok"] = False
# quadrature cross-check of Sheppard's formula at r = 0.01, s = b(r), and an even mean-zero example
sg = np.sign(x); G = PSI @ sg
r0 = 0.01; s0 = float(b(mp.mpf(r0)))
chk = float(np.sum(s0 ** n * G * G)) / math.sqrt(float(np.sum(r0 ** n * G * G)))
c0 = 0.6744897501960817                                   # median of |X|
ge = np.sign(x * x - c0 * c0); Ge = PSI @ ge
ev = float(np.sum(s0 ** n * Ge * Ge)) / math.sqrt(float(np.sum(r0 ** n * Ge * Ge)))
out["sign(x) ratio by quadrature at r=0.01"] = chk; out["even mean-zero sign(x^2-c^2) ratio at r=0.01"] = ev
if abs(chk / sgn[2]["ratio for g = sign(x)"] - 1) > 2e-3 or ev <= 1.05: out["ok"] = False
for row in sgn: print(row)
print("quadrature:", chk, " even example:", ev)

json.dump(out, open(sys.argv[1], "w"), indent=1)
for p in out["pairs"]: print(p)
for t in tests: print(t)
print({k: v for k, v in out.items() if k not in ("pairs", "inequality_tests")})
print("ALL OK" if out["ok"] else "FAILED"); sys.exit(0 if out["ok"] else 1)
