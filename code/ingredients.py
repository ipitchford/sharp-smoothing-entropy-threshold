#!/usr/bin/env python3
"""Checks of the closed-form ingredients used in the endpoint proof (q = 4/7).

Everything here is high-precision floating point (mpmath, 60 digits), not interval
arithmetic.  The script fails (non-zero exit) if any checked inequality fails.

Notation follows the note:
  r = 2/7, q = 4/7, a = (1-r^2)^(-1/2) = 7/(3 sqrt 5), kappa = (3 - sqrt 5)/2
  h_k(x) = sqrt2 exp(-a x^2/2) He_k(sqrt(2a) x)/sqrt(k!)
  S_k(u) = || T_{sqrt u} h_k ||^2_{L^2(phi)} = 2 P_k(u) / sqrt(Pi(u))
"""
import json, sys, math
import numpy as np
import mpmath as mp

mp.mp.dps = 60
r = mp.mpf(2) / 7
q = mp.mpf(4) / 7
a = 1 / mp.sqrt(1 - r * r)
kappa = (3 - mp.sqrt(5)) / 2
KMAX = 400

def tau(u): return (1 - a + a * u) / (1 + a - a * u)
def nu(u): return (a - 1 + a * u) / (1 + a + a * u)
def Pi(u): return (1 + a) ** 2 - a * a * u * u

cm = [mp.mpf(1)]
for m in range(1, KMAX + 2):
    cm.append(cm[-1] * (2 * m - 1) / (2 * m))

def P(k, u):
    t, n = tau(u), nu(u)
    return mp.fsum(cm[j] * cm[k - j] * t ** j * n ** (k - j) for j in range(k + 1))

def S(k, u):
    return 2 * P(k, u) / mp.sqrt(Pi(u))

out = {"checks": [], "ok": True}
def check(name, cond, value):
    out["checks"].append({"name": name, "ok": bool(cond), "value": str(value)})
    if not cond:
        out["ok"] = False
        print("FAIL", name, value)

# ---- 0. algebraic identities at the endpoint
check("a^2(1-r^2)=1", abs(a * a * (1 - r * r) - 1) < mp.mpf(10) ** -50, a)
check("tau(r)=nu(r)=kappa^2", abs(tau(r) - kappa ** 2) + abs(nu(r) - kappa ** 2) < mp.mpf(10) ** -50, tau(r))
check("tau(q)=kappa", abs(tau(q) - kappa) < mp.mpf(10) ** -50, tau(q))
check("0<nu(q)<tau(q)", 0 < nu(q) < tau(q), nu(q) / tau(q))
check("b(2/7)=4/7", abs(mp.sqrt(2 * r * (1 + r)) - r - q) < mp.mpf(10) ** -50, mp.sqrt(2 * r * (1 + r)) - r)

# ---- 1. closed form S_k(u) against direct two-dimensional Gauss-Hermite quadrature
from numpy.polynomial.hermite_e import hermegauss, hermeval
xg, wg = hermegauss(180)
wg = wg / math.sqrt(2 * math.pi)
af = float(a)
def h_np(k, x):
    c = np.zeros(k + 1); c[k] = 1.0
    return math.sqrt(2) * np.exp(-af * x * x / 2) * hermeval(math.sqrt(2 * af) * x, c) / math.sqrt(math.factorial(k))
worst = 0.0
for k in (0, 2, 4, 6, 10, 16):
    for u in (0.1, 2 / 7, 0.45, 4 / 7, 0.7):
        X = xg[:, None]; Z = xg[None, :]
        Y = u * X + math.sqrt(1 - u * u) * Z
        val = float(np.sum(wg[:, None] * wg[None, :] * h_np(k, X) * h_np(k, Y)))
        ref = float(S(k, mp.mpf(u)))
        worst = max(worst, abs(val - ref) / (abs(ref) * 1e-9 + 1e-13))
check("S_k closed form vs quadrature (|err| <= 1e-13 + 1e-9|ref|)", worst < 1, worst)

# ---- 2. diagonal bound P_k(q) <= C0 kappa^k (k+1)^(-1/2)   [(U_0)]
#         and weighted versions (U_{1/2}), (U_1)
def theta_derivs(k, u):
    f = lambda y: S(k, mp.e ** y)
    y0 = mp.log(u)
    return f(y0), mp.diff(f, y0, 1), mp.diff(f, y0, 2)

ks = [k for k in range(0, KMAX + 1, 2)]
U0, U12, U1, om0, om1 = [], [], [], {}, {}
for k in ks:
    s0, s1, s2 = theta_derivs(k, q)
    base = kappa ** k
    U0.append(s0 / (base * mp.mpf(k + 1) ** (-0.5)))
    U12.append(s1 / (base * mp.mpf(k + 1) ** 0.5))
    U1.append(s2 / (base * mp.mpf(k + 1) ** 1.5))
    om0[k] = mp.sqrt(s0); om1[k] = mp.sqrt(s2)
check("U_0: sup_k S_k(q)(k+1)^(1/2)/kappa^k finite", max(U0) < 2, max(U0))
check("U_1/2: sup_k theta S_k(q)/((k+1)^(1/2) kappa^k) finite", max(U12) < 5, max(U12))
check("U_1: sup_k theta^2 S_k(q)/((k+1)^(3/2) kappa^k) finite", max(U1) < 10, max(U1))
lim1 = (q * 2 * a / (1 + a - a * q) ** 2 / tau(q)) ** 2 * 2 / mp.sqrt(Pi(q)) / mp.sqrt(mp.pi * (1 - nu(q) / tau(q)))
check("U_1 values below their analytic limit * 1.01", max(U1) < lim1 * mp.mpf('1.01'), lim1)
out["U1_limit"] = str(lim1)
out["U0_sup"], out["U12_sup"], out["U1_sup"] = str(max(U0)), str(max(U12)), str(max(U1))
out["U0_last"], out["U1_last"] = str(U0[-1]), str(U1[-1])

# ---- 3. lower-bound constant  c_L(k)^2 = (k+1)^2 * int_0^1 S_k(r v) log(1/v) dv / S_k(r)   [(L1)]
cL = {}
for k in ks:
    f = lambda v: S(k, r * v) * mp.log(1 / v)
    val = mp.quad(f, mp.linspace(0, 1, 9))
    cL[k] = mp.sqrt(val * (k + 1) ** 2 / S(k, r))
check("L1: sup_k c_L(k) finite", max(cL.values()) < 3, max(cL.values()))
out["cL_sup"] = str(max(cL.values())); out["cL_last"] = str(cL[KMAX])
check("S_k(r) exact = 2 kappa^(2k)/sqrt(2(1+a))",
      max(abs(S(k, r) / (2 * kappa ** (2 * k) / mp.sqrt(2 * (1 + a))) - 1) for k in ks) < mp.mpf(10) ** -45, 0)

# ---- 4. L^1(dx) norms of h_k  (c_S(k) = ||h_k||_1 / (||h_k||_2^2 (k+1)^(1/4)))
#         orthonormal Hermite functions by the stable three-term recurrence
def hermite_functions(kmax, y):
    """psi_k(y) = He_k(y) exp(-y^2/4) / sqrt(k! sqrt(2 pi)), orthonormal in L^2(dy)."""
    psi = np.zeros((kmax + 1, y.size))
    psi[0] = np.exp(-y * y / 4) / (2 * math.pi) ** 0.25
    if kmax >= 1:
        psi[1] = y * psi[0]
    for k in range(1, kmax):
        psi[k + 1] = (y * psi[k] - math.sqrt(k) * psi[k - 1]) / math.sqrt(k + 1)
    return psi
y = np.linspace(-70, 70, 280001)
psi = hermite_functions(KMAX, y)
dy = y[1] - y[0]
# h_k(x) = sqrt2 (2 pi)^(1/4) psi_k(sqrt(2a) x);  dx = dy / sqrt(2a)
scale = math.sqrt(2) * (2 * math.pi) ** 0.25
L2sq = 2 * math.sqrt(math.pi / af)
cS = {}
worst_norm = 0.0
for k in ks:
    l1 = scale * np.sum(np.abs(psi[k])) * dy / math.sqrt(2 * af)
    l2 = scale ** 2 * np.sum(psi[k] ** 2) * dy / math.sqrt(2 * af)
    worst_norm = max(worst_norm, abs(l2 / L2sq - 1))
    cS[k] = l1 / (L2sq * (k + 1) ** 0.25)
check("||h_k||_2^2 = 2 sqrt(pi/a) (rel err)", worst_norm < 1e-9, worst_norm)
rig = math.sqrt(math.pi) * math.sqrt(1 + 1 / af) / math.sqrt(L2sq)   # B_a of Lemma 3.2
check("c_S(k) <= rigorous Cauchy-Schwarz constant", max(cS.values()) <= rig, rig)
out["cS_sup"] = max(cS.values()); out["cS_rigorous_bound"] = rig
out["sup_h_k_over_k"] = float(max(scale * np.max(np.abs(psi[k])) for k in ks))

# ---- 5. the two-sided summation: constants in (E0), (E1)
#   |c_k| <= min(alpha_k D1, beta_k delta),
#   alpha_k = (5/4) c_L(k) (k+1)^(-1) S_k(r)^(-1/2),  beta_k = ||h_k||_1/||h_k||_2^2
alpha = {k: mp.mpf(5) / 4 * cL[k] / (k + 1) / mp.sqrt(S(k, r)) for k in ks}
beta = {k: mp.mpf(cS[k]) * mp.mpf(k + 1) ** 0.25 for k in ks}
def W(x):
    w0 = mp.fsum(min(alpha[k] * x, beta[k]) * om0[k] for k in ks)
    w1 = mp.fsum(min(alpha[k] * x, beta[k]) * om1[k] for k in ks)
    return w0, w1
E0, E1, grid = [], [], []
for e in np.arange(-120.0, 1.01, 0.25):
    x = mp.mpf(10) ** e
    w0, w1 = W(x)
    E0.append(w0 * w0 / x); E1.append(w0 * w1 / x); grid.append(e)
CE0, CE1 = max(E0), max(E1)
check("E0 constant finite on grid", CE0 < 50, CE0)
check("E1 constant finite on grid", CE1 < 200, CE1)
check("E1 profile decreasing at small x (log gain)", E1[0] < E1[len(E1) // 2], (E1[0], E1[len(E1) // 2]))
out["C_E0"], out["C_E1"] = float(CE0), float(CE1)
out["E1_profile"] = [[float(g), float(v)] for g, v in list(zip(grid, E1))[::40]]

# ---- 5b. plain endpoint lemma (all k, no parity, no weights):
#   ||T_sqrt(q) g||^2 <= J_end ||g||_inf ||T_sqrt(r) g||,  |c_k| <= min(S_k(r)^(-1/2) D0, beta_k delta)
allk = list(range(0, KMAX + 1))
om0_all = {k: mp.sqrt(S(k, q)) for k in allk}
beta_all = {}
for k in allk:
    l1 = scale * np.sum(np.abs(psi[k])) * dy / math.sqrt(2 * af)
    beta_all[k] = mp.mpf(l1 / L2sq)
a0 = {k: 1 / mp.sqrt(S(k, r)) for k in allk}
Jg = []
for e in np.arange(-120.0, 1.01, 0.25):
    x = mp.mpf(10) ** e
    w = mp.fsum(min(a0[k] * x, beta_all[k]) * om0_all[k] for k in allk)
    Jg.append(w * w / x)
out["J_end"] = float(max(Jg))
check("plain endpoint constant finite on grid", max(Jg) < 50, max(Jg))

out["per_k"] = {str(k): {"cL": float(cL[k]), "beta": float(beta[k]), "S_r": mp.nstr(S(k, r), 20),
                          "om0": mp.nstr(om0[k], 20), "om1": mp.nstr(om1[k], 20)} for k in ks}

# ---- 6. indicative radius, uniform in t in [0,1] and in q in (0, 4/7]
#   ||N M||        <= 2 ||N h|| ||h||      <= 2 C_E1 delta D_1 <= C1 delta D,   C1 = C_E1   (D >= 2 D_1)
#   ||N^-1 M_tt||  <= 2.97 ||N h|| ||h||   <= C2 delta D,                        C2 = 2.97 C_E1 / 2
#   F'' <= -D^2/8 + [C1 + 20 (1 + C1 d) + C2 (1 + 20 d)(1 + C1 d)] d D^2
C1 = float(CE1)
C2 = 2.97 * float(CE1) / 2
def A_in(d): return C1 + 20 * (1 + C1 * d) + C2 * (1 + 20 * d) * (1 + C1 * d)
lo, hi = 0.0, 1 / 12
for _ in range(200):
    mid = (lo + hi) / 2
    if A_in(mid) * mid <= 1 / 16: lo = mid
    else: hi = mid
out["delta_star_indicative"] = lo; out["C1"] = C1; out["C2"] = C2
# log gain at the endpoint: J(x) = (level sum)^2 / x as a function of x = ||T_sqrt(r) g|| / ||g||_inf
out["plain_profile_J(x)"] = [[float(e), float(v)] for e, v in list(zip(np.arange(-120.0, 1.01, 0.25), Jg))[::40]]
check("plain endpoint profile decreases as x -> 0 (log gain)", Jg[0] < Jg[len(Jg) // 2] < max(Jg), (Jg[0], max(Jg)))

json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "ingredients.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("checks", "E1_profile", "per_k", "plain_profile_J(x)")}, indent=1))
print("ALL OK" if out["ok"] else "FAILED")
sys.exit(0 if out["ok"] else 1)
