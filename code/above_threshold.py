#!/usr/bin/env python3
"""Two copies above the threshold: checks for the self-contained proof that concavity fails for q > 4/7.

For g = h_N (the oscillator function of level N for a = a_{q/2}) the proof uses

 (S)  the square-sum formula for ||T_sqrt(u) v||^2, v = sum_k beta_k h_k;
 (R)  the 45-degree rotation coefficients of h_N (x) h_N, hence the coefficients beta of
      u = E[h_N(G1) h_N(G2) | (G1+G2)/sqrt2], with sign (-1)^s on level 2N-2s;
 (L)  the lower bound  ||M_{1/2}||^2 >= (2/sqrt(Pi_a(q))) (2/(1+a)) (C(2N,N)/4^N) tau^(2N) / (2(2N+1));
 (U)  the upper bound  ||M_t||^2 <= (4/Delta) (N!)^2 e^(2N) N^(-2N) Lambda(t)^(2N) near t = 1/2,
      Lambda(t) = sqrt(U^2 t(1-t) - b tau) + U sqrt(t(1-t)),  U = tau + b,  Lambda(1/2) = tau;
 (P)  positivity of the entries of the 4x4 matrix B(t) near t = 1/2, i.e. tau - nu > 2b.

Part 1 (mpmath, 60 digits) checks (S), (R), (L), (U), (P) against ||M_t||^2 computed directly from
the Hermite coefficients of h_N.  Part 2 (double precision) evaluates the entropy itself for
g = eps * (admissible correction of h_N) and exhibits F(t_0) > F(1/2) and F''(1/2) > 0.
"""
import json, math, sys
import numpy as np
import mpmath as mp
mp.mp.dps = 60
out = {"ok": True, "part1": [], "part2": []}

def params(q):
    q = mp.mpf(q); r = q / 2; a = 1 / mp.sqrt(1 - r * r)
    tau = (1 - a + a * q) / (1 + a - a * q); nu = (a - 1 + a * q) / (1 + a + a * q)
    b = (a - 1) / (a + 1); Pi = (1 + a) ** 2 - a * a * q * q
    return q, r, a, tau, nu, b, Pi

def hermite_coeffs(N, a, nmax):
    """<h_N, e_n>_phi for n <= nmax from the generating function
       sqrt(2/(1+a)) exp(b z^2/2 - a w^2/(2(1+a)) + sqrt(2a) z w/(1+a))."""
    b = (a - 1) / (a + 1); cw = -a / (2 * (1 + a)); cz = mp.sqrt(2 * a) / (1 + a)
    pref = mp.sqrt(2 / (1 + a)); res = [mp.mpf(0)] * (nmax + 1)
    fN = mp.sqrt(mp.factorial(N))
    for n in range(N % 2, nmax + 1, 2):
        s = mp.mpf(0)
        for i in range(N % 2, min(N, n) + 1, 2):
            s += cz ** i / mp.factorial(i) * (b / 2) ** ((N - i) // 2) / mp.factorial((N - i) // 2) \
                 * cw ** ((n - i) // 2) / mp.factorial((n - i) // 2)
        res[n] = pref * fN * mp.sqrt(mp.factorial(n)) * s
    return res

def M_norm2(hq, t, nmax):
    lf = [mp.log(mp.factorial(n)) for n in range(2 * nmax + 1)]
    idx = [j for j in range(nmax + 1) if hq[j] != 0]
    pt = {j: mp.mpf(t) ** (mp.mpf(j) / 2) for j in idx}; p1 = {j: (1 - mp.mpf(t)) ** (mp.mpf(j) / 2) for j in idx}
    tot = mp.mpf(0)
    for n in range(0, 2 * nmax + 1, 2):
        s = mp.mpf(0)
        for j in idx:
            k = n - j
            if k < 0 or k > nmax or hq[k] == 0: continue
            s += hq[j] * hq[k] * mp.exp((lf[n] - lf[j] - lf[k]) / 2) * pt[j] * p1[k]
        tot += s * s
    return tot

def square_sum(beta, a, u):
    """(2/sqrt(Pi_a(u))) sum_m sigma^m ( sum_l gamma^l/l! sqrt((m+2l)!/m!) beta_{m+2l} )^2."""
    tau = (1 - a + a * u) / (1 + a - a * u); nu = (a - 1 + a * u) / (1 + a + a * u)
    sig = (tau + nu) / 2; gam = (nu - tau) / 4; Pi = (1 + a) ** 2 - a * a * u * u
    K = len(beta) - 1; tot = mp.mpf(0)
    for m in range(K + 1):
        I = mp.fsum(gam ** l / mp.factorial(l) * mp.sqrt(mp.factorial(m + 2 * l) / mp.factorial(m)) * beta[m + 2 * l]
                    for l in range((K - m) // 2 + 1))
        tot += sig ** m * I * I
    return 2 / mp.sqrt(Pi) * tot

def rotation_beta(N, a):
    """beta_{2N-2s} = D_{2N-2s} m_{2s}: coefficients of u = E[h_N h_N | Y] on h_{2N-2s}."""
    b = (a - 1) / (a + 1); beta = [mp.mpf(0)] * (2 * N + 1)
    for s in range(N + 1):
        D = (-1) ** s * mp.factorial(N) * mp.binomial(2 * N, N) * mp.binomial(N, s) / mp.binomial(2 * N, 2 * s) \
            / (mp.mpf(2) ** N * mp.sqrt(mp.factorial(2 * N - 2 * s) * mp.factorial(2 * s)))
        m = mp.sqrt(2 / (1 + a)) * mp.sqrt(mp.factorial(2 * s)) * b ** s / (mp.mpf(2) ** s * mp.factorial(s))
        beta[2 * N - 2 * s] = D * m
    return beta

# ---------------- Part 1
for qv in ("0.58", "0.6", "0.75", "0.9", "1"):
    q, r, a, tau, nu, b, Pi = params(qv)
    lam = r / (1 + mp.sqrt(1 - r * r))
    check_P = tau - nu > 2 * b
    gap = tau ** 2 > lam                       # equivalent to q > 4/7
    for N in (8, 16, 24):
        nmax = 260 if q < 1 else 420
        hN = hermite_coeffs(N, a, nmax)
        hq = [hN[n] * q ** (mp.mpf(n) / 2) for n in range(nmax + 1)]
        # sanity: S_N(q) from coefficients against the closed form of Fact 3.1
        SNq = mp.fsum(x * x for x in hq)
        cm = [mp.mpf(1)]
        for m in range(1, N + 1): cm.append(cm[-1] * (2 * m - 1) / (2 * m))
        closed = 2 / mp.sqrt(Pi) * mp.fsum(cm[j] * cm[N - j] * tau ** j * nu ** (N - j) for j in range(N + 1))
        M_half = M_norm2(hq, mp.mpf(1) / 2, nmax)
        beta = rotation_beta(N, a)
        M_half_formula = square_sum(beta, a, q)
        lower = 2 / mp.sqrt(Pi) * 2 / (1 + a) * mp.binomial(2 * N, N) / mp.mpf(4) ** N * tau ** (2 * N) / (2 * (2 * N + 1))
        Delta = (1 + a) * mp.sqrt(Pi); U = tau + b
        two_gam = (tau - nu) / 2                    # 2|gamma|
        t_lo, t_hi = b / (two_gam + b), two_gam / (two_gam + b)   # interval of Proposition 7.2(ii)
        rows = []
        for t in ("0.5", "0.53", "0.56", "0.6"):
            t = mp.mpf(t); w = t * (1 - t)
            Lam = mp.sqrt(U * U * w - b * tau) + U * mp.sqrt(w)
            upper = 4 / Delta * mp.factorial(N) ** 2 * mp.e ** (2 * N) / mp.mpf(N) ** (2 * N) * Lam ** (2 * N)
            Mt = M_norm2(hq, t, nmax) if t != mp.mpf("0.5") else M_half
            rows.append({"t": float(t), "inside interval": bool(t_lo < t < t_hi), "Lambda/tau": float(Lam / tau),
                         "||M_t||^2/upper": float(Mt / upper), "||M_t||^2/||M_1/2||^2": float(Mt / M_half)})
        rec = {"q": float(q), "N": N, "tau-nu-2b": float(tau - nu - 2 * b), "tau^2/lambda": float(tau ** 2 / lam),
               "interval of Prop 7.2(ii)": [float(t_lo), float(t_hi)],
               "lower bound of Prop 7.2(i) in closed form / lower": float(tau ** (2 * N) / (Delta * mp.sqrt(N) * (2 * N + 1)) / (2 / mp.sqrt(Pi) * 2 / (1 + a) * mp.binomial(2 * N, N) / mp.mpf(4) ** N * tau ** (2 * N) / (2 * (2 * N + 1)))),
               "S_N(q) coefficients/closed form - 1": float(SNq / closed - 1),
               "||M_1/2||^2 direct/square-sum formula - 1": float(M_half / M_half_formula - 1),
               "||M_1/2||^2 / tau^(2N)": float(M_half / tau ** (2 * N)),
               "lower bound / ||M_1/2||^2": float(lower / M_half), "upper": rows}
        ok = (check_P and gap and t_lo < mp.mpf("0.5") < t_hi and rec["lower bound of Prop 7.2(i) in closed form / lower"] <= 1
              and abs(SNq / closed - 1) < mp.mpf(10) ** -30 and abs(M_half / M_half_formula - 1) < mp.mpf(10) ** -25
              and lower <= M_half and all(x["||M_t||^2/upper"] <= 1 for x in rows)
              and abs(rows[0]["Lambda/tau"] - 1) < 1e-12 and all(x["Lambda/tau"] < 1 for x in rows[1:]))
        rec["ok"] = bool(ok); out["ok"] = out["ok"] and bool(ok); out["part1"].append(rec)
        print(f"q={float(q):.2f} N={N:2d} M/tau^2N={rec['||M_1/2||^2 / tau^(2N)']:.3e} lower/M={rec['lower bound / ||M_1/2||^2']:.3f} "
              f"M_t/upper=" + " ".join(f"{x['||M_t||^2/upper']:.3f}" for x in rows) + f" | {'ok' if ok else 'FAIL'}")
# square-sum formula on random coefficients (independent of the rotation formula)
rng = np.random.default_rng(5)
q, r, a, tau, nu, b, Pi = params("0.7")
worst = mp.mpf(0)
for trial in range(5):
    K = 12; beta = [mp.mpf(float(x)) for x in rng.normal(size=K + 1)]
    cols = [hermite_coeffs(k, a, 200) for k in range(K + 1)]
    direct = mp.fsum((mp.fsum(beta[k] * cols[k][n] for k in range(K + 1))) ** 2 * q ** n for n in range(201))
    worst = max(worst, abs(direct / square_sum(beta, a, q) - 1))
out["square-sum formula, random beta, max rel err"] = float(worst)
if worst > mp.mpf(10) ** -25: out["ok"] = False

# ---------------- Part 2: the entropy itself (double precision)
from scipy.special import gammaln
NG, NOUT = 500, 300
xg = np.linspace(-16, 16, 32001); dxg = xg[1] - xg[0]
x = np.linspace(-9, 9, 7201); dx = x[1] - x[0]
phi = np.exp(-x * x / 2) / math.sqrt(2 * math.pi); phig = np.exp(-xg * xg / 2) / math.sqrt(2 * math.pi)
def hf(kmax, y):
    psi = np.zeros((kmax + 1, y.size)); psi[0] = np.exp(-y * y / 4) / (2 * math.pi) ** 0.25
    psi[1] = y * psi[0]
    for k in range(1, kmax): psi[k + 1] = (y * psi[k] - math.sqrt(k) * psi[k - 1]) / math.sqrt(k + 1)
    return psi
PSIG = hf(NG - 1, xg) * np.sqrt(phig) * dxg; PSIX = hf(NOUT - 1, x) / np.sqrt(phi)
lg = gammaln(np.arange(2 * NOUT + 2) + 1.0); nn = np.arange(NOUT)
def series(gh, t):
    R = np.zeros(NOUT); R2 = np.zeros(NOUT); L = np.zeros(NOUT)
    ev = np.arange(4, NOUT, 2); m = ev // 2
    L[ev] = gh[ev] * (t ** m + (1 - t) ** m); R += L
    R2[ev] = gh[ev] * m * (m - 1) * (t ** (m - 2) + (1 - t) ** (m - 2))
    A = np.zeros(NOUT); B = np.zeros(NOUT); A2 = np.zeros(NOUT); B2 = np.zeros(NOUT); A1 = np.zeros(NOUT); B1 = np.zeros(NOUT)
    A[ev] = gh[ev] * t ** m; B[ev] = gh[ev] * (1 - t) ** m
    A1[ev] = gh[ev] * m * t ** (m - 1); B1[ev] = gh[ev] * m * (1 - t) ** (m - 1)
    A2[ev] = gh[ev] * m * (m - 1) * t ** (m - 2); B2[ev] = gh[ev] * m * (m - 1) * (1 - t) ** (m - 2)
    Mq = np.zeros(NOUT); M2 = np.zeros(NOUT)
    for n in range(8, NOUT, 2):
        j = np.arange(4, n - 3, 2); k = n - j
        w = np.exp(0.5 * (lg[n] - lg[j] - lg[k]))
        Mq[n] = np.sum(w * A[j] * B[k]); M2[n] = np.sum(w * (A2[j] * B[k] - 2 * A1[j] * B1[k] + A[j] * B2[k]))
    return R + Mq, R2 + M2, L, Mq
def F_of(gh, t):
    R = series(gh, t)[0] @ PSIX
    return -float(np.sum((1 + R) * np.log1p(R) * phi) * dx)
def Fpp_half(gh):
    R, R2, L, Mq = series(gh, 0.5)
    Rx, R2x = R @ PSIX, R2 @ PSIX
    return -float(np.sum(R2x * np.log1p(Rx) * phi) * dx), float(np.sum(L * L)), float(np.sum(Mq * Mq))
# closer to 4/7 the dip is far below double precision (about 1e-22 at q = 0.7), so only q = 1 and q = 0.9 are shown
for qv, N, eps in ((1.0, 12, 0.05), (1.0, 16, 0.05), (0.9, 16, 0.04), (0.9, 20, 0.04)):
    a = 1 / math.sqrt(1 - (qv / 2) ** 2)
    hN = math.sqrt(2) * (2 * math.pi) ** 0.25 * hf(N, math.sqrt(2 * a) * xg)[N]
    c1, c2 = np.exp(-xg * xg / 2), xg * xg * np.exp(-xg * xg / 2)
    mom = lambda f: np.array([(PSIG @ f)[0], (PSIG @ f)[2]])
    co = np.linalg.solve(np.column_stack([mom(c1), mom(c2)]), mom(hN))
    g = eps * (hN - co[0] * c1 - co[1] * c2)
    gc = PSIG @ g
    gh = np.zeros(NOUT); gh[:] = (gc[:NOUT] * qv ** (nn / 2)); gh[1::2] = 0.0; gh[:4] = 0.0
    F2, L2, M2 = Fpp_half(gh)
    F_half = F_of(gh, 0.5)
    ts = np.linspace(0.5, 0.75, 51)
    diffs = [F_of(gh, t) - F_half for t in ts]
    i = int(np.argmax(diffs))
    curve_t = np.linspace(0.25, 0.75, 101); curve = [F_of(gh, t) - F_half for t in curve_t]
    rec = {"q": qv, "N": N, "sup norm of g": float(np.max(np.abs(g))), "low-degree residual": float(max(abs(gc[0]), abs(gc[2]))),
           "F''(1/2)": F2, "||L_1/2||^2": L2, "||M_1/2||^2": M2, "t_0": float(ts[i]), "F(t_0) - F(1/2)": float(diffs[i]),
           "eps": eps, "curve_t": [float(v) for v in curve_t], "curve_F_minus_F_half": [float(v) for v in curve]}
    ok = rec["F(t_0) - F(1/2)"] > 1e-14 and F2 > 0 and rec["low-degree residual"] < 1e-9
    rec["ok"] = bool(ok); out["ok"] = out["ok"] and bool(ok); out["part2"].append(rec)
    print(f"q={qv} N={N} sup={rec['sup norm of g']:.3f} F''(1/2)={F2:+.3e} t0={rec['t_0']:.3f} F(t0)-F(1/2)={rec['F(t_0) - F(1/2)']:+.3e} | {'ok' if ok else 'FAIL'}")
json.dump(out, open(sys.argv[1], "w"), indent=1)
print("square-sum random check:", out["square-sum formula, random beta, max rel err"])
print("ALL OK" if out["ok"] else "FAILED"); sys.exit(0 if out["ok"] else 1)
