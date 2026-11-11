#!/usr/bin/env python3
"""An analytic lower bound for the radius of Theorem C (Proposition 6.6 of the paper).

Every constant of the proof at (r, s) = (2/7, 4/7) is bounded here by a closed-form expression that
is valid for ALL oscillator levels k (no truncation, no grid).  The script evaluates those closed
forms (mpmath, 50 digits), checks that each dominates the corresponding finite maximum recorded by
ingredients.py, and evaluates the absorption inequality A(delta) delta <= 1/16 at delta = 1e-8.

Usage: python3 code/radius_bound.py results/ingredients.json results/radius_bound.json
The closed forms are elementary and have large margins; evaluating them in floating point is not an
interval certificate, but each line can be checked by hand from the rounded values printed.
"""
import json, sys
import mpmath as mp
mp.mp.dps = 50
ing = json.load(open(sys.argv[1]))
r = mp.mpf(2) / 7; s = mp.mpf(4) / 7
a = 1 / mp.sqrt(1 - r * r); lam = r / (1 + mp.sqrt(1 - r * r)); kap = mp.sqrt(lam)
tau = lambda u: (1 - a + a * u) / (1 + a - a * u)
nu = lambda u: (a - 1 + a * u) / (1 + a + a * u)
Pi = lambda u: (1 + a) ** 2 - a * a * u * u
out = {"ok": True}
def chk(name, cond, val):
    out[name] = {"ok": bool(cond), "value": float(val)}
    out["ok"] = out["ok"] and bool(cond)
    print(f"{'ok  ' if cond else 'FAIL'} {name}: {mp.nstr(val, 8)}")

chk("tau_a(s) = kappa", abs(tau(s) - kap) < mp.mpf(10) ** -40, tau(s))
# Lemma 4.1 and Fact 3.2: exact closed forms
B_a = (mp.pi * a) ** 0.25 * mp.sqrt((1 + 1 / a) / 2)
vs = ((1 + a) / 2) ** 0.25
chk("B_a < 1.34", B_a < mp.mpf("1.34"), B_a)
chk("varsigma < 1.01", vs < mp.mpf("1.01"), vs)
chk("B_a dominates the finite maximum of ||h_k||_1/(||h_k||^2 (k+1)^(1/4))", B_a >= ing["cS_sup"], ing["cS_sup"])
# Lemma 5.1 with tau = kappa:  C_tau <= kappa^(-1/2) (1+x)/(1-x)^3 + (1-x)^(-2),  x = sqrt(kappa)
x = mp.sqrt(kap)
C_tau = (1 + x) / (1 - x) ** 3 / x + (1 - x) ** -2
chk("C_tau < 54", C_tau < 54, C_tau)
# Lemma 4.2:  C_L^2 <= eps0^(-2) + (2/(e eps0 log 2))^2 2^eps0,  eps0 = elasticity of nu_a at r/2
u0 = r / 2
eps0 = 2 * a * u0 / (a * a * (1 + u0) ** 2 - 1)
chk("eps0 < 1", eps0 < 1, eps0)
CL2 = eps0 ** -2 + (2 / (mp.e * eps0 * mp.log(2))) ** 2 * 2 ** eps0
chk("C_L < 2.4", mp.sqrt(CL2) < mp.mpf("2.4"), mp.sqrt(CL2))
chk("C_L dominates the finite maximum", mp.sqrt(CL2) >= mp.mpf(ing["cL_sup"]), mp.mpf(ing["cL_sup"]))
# Lemma 4.4 at s = 4/7
th = nu(s) / tau(s); Lg = mp.log(1 / th)
brk = (1 - th) ** -0.5 + (3 / (mp.e * Lg)) ** 1.5 * th ** -0.5       # (1-theta)^(-1/2) + sup (x+1)^(3/2) theta^(x/2)
U_p0 = 2 / mp.sqrt(Pi(s)) * brk                                      # p = 0, all k
chk("a^2 (1-s)^2 = 1/5", abs(a * a * (1 - s) ** 2 - mp.mpf(1) / 5) < mp.mpf(10) ** -40, a * a * (1 - s) ** 2)
U_p1 = 2 * 2 / mp.sqrt(Pi(s * mp.e ** 0.5)) * mp.e ** (s * mp.mpf(5) / 2 * a) * brk   # p = 1, k >= 1
c0 = mp.log(1 / s) / 2
U_p1_k0 = 2 / c0 ** 2 * 2 / mp.sqrt(Pi(mp.sqrt(s)))                  # p = 1, k = 0
CU2 = max(U_p0, U_p1, U_p1_k0)
chk("C_U^2 (p = 0) < 6.6", U_p0 < mp.mpf("6.6"), U_p0)
chk("C_U^2 (p = 1, k >= 1) < 64", U_p1 < 64, U_p1)
chk("C_U^2 (p = 1, k = 0) < 28", U_p1_k0 < 28, U_p1_k0)
chk("p = 0 bound dominates the finite maximum", U_p0 >= mp.mpf(ing["U0_sup"]), mp.mpf(ing["U0_sup"]))
chk("p = 1 bound dominates the finite maximum and the limit", U_p1 >= max(mp.mpf(ing["U1_sup"]), mp.mpf(ing["U1_limit"])), mp.mpf(ing["U1_limit"]))
# Upper-rounded values displayed in the proof of Proposition 6.6; each displayed line must be a valid inequality
chk("B_a < 1.3315 (displayed)", B_a < mp.mpf("1.3315"), B_a)
chk("varsigma < 1.0054 (displayed)", vs < mp.mpf("1.0054"), vs)
Ct1 = (1 + x) / (1 - x) ** 3 / x; Ct2 = (1 - x) ** -2
chk("(1+x)/(x(1-x)^3) < 46.98 (displayed)", Ct1 < mp.mpf("46.98"), Ct1)
chk("(1-x)^(-2) < 6.86 (displayed)", Ct2 < mp.mpf("6.86"), Ct2)
CL2a = eps0 ** -2; CL2b = (2 / (mp.e * eps0 * mp.log(2))) ** 2 * 2 ** eps0
chk("eps0^(-2) < 2.006 (displayed)", CL2a < mp.mpf("2.006"), CL2a)
chk("(2/(e eps0 log 2))^2 2^eps0 < 3.687 (displayed)", CL2b < mp.mpf("3.687"), CL2b)
chk("C_L^2 < 5.693 (displayed)", CL2 < mp.mpf("5.693"), CL2)
chk("C_L < 2.39 (displayed)", mp.sqrt(CL2) < mp.mpf("2.39"), mp.sqrt(CL2))
P0 = 2 / mp.sqrt(Pi(s)); one_m = (1 - th) ** -0.5; supterm = (3 / (mp.e * Lg)) ** 1.5 * th ** -0.5
chk("2/sqrt(Pi_a(s)) < 1.0233 (displayed)", P0 < mp.mpf("1.0233"), P0)
chk("(1-theta)^(-1/2) < 1.6542 (displayed)", one_m < mp.mpf("1.6542"), one_m)
chk("sup (y+1)^(3/2) theta^(y/2) < 4.7441 (displayed)", supterm < mp.mpf("4.7441"), supterm)
chk("bracket < 6.399 (displayed)", brk < mp.mpf("6.399"), brk)
chk("1.0233 * 6.399 < 6.55 (displayed)", mp.mpf("1.0233") * mp.mpf("6.399") < mp.mpf("6.55"), mp.mpf("1.0233") * mp.mpf("6.399"))
P1 = 2 / mp.sqrt(Pi(s * mp.e ** 0.5)); ex = mp.e ** (s * mp.mpf(5) / 2 * a)
chk("2/sqrt(Pi_a(s e^(1/2))) < 1.1164 (displayed)", P1 < mp.mpf("1.1164"), P1)
chk("exp(5 a s/2) < 4.441 (displayed)", ex < mp.mpf("4.441"), ex)
prod = 2 * mp.mpf("1.1164") * mp.mpf("4.441") * mp.mpf("6.399")
chk("2 * 1.1164 * 4.441 * 6.399 < 63.46 < 64 (displayed)", prod < mp.mpf("63.46"), prod)
chk("C_U^2 (p = 1, k = 0) < 27.11 (displayed)", U_p1_k0 < mp.mpf("27.11"), U_p1_k0)
# Optional sharpening (not used for the stated radius): the proof of Theorem B multiplies the p = 0 and p = 1 bounds,
# so C_U^2 may be replaced by sqrt(C_U^2(p=0) * C_U^2(p=1)).
CUgeo = mp.sqrt(mp.mpf("6.55") * mp.mpf("63.46"))
out["sqrt(6.55 * 63.46) (product of the p = 0 and p = 1 constants)"] = float(CUgeo)
# Theorem B and Section 6.5, with the rounded-up constants
C_r = 64 * 54 ** 2 * mp.mpf("2.4") * mp.mpf("1.01") * mp.mpf("1.34")
chk("C_{2/7} <= 64 * 54^2 * 2.4 * 1.01 * 1.34 < 6.1e5", C_r < mp.mpf("6.1e5"), C_r)
C1 = mp.mpf(5) / 4 * mp.mpf("6.1e5"); C2 = mp.mpf("2.97") * mp.mpf(5) / 8 * mp.mpf("6.1e5")
A = lambda d: C1 + 20 * (1 + C1 * d) + C2 * (1 + 20 * d) * (1 + C1 * d)
d0 = mp.mpf(10) ** -8
chk("A(1e-8) * 1e-8 < 0.02 < 1/16", A(d0) * d0 < mp.mpf("0.02"), A(d0) * d0)
lo, hi = mp.mpf(0), mp.mpf(1) / 12
for _ in range(200):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if A(mid) * mid <= mp.mpf(1) / 16 else (lo, mid)
out["largest delta allowed by these constants"] = float(lo)
Cg = CUgeo * 54 ** 2 * mp.mpf("2.4") * mp.mpf("1.01") * mp.mpf("1.34")
C1g = mp.mpf(5) / 4 * Cg; C2g = mp.mpf("2.97") * mp.mpf(5) / 8 * Cg
Ag = lambda d: C1g + 20 * (1 + C1g * d) + C2g * (1 + 20 * d) * (1 + C1g * d)
lo2, hi2 = mp.mpf(0), mp.mpf(1) / 12
for _ in range(200):
    mid = (lo2 + hi2) / 2
    lo2, hi2 = (mid, hi2) if Ag(mid) * mid <= mp.mpf(1) / 16 else (lo2, mid)
out["largest delta with the product constant (optional sharpening, not claimed)"] = float(lo2)
print("optional: sqrt(6.55*63.46) =", mp.nstr(CUgeo, 4), "; largest delta with it:", mp.nstr(lo2, 4))
out["indicative radius from ingredients.py (not certified)"] = ing["delta_star_indicative"]
print("largest delta allowed by the rounded analytic constants:", mp.nstr(lo, 4))
json.dump(out, open(sys.argv[2], "w"), indent=1)
print("ALL OK" if out["ok"] else "FAILED"); sys.exit(0 if out["ok"] else 1)
