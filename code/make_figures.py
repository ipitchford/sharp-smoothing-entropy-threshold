#!/usr/bin/env python3
"""Figures for the paper.  Usage: python3 code/make_figures.py results figures

 fig1_exponent.png      the sharp exponent E(2/7, s) of Theorem A with the exponents observed along h_n
 fig2_three_copies.png  three copies: F''(0) on the segment, divided by its predicted limit, against eps
 fig3_above.png         two copies, q = 1: F(t) - F(1/2) for a small multiple of h_16 (read from above_threshold.json)
Colours are the first four slots of a validated categorical palette; every series is also labelled.
"""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res, figdir = sys.argv[1], sys.argv[2]
os.makedirs(figdir, exist_ok=True)
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})

# ---------------- figure 1
r = 2 / 7; a = 1 / math.sqrt(1 - r * r); lam = r / (1 + math.sqrt(1 - r * r)); b = math.sqrt(2 * r * (1 + r)) - r
s = np.linspace(r + 1e-4, 0.999, 600)
E = 2 * np.log((1 - a * (1 - s)) / (1 + a * (1 - s))) / math.log(lam)
obs = [(p["s"], p["observed_exponent_n=100,400,1200"][2]) for p in json.load(open(os.path.join(res, "holder.json")))["pairs"]
       if abs(p["r"] - r) < 1e-9]
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.axhline(1, color=GRID, lw=1); ax.axvline(b, color=GRID, lw=1)
ax.plot(s, E, color=C[0], lw=2, solid_capstyle="round", label="sharp exponent $E(2/7,s)$ (Theorem A)")
ax.plot([x for x, _ in obs], [y for _, y in obs], "o", ms=7, color=C[1], mec=SURF, mew=2,
        label="observed along $h_n$, $n=1200$")
ax.annotate("critical level $s=b(r)=4/7$, $E=1$", xy=(b, 1), xytext=(b + 0.03, 1.32), color=INK, fontsize=9,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.text(0.30, 0.62, "$E\\geq 1$: the linear\ninequality holds", color=INK2, fontsize=9)
ax.text(0.74, 0.86, "$E<1$: it fails", color=INK2, fontsize=9)
ax.set_xlabel("smoothing level $s$  (with $r=2/7$)"); ax.set_ylabel("exponent $E$")
ax.set_xlim(r, 1); ax.set_ylim(0, 2.05); ax.grid(axis="y", color=GRID, lw=1)
ax.legend(frameon=False, loc="upper right", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(figdir, "fig1_exponent.png"), dpi=200); plt.close(fig)

# ---------------- figure 2: two copies above the threshold (curve stored by above_threshold.py)
rec = [p for p in json.load(open(os.path.join(res, "above_threshold.json")))["part2"] if p["q"] == 1.0 and p["N"] == 16][0]
tt = np.array(rec["curve_t"]); dd = np.array(rec["curve_F_minus_F_half"]) * 1e8
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.axhline(0, color=INK2, lw=1)
ax.plot(tt, dd, color=C[0], lw=2, solid_capstyle="round")
i0 = int(np.argmax(dd))
ax.plot([0.5], [0.0], "o", ms=7, color=C[1], mec=SURF, mew=2)
ax.plot([tt[i0], 1 - tt[i0]], [dd[i0], dd[i0]], "o", ms=7, color=C[0], mec=SURF, mew=2)
ax.annotate("local minimum at $t=1/2$", xy=(0.5, 0), xytext=(0.5, -0.9 * dd[i0]), ha="center", color=INK, fontsize=9,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.annotate(f"$F(t_0)>F(1/2)$ at $t_0={tt[i0]:.2f}$", xy=(tt[i0], dd[i0]), xytext=(tt[i0] - 0.02, dd[i0] * 1.35), ha="right",
            color=INK, fontsize=9, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.set_xlabel("share $t$ of the signal carried by the first copy")
ax.set_ylabel("$F(t)-F(1/2)$   ($\\times10^{-8}$)")
ax.set_xlim(0.25, 0.75); ax.set_ylim(min(dd.min(), -1.2 * dd[i0]) * 1.05, dd[i0] * 1.7); ax.grid(axis="y", color=GRID, lw=1)
ax.text(0.255, dd[i0] * 1.5, f"two copies, $q=1$, $g=0.05\\,\\tilde h_{{16}}$, $\\|g\\|_\\infty={rec['sup norm of g']:.3f}$", color=INK2, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(figdir, "fig3_above.png"), dpi=200); plt.close(fig)

# ---------------- figure 3 (recomputed on a finer grid of eps with the functions of three_copies.py)
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "three_copies.py")).read()
exec(src[:src.index('out = {"ok": True, "cases": []}')])
g0t = corrected(np.cos(2.5 * xg), [0, 2, 4, 6]); g1t = corrected(np.cos(1.5 * xg), [0, 2, 6])
eps_grid = np.logspace(-9, -1.5, 31)
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.axhline(1, color=GRID, lw=1); ax.axhline(0, color=INK2, lw=1)
table = {}
for ci, q in enumerate((0.9, 4 / 7, 0.5, 0.2)):
    d = 0.25
    h0 = coeffs(d * g0t, q); h1 = coeffs(d * g1t, q)
    theta = -0.5 * K4(h0)[0] / h1[4]
    g = d * (g0t + theta * g1t); h = coeffs(g, q)
    pred = -4 * h[4] * K4(h)[0]
    ratio = np.array([F_and_curv(h, e, 0.0)[1] / pred for e in eps_grid])
    table[f"{q:.4f}"] = [[float(e), float(v)] for e, v in zip(eps_grid, ratio)]
    lab = {0.9: "$q=0.9$", 4 / 7: "$q=4/7$", 0.5: "$q=0.5$", 0.2: "$q=0.2$"}[q]
    ax.plot(eps_grid, np.clip(ratio, -1.6, 2), color=C[ci], lw=2, solid_capstyle="round")
    lx, ly, ha = {0.9: (1.15e-2, 0.84, "left"), 4 / 7: (9.5e-3, -1.1, "left"), 0.5: (1.5e-3, -0.75, "right"),
                  0.2: (1.6e-5, -0.75, "right")}[q]
    ax.plot([lx * (0.82 if ha == "left" else 1.22)], [ly], "s", ms=5, color=C[ci])
    ax.text(lx, ly, lab, color=INK, fontsize=9, ha=ha, va="center")
ax.set_xscale("log"); ax.set_ylim(-1.5, 1.25); ax.set_xlim(1e-9, 10 ** -1.5)
ax.set_xlabel("distance $\\varepsilon$ from the face of the simplex")
ax.set_ylabel("$F''(0)$ / predicted limit")
ax.text(1.3e-9, 1.06, "limit $-4\\hat h_4\\mathcal{K}_4>0$: convex", color=INK2, fontsize=9)
ax.text(1.3e-9, -1.38, "below zero: concave", color=INK2, fontsize=9)
ax.grid(axis="y", color=GRID, lw=1); ax.text(1.3e-9, 0.45, "three copies, amplitude $d=1/4$ (exploratory)", color=INK2, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(figdir, "fig2_three_copies.png"), dpi=200); plt.close(fig)
json.dump(table, open(os.path.join(res, "fig2_three_copies_table.json"), "w"), indent=1)
print("figures written")
