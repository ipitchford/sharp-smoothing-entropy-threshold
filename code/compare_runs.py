#!/usr/bin/env python3
"""Convergence table for three_copies.py: standard against fine numerical parameters.

Usage: python3 code/compare_runs.py results/three_copies.json results/three_copies_fine.json results/three_copies_convergence.json
"""
import json, sys
a, b = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
out = {"standard": a["numerics"], "fine": b["numerics"], "cases": [], "ok": True}
worst = 0.0
matched = len(a["cases"]) == len(b["cases"]) and len(a["cases"]) > 0
for ca, cb in zip(a["cases"], b["cases"]):
    if (ca["q"], ca["d"]) != (cb["q"], cb["d"]):        # explicit check (no assert, so it also runs under python -O)
        matched = False; break
    rel = abs(cb["prediction -4 h_4 K_4"] / ca["prediction -4 h_4 K_4"] - 1)
    ratio_shift = max(abs(ra["F''(0)/prediction"] - rb["F''(0)/prediction"]) for ra, rb in zip(ca["segment"], cb["segment"]))
    worst = max(worst, rel)
    out["cases"].append({"q": ca["q"], "d": ca["d"], "sup norm of g": ca["sup norm of g"],
                         "predicted limit (standard)": ca["prediction -4 h_4 K_4"], "predicted limit (fine)": cb["prediction -4 h_4 K_4"],
                         "relative change": rel, "largest change in F''(0)/prediction over the five eps": ratio_shift,
                         "passes (standard, fine)": [ca["ok"], cb["ok"]]})
    print(f"q={ca['q']:.4f} d={ca['d']:.4f} sup={ca['sup norm of g']:.3f} limit={ca['prediction -4 h_4 K_4']:+.4e} "
          f"rel.change={rel:.2e} ratio shift={ratio_shift:.2e} pass={ca['ok']},{cb['ok']}")
out["largest relative change in the predicted limit"] = worst
out["cases matched"] = matched
out["ok"] = bool(matched and a["ok"] and b["ok"] and worst < 1e-5)
json.dump(out, open(sys.argv[3], "w"), indent=1)
print("largest relative change:", f"{worst:.3e}")
print("ALL OK" if out["ok"] else "FAILED"); sys.exit(0 if out["ok"] else 1)
