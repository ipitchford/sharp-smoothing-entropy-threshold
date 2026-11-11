#!/bin/sh
# Runs every diagnostic and the figure script into a fresh directory, then compares the regenerated
# JSON files with the ones shipped in results/.
#   sh run_all.sh OUTDIR          (about five minutes on one core of the machine in ENVIRONMENT.txt)
set -e
OUT=${1:?usage: sh run_all.sh OUTDIR}
[ -e "$OUT" ] && { echo "$OUT already exists; give a fresh directory"; exit 1; }
mkdir -p "$OUT"
step() { name=$1; shift; echo "== $name"; "$@" > "$OUT/$name.log" 2>&1 || { echo "FAILED: $name (see $OUT/$name.log)"; exit 1; }; tail -n 1 "$OUT/$name.log"; }
step ingredients           python3 code/ingredients.py "$OUT/ingredients.json"
step radius_bound          python3 code/radius_bound.py "$OUT/ingredients.json" "$OUT/radius_bound.json"
step holder                python3 code/holder.py "$OUT/holder.json"
step baseline              python3 code/baseline.py "$OUT/baseline.json"
step best_constant         python3 code/best_constant.py "$OUT/best_constant.json"
step adversarial           python3 code/adversarial.py "$OUT/ingredients.json" "$OUT/adversarial.json"
step bilinear              python3 code/bilinear.py "$OUT/bilinear.json"
step entropy_curvature     python3 code/entropy_curvature.py "$OUT/ingredients.json" "$OUT/entropy_curvature.json"
step positive_control      python3 code/positive_control.py "$OUT/positive_control.json"
step above_threshold       python3 code/above_threshold.py "$OUT/above_threshold.json"
step three_copies          python3 code/three_copies.py "$OUT/three_copies.json"
step three_copies_fine     python3 code/three_copies.py "$OUT/three_copies_fine.json" --fine
step three_copies_convergence python3 code/compare_runs.py "$OUT/three_copies.json" "$OUT/three_copies_fine.json" "$OUT/three_copies_convergence.json"
step figures               python3 code/make_figures.py "$OUT" "$OUT/figures"
echo "== comparison with results/"
status=0
for f in "$OUT"/*.json; do
  b=$(basename "$f")
  if [ -f "results/$b" ]; then
    if cmp -s "$f" "results/$b"; then echo "identical  $b"; else echo "DIFFERS    $b"; status=1; fi
  else echo "new        $b"; fi
done
[ $status -eq 0 ] && echo "ALL REPRODUCED" || echo "SOME FILES DIFFER (floating-point libraries of another version can change last digits)"
exit $status
