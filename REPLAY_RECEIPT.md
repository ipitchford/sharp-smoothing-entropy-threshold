# Replay receipt (release 0.1.0-candidate, 11 November 2026)

Clean replay before freezing, from two fresh copies of the approved package (never the working directory), on
macOS 27.2 arm64 with CPython 3.13.5, numpy 2.5.3, scipy 1.18.1, mpmath 1.4.1 and matplotlib 3.11.2, started with

    PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.13 --with numpy==2.5.3 --with scipy==1.18.1 \
        --with mpmath==1.4.1 --with matplotlib==3.11.2 sh run_all.sh OUTDIR

| Mode | Every script | `run_all.sh` comparison | Wall time |
|---|---|---|---|
| normal | 13 steps printed `ALL OK` (figures: `figures written`), exit 0 | ALL REPRODUCED: all 14 JSON files byte-identical to `results/` | about 2 minutes |
| `PYTHONOPTIMIZE=1` (assertions removed) | the same, exit 0 | ALL REPRODUCED | about 4 minutes |

In both modes every log file and all three figures were also byte-identical to the shipped ones in `results/` and
`figures/`. The scripts use explicit checks rather than `assert`, so they stay active under optimisation.

The only later changes to the package are text and metadata: the DOI in the paper's front matter (PDF rebuilt with
pandoc 3.9 and Tectonic 0.17.0), README "How to cite" and "Status" sections, and the release files (`AI_INDEX.md`,
`REVIEW_SUMMARY.md`, `RESEARCH_GATES.json`, `refs.bib`, `LICENSE.md`, `CITATION.cff`, `.zenodo.json`, this receipt,
the CI workflow and `MANIFEST.sha256`). The code, `run_all.sh`, `results/` and `figures/` are byte-identical to the
replayed copies, so the replay was not repeated after those edits. GitHub CI reruns every script on Linux, in both
modes; there the byte comparison is informative only, because other platforms may change last digits.

This receipt records execution only. The computations are diagnostics, not proofs; the theorems rest on the written
proofs in the paper.
