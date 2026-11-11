#!/usr/bin/env bash
# Usage: ./run_all.sh            full pipeline: setup, downloads, smoke test, 3,840 designs, timing, repro check, analysis, report
#        ./run_all.sh report     re-run analysis + report from the existing results/
# Run the full pipeline on an otherwise idle machine: the timing benchmark is only meaningful that way (about 2 hours total on an M2 Pro).
# selection/selected.csv (the structure panel) is committed; delete it to re-select from live RCSB queries (results may then differ).
set -euo pipefail
cd "$(dirname "$0")"
MODE="${1:-full}"
PY="${PYTHON:-python3.11}"

[ -d .venv ] || "$PY" -m venv .venv
. .venv/bin/activate
pip install -q -r requirements.txt

if [ "$MODE" != "report" ]; then
  mkdir -p data results vendor/ProteinMPNN/vanilla_model_weights
  MPNN=https://raw.githubusercontent.com/dauparas/ProteinMPNN/main
  for f in protein_mpnn_utils.py protein_mpnn_run.py vanilla_model_weights/v_48_020.pt; do
    [ -f "vendor/ProteinMPNN/$f" ] || curl -sfL "$MPNN/$f" -o "vendor/ProteinMPNN/$f"
  done
  ESM=data/esm_if1_gvp4_t16_142M_UR50.pt
  [ -f "$ESM" ] || curl -fL --progress-bar https://dl.fbaipublicfiles.com/fair-esm/models/esm_if1_gvp4_t16_142M_UR50.pt -o "$ESM"
  echo "be4ba36edec22a9bfaa4946ff6b2815f1f19d8a3d7e0eada8b796d5a0eae9fd4  $ESM" | shasum -a 256 -c -
  echo "c9cb4a671d79604111231f8dbfc7c590e06f1197453b7a6854ac6661a642f5bd  vendor/ProteinMPNN/vanilla_model_weights/v_48_020.pt" | shasum -a 256 -c -
  [ -f selection/selected.csv ] || python selection/select_structures.py
  python scripts/smoke_test.py
  rm -rf results/raw && mkdir -p results/raw
  python scripts/run_experiment.py        # 3,840 designs; downloads the panel's PDB files; ~100 min
  python scripts/benchmark_timing.py      # idle-machine timing reference
  python scripts/retime.py                # re-measures settings slowed by machine load (sequences verified identical)
  python scripts/repro_check.py
fi
python scripts/analyze.py
python scripts/build_report.py
echo "Open report/report.html"
