#!/usr/bin/env bash
# Run the whole demo end to end, teeing every step into demo/logs/.
# Usage: bash demo/run_all.sh              (build + synthetic + scaling + real)
#        SKIP_SCALING=1 bash demo/run_all.sh   (skip the scaling sweep)
#        SKIP_REAL=1 bash demo/run_all.sh      (synthetic only - no downloads)
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
LOGS="$DEMO_DIR/logs"; mkdir -p "$LOGS"

bash "$DEMO_DIR/02_build_g2g.sh"        2>&1 | tee "$LOGS/02_build_g2g.log"
bash "$DEMO_DIR/03_synthetic_demo.sh"   2>&1 | tee "$LOGS/03_synthetic_demo.log"

if [ "${SKIP_SCALING:-0}" != "1" ]; then
    bash "$DEMO_DIR/05_scaling_demo.sh" 2>&1 | tee "$LOGS/05_scaling_demo.log"
fi

if [ "${SKIP_REAL:-0}" = "1" ]; then
    echo "SKIP_REAL=1 - skipping the real-data demo."
    exit 0
fi
bash "$DEMO_DIR/01_setup_tools.sh"      2>&1 | tee "$LOGS/01_setup_tools.log"
bash "$DEMO_DIR/04_real_hla_demo.sh"    2>&1 | tee "$LOGS/04_real_hla_demo.log"
echo; echo "All demos complete. Logs in $LOGS"
