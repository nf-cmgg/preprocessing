#!/usr/bin/env bash
set -euo pipefail

export NXF_DISABLE_PIPELINE_PLUGINS=true
export NXF_VER=25.10.4

if [[ $# -eq 0 ]]; then
  echo "Usage: tests/run-nf-test-local.sh <nf-test args>"
  echo "Example: tests/run-nf-test-local.sh test tests/modules/local/umi_consensus/main.nf.test"
  exit 1
fi

nf-test "$@"
