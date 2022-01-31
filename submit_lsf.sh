#!/usr/bin/env bash
set -eux

JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')
LOG_DIR=logs/

if [[ ! -d "$LOG_DIR" ]]; then
  echo "Error: Log directory $LOG_DIR does not exist"
  exit 1
fi

MEMORY=4000
THREADS=4
PROFILE="lsf"
SINGULARITY_ARGS="'--nv'"

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
  -q long \
  -M "$MEMORY" \
  -n "$THREADS" \
  -o "$LOG_DIR"/"$JOB_NAME".o \
  -e "$LOG_DIR"/"$JOB_NAME".e \
  -J "$JOB_NAME" \
  snakemake --profile "$PROFILE" \
  --local-cores "$THREADS" \
  "$@" \
  --singularity-args "$SINGULARITY_ARGS"

exit 0
