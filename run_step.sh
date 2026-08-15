#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 {alignment|small_variants|structural_variants|phasing|cnv|methylation}" >&2
    exit 2
fi

stage="$1"
case "$stage" in
    alignment|small_variants|structural_variants|phasing|cnv|methylation) ;;
    *)
        echo "Invalid stage: $stage" >&2
        echo "Choose: alignment, small_variants, structural_variants, phasing, cnv, or methylation" >&2
        exit 2
        ;;
esac

params_file="${PARAMS_FILE:-params.local.yml}"
profile="${NXF_PROFILE:-workstation,docker}"

if [[ ! -f "$params_file" ]]; then
    echo "Parameter file not found: $params_file" >&2
    exit 2
fi

exec nextflow run main.nf \
    -profile "$profile" \
    -params-file "$params_file" \
    --stage "$stage" \
    -resume \
    -ansi-log false
