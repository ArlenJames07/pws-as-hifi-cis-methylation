#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./import_existing_results.sh [--dry-run] [--results-dir PATH]

Create collision-safe symbolic links from the legacy analysis outputs into the
numbered results layout. No bioinformatics program is executed and no existing
destination file is overwritten.
EOF
}

dry_run=false
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
results_dir="${RESULTS_DIR:-${project_dir}/results}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            dry_run=true
            shift
            ;;
        --results-dir)
            [[ $# -ge 2 ]] || { echo "Missing value for --results-dir" >&2; exit 2; }
            results_dir="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

alignment_dir="${LEGACY_ALIGNMENT_DIR:-/mnt/diskrare/arlenb/08/aligned_reads/t2t}"
small_raw_dir="${LEGACY_SMALL_RAW_DIR:-/mnt/diskrare/arlenb/08/small_variants/t2t}"
small_filtered_dir="${LEGACY_SMALL_FILTERED_DIR:-/home/rare/arlen/outputs/Variants/small_variants/SNVs_filtered}"
sv_signatures_dir="${LEGACY_SV_SIGNATURES_DIR:-/mnt/diskrare/arlenb/08/structural_variants/signatures/t2t}"
sv_calls_dir="${LEGACY_SV_CALLS_DIR:-/mnt/diskrare/arlenb/08/structural_variants/vcfs_filtered/t2t}"
phased_bam_dir="${LEGACY_PHASED_BAM_DIR:-/mnt/diskrare/arlenb/08/hiphase_results/bamfiles}"
phased_variants_dir="${LEGACY_PHASED_VARIANTS_DIR:-/mnt/diskrare/arlenb/08/hiphase_results/variants}"
cnv_dir="${LEGACY_CNV_DIR:-/home/rare/arlen/outputs/Variants/Structural_variants/hifi_cnv}"
methylation_dir="${LEGACY_METHYLATION_DIR:-/home/rare/arlen/outputs/methylation/genomes_2}"

linked=0
unchanged=0
conflicts=0
missing_sources=0

ensure_dir() {
    local directory="$1"
    if [[ "$dry_run" == false ]]; then
        mkdir -p "$directory"
    fi
}

safe_link() {
    local source_file="$1"
    local destination_dir="$2"
    local destination_file="${destination_dir}/$(basename "$source_file")"

    ensure_dir "$destination_dir"

    if [[ -L "$destination_file" ]]; then
        if [[ "$(readlink -f "$destination_file")" == "$(readlink -f "$source_file")" ]]; then
            ((unchanged += 1))
            return
        fi
        echo "CONFLICT: keeping existing symlink: $destination_file" >&2
        ((conflicts += 1))
        return
    fi
    if [[ -e "$destination_file" ]]; then
        echo "CONFLICT: keeping existing file: $destination_file" >&2
        ((conflicts += 1))
        return
    fi

    if [[ "$dry_run" == true ]]; then
        printf 'LINK\t%s\t->\t%s\n' "$destination_file" "$source_file"
    else
        ln -s "$source_file" "$destination_file"
    fi
    ((linked += 1))
}

link_flat_directory() {
    local source_dir="$1"
    local destination_dir="$2"
    local label="$3"

    if [[ ! -d "$source_dir" ]]; then
        echo "MISSING SOURCE ($label): $source_dir" >&2
        ((missing_sources += 1))
        return
    fi

    echo "Importing $label"
    while IFS= read -r -d '' source_file; do
        safe_link "$source_file" "$destination_dir"
    done < <(find "$source_dir" -maxdepth 1 -type f -print0 | sort -z)
}

sample_from_filename() {
    local filename="$1"
    local primary_prefix="${filename%%.*}"
    if [[ "$primary_prefix" =~ .*_([0-9]{3}[[:alpha:]])$ ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
    elif [[ "$primary_prefix" =~ ^([0-9]{3}[[:alpha:]])$ ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
    else
        printf '%s\n' "legacy_unassigned"
    fi
}

link_per_sample_directory() {
    local source_dir="$1"
    local destination_root="$2"
    local label="$3"
    local sample

    if [[ ! -d "$source_dir" ]]; then
        echo "MISSING SOURCE ($label): $source_dir" >&2
        ((missing_sources += 1))
        return
    fi

    echo "Importing $label"
    while IFS= read -r -d '' source_file; do
        sample="$(sample_from_filename "$(basename "$source_file")")"
        safe_link "$source_file" "${destination_root}/${sample}"
    done < <(find "$source_dir" -maxdepth 1 -type f -print0 | sort -z)
}

link_flat_directory "$alignment_dir" "${results_dir}/01_alignment" "alignment"
link_flat_directory "$small_raw_dir" "${results_dir}/02_small_variants/raw" "raw small variants"
link_flat_directory "$small_filtered_dir" "${results_dir}/02_small_variants/filtered" "filtered small variants"
link_flat_directory "$sv_signatures_dir" "${results_dir}/03_structural_variants/signatures" "structural-variant signatures"
link_flat_directory "$sv_calls_dir" "${results_dir}/03_structural_variants/calls" "structural-variant calls"
link_flat_directory "$phased_bam_dir" "${results_dir}/04_phasing" "phased BAM files"
link_flat_directory "$phased_variants_dir" "${results_dir}/04_phasing" "phased variant and statistics files"
link_per_sample_directory "$cnv_dir" "${results_dir}/05_cnv" "HiFiCNV outputs"
link_per_sample_directory "$methylation_dir" "${results_dir}/06_methylation" "haplotype methylation outputs"

printf '\nSummary: new_links=%d unchanged_links=%d conflicts_kept=%d missing_sources=%d\n' \
    "$linked" "$unchanged" "$conflicts" "$missing_sources"

if ((missing_sources > 0 || conflicts > 0)); then
    exit 1
fi
