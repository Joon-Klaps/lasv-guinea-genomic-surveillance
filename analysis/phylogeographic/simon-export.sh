#!/usr/bin/env bash
set -euo pipefail

# Export the files for Simon to a tarball, with:
# - top-level analysis files included by base name only
# - Environmental_files directory included recursively (preserving structure)
# - no parent-directory (../) paths in the archive
# - macOS extended attributes and AppleDouble (._*) entries excluded

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/../.." && pwd)

staging_dir="$script_dir/simon_export_staging"
archive_path="$script_dir/simon_export.tar.gz"

echo "Preparing staging directory: $staging_dir"
rm -rf "$staging_dir"
mkdir -p "$staging_dir"

# 1) Files that should appear in the tarball root by base name only (all live in phylogeographic dir)
base_files=(
    "LASV_GENBANK_20250725_curated.L.noNIG-CIV.beast2.log"
    "LASV_GENBANK_20250725_curated.L.noNIG-CIV.beast2.metadata_L.csv"
    "LASV_GENBANK_20250725_curated.L.noNIG-CIV.beast2.trees"
    "LASV_GENBANK_20250725_curated.L.noNIG-CIV.beast2.xml"
    "LASV_GENBANK_20250725_curated.S.noNIG-CIV.beast2.log"
    "LASV_GENBANK_20250725_curated.S.noNIG-CIV.beast2.metadata_S.csv"
    "LASV_GENBANK_20250725_curated.S.noNIG-CIV.beast2.trees"
    "LASV_GENBANK_20250725_curated.S.noNIG-CIV.beast2.xml"
)

for f in "${base_files[@]}"; do
    src=$(find "$script_dir/runs" -name "$f" -type f 2>/dev/null | head -1)
    if [[ -n "$src" ]]; then
        # Preserve relative path under runs/ in the staging dir
        rel="${src#$script_dir/}"
        dest_dir="$staging_dir/$(dirname "$rel")"
        mkdir -p "$dest_dir"
        echo "Copying file: $rel"
        cp -p "$src" "$dest_dir/"
    else
        echo "Warning: missing file '$f' — skipping" >&2
    fi
done

# 2) Copy Environmental_files directory recursively, preserving its internal structure
supp_dir="$script_dir/seraphim-scripts"
env_dir="$repo_root/analysis/tree-visualisation/Environmental_files"
if [[ -d "$env_dir" ]]; then
    echo "Copying directory recursively: Environmental_files"
    mkdir -p "$staging_dir/Environmental_files"
    # Use rsync if available for robust copying; fallback to cp -R otherwise
    if command -v rsync >/dev/null 2>&1; then
        rsync -a "$env_dir/" "$staging_dir/Environmental_files/"
    else
        cp -R "$env_dir/" "$staging_dir/Environmental_files/"
    fi
else
    echo "Warning: missing directory '$env_dir' — skipping" >&2
fi

# 3) Copy all *.r files from supplementary_files to the tarball root
shopt -s nullglob
r_scripts=("$supp_dir"/*.r)
for r in "${r_scripts[@]}"; do
    bname=$(basename "$r")
    echo "Copying R script: $bname"
    cp -p "$r" "$staging_dir/$bname"
done
shopt -u nullglob

# 4) Copy the Excel file to the tarball root
xlsx_src="$repo_root/data/new-guinea-data/GIN-LASV-final-ViMOP-seq-for-Leuven_V02.xlsx"
if [[ -e "$xlsx_src" ]]; then
    echo "Copying Excel file: $(basename "$xlsx_src")"
    cp -p "$xlsx_src" "$staging_dir/"
else
    echo "Warning: missing file '$xlsx_src' — skipping" >&2
fi

# 5) Copy the runs directory contents (treeannotater.slurm etc.) to a subdirectory in the tarball
runs_src="$script_dir/runs"
runs_dest="$staging_dir/runs"
if [[ -d "$runs_src" ]]; then
    echo "Copying runs directory contents"
    mkdir -p "$runs_dest"
    # Use rsync if available for robust copying; fallback to cp -R otherwise
    if command -v rsync >/dev/null 2>&1; then
        rsync -a "$runs_src/" "$runs_dest/"
    else
        cp -R "$runs_src/"* "$runs_dest/"
    fi
else
    echo "Warning: missing directory '$runs_src' — skipping" >&2
fi

# Create tar.gz from the staging directory content
echo "Creating archive: $archive_path"
rm -f "$archive_path"

# On macOS, avoid AppleDouble (._*) and xattrs by disabling copyfile metadata
export COPYFILE_DISABLE=1
tar -czf "$archive_path" -C "$staging_dir" .

echo "Done. Contents of archive root:"
tar -tzf "$archive_path" | sed -e 's/^\.\///' | awk -F/ '{print $1}' | sort -u

echo "Archive created at: $archive_path"