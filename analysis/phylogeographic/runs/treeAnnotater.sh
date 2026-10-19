#!/usr/bin/env bash
set -euo pipefail

# Minimal editable LogCombiner wrapper
# Edit variables at top. If you want the script to attempt locating beast.jar you can call:
#   BEAST_JAR="$(find_beast_jar)"  or export BEAST_JAR before running.
find_beast_jar() {
  # 1) BEAST_JAR env
  [ -n "${BEAST_JAR:-}" ] && [ -f "$BEAST_JAR" ] && { echo "$BEAST_JAR"; return 0; }
  # 2) BEAST env (official wrapper style)
  if [ -n "${BEAST:-}" ] && [ -f "$BEAST/lib/beast.jar" ]; then
    echo "$BEAST/lib/beast.jar"; return 0
  fi
  # 3) 'beast' on PATH and resolve symlinks
  if command -v beast >/dev/null 2>&1; then
    PRG=$(command -v beast)
    while [ -h "$PRG" ]; do
      ls=$(ls -ld "$PRG")
      link=$(expr "$ls" : '.*-> \(.*\)$')
      if expr "$link" : '/.*' >/dev/null; then
        PRG="$link"
      else
        PRG="$(dirname "$PRG")/$link"
      fi
    done
    BEAST0="$(dirname "$PRG")/.."
    BEAST_HOME="$(cd "$BEAST0" && pwd 2>/dev/null || true)"
    [ -f "$BEAST_HOME/lib/beast.jar" ] && { echo "$BEAST_HOME/lib/beast.jar"; return 0; }
  fi
  # 4) common locations
  for p in "/usr/local/lib/beast.jar" "/usr/local/Cellar/beast"/*/lib/beast.jar "$HOME/bin/BEAST_X"*/lib/beast.jar; do
    [ -f "$p" ] && { echo "$p"; return 0; }
  done
  return 1
}

# Default values (edit as needed)
BEAST_JAR="${BEAST_JAR:-/Users/joonklaps/bin/BEAST_X_v10.5.0-beta3/lib/beast.jar}"  # change this to your beast.jar
JAVA_XMX=2048m

# Edit these input files and output filename as required
INPUT=$(find "./$1" -name "*combined.trees")
OUTPUT="./${1}/${1}.12.combined.ann.tree"
type=hipstr

echo "$INPUT -> $OUTPUT"

[ -f "$BEAST_JAR" ] || { echo "ERROR: beast.jar not found at $BEAST_JAR" >&2; exit 1; }

java -Xmx${JAVA_XMX} -cp "$BEAST_JAR" dr.app.tools.treeannotator.TreeAnnotator \
  --help
    # -type $type \
    # -heights mean \
    # $INPUT \
    # "$OUTPUT"




# #!/bin/bash
# #SBATCH --clusters=genius
# #SBATCH --account="lp_phylogeo_inf_gpu"
# #SBATCH --mail-type="END,FAIL,TIME_LIMIT"
# #SBATCH --mail-user="joon.klaps@kuleuven.be"
# #SBATCH --partition=batch
# #SBATCH --exclusive
# #SBATCH --time="24:00:00"

# abort() {
#     echo >&2 '
# ***************
# *** ABORTED ***
# ***************
# '
#     echo "An error occurred. Exiting..." >&2
#     exit 1
# }

# set -e
# trap 'abort' 0

# hostname >&2
# echo $0 >&2
# printf START >&2; uptime >&2
# date >&2

# module load cluster/genius/gpu_p100_long
# module load CUDA
# module load Beast/10.5.0-GCC-12.3.0-CUDA-12.1.1

# java -Xmx180g -cp /apps/leuven/rocky8/skylake/2023a/software/Beast/10.5.0-GCC-12.3.0-CUDA-12.1.1/lib/beast.jar dr.app.tools.treeannotator.TreeAnnotator LASV.S.combined1-5.trees LASV.S.combined1-5.tree

# printf END >&2; uptime >&2

# trap : 0

# echo >&2 '
# ************
# *** DONE ***
# ************
# '