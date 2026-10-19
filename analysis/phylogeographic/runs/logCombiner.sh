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
BURNIN=10000000
JAVA_XMX=2048m

# Edit these input files and output filename as required
INPUT=$(find "$1"/run* -type f -name "*.trees" -print)
OUTPUT="./${1}/${1}.12.combined.trees"

echo "$INPUT -> $OUTPUT"

[ -f "$BEAST_JAR" ] || { echo "ERROR: beast.jar not found at $BEAST_JAR" >&2; exit 1; }

java -Xmx${JAVA_XMX} -cp "$BEAST_JAR" dr.app.tools.logcombiner.LogCombiner \
    -burnin "$BURNIN" \
    -trees \
    $INPUT \
    "$OUTPUT"


