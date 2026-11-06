#!/usr/bin/env bash
# Build the release zip of the Compounds Chromatograms Compound Discoverer node.
#
#   tools/cd-node/build-zip.sh [output directory, default dist/]
#
# The zip holds the node folder (with node-info.json, which records the node
# version, the lcmsPlot version its installer sets up and the commit) and the
# installer: install.bat, install.ps1, install-packages.R and INSTALL.txt.
# Writes CompoundsChromatograms-<node version>.zip and SHA256SUMS.txt, and
# prints the node and lcmsPlot versions.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
OUT=$(mkdir -p "${1:-$ROOT/dist}" && cd "${1:-$ROOT/dist}" && pwd)
NODE_DIR="$ROOT/inst/compound-discoverer-nodes/CompoundsChromatograms"
TOOLS="$ROOT/tools/cd-node"

node_version=$(tr -d '[:space:]' < "$NODE_DIR/VERSION")
lcmsplot_version=$(sed -n 's/^Version: *//p' "$ROOT/DESCRIPTION" | tr -d '[:space:]')
commit=${GITHUB_SHA:-$(git -C "$ROOT" rev-parse HEAD)}
for v in "$node_version" "$lcmsplot_version"; do
  [[ $v =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || { echo "Invalid version: '$v'" >&2; exit 1; }
done

stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
mkdir "$stage/CompoundsChromatograms"
cp "$NODE_DIR"/* "$stage/CompoundsChromatograms/"
cp "$TOOLS/install.bat" "$TOOLS/install.ps1" "$TOOLS/install-packages.R" \
   "$TOOLS/INSTALL.txt" "$stage/"

cat > "$stage/CompoundsChromatograms/node-info.json" <<EOF
{
  "node_version": "$node_version",
  "lcmsplot_version": "$lcmsplot_version",
  "commit": "$commit"
}
EOF

# Windows line endings for the files read on Windows.
for f in "$stage"/install.bat "$stage"/install.ps1 "$stage"/INSTALL.txt; do
  perl -pi -e 's/\r?\n/\r\n/' "$f"
done

zip_name="CompoundsChromatograms-$node_version.zip"
rm -f "$OUT/$zip_name"
(cd "$stage" && zip -qrX "$OUT/$zip_name" .)
(cd "$OUT" && shasum -a 256 "$zip_name" > SHA256SUMS.txt)

echo "node_version=$node_version"
echo "lcmsplot_version=$lcmsplot_version"
echo "zip=$OUT/$zip_name"
