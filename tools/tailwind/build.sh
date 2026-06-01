#!/usr/bin/env bash
#
# Recompile Tailwind CSS for the lcmsPlot Shiny app.
#
# Run from this directory (tools/tailwind/) or from anywhere — the script
# resolves paths relative to its own location. Output is written to
# inst/www/lcmsPlot.css, which IS committed to git (the installed package
# ships only the precompiled CSS, never the Tailwind sources).
#
# Requires: npx (Node >= 18).

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pkg_root="$(cd "${here}/../.." && pwd)"

mkdir -p "${pkg_root}/inst/www"

cd "${here}"

npx --yes tailwindcss@3 \
    -i ./input.css \
    -o "${pkg_root}/inst/www/lcmsPlot.css" \
    --minify
