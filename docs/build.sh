#!/bin/bash
# Render a findings document to Word and PDF.
#
# There is no LaTeX on this machine, so the PDF goes through weasyprint rather
# than pdflatex. Both need the repository root on the resource path, because the
# figures are referenced as ../charts/inv009/ from inside docs/.
set -euo pipefail
cd "$(dirname "$0")/.."
SRC="${1:?usage: docs/build.sh docs/<name>.md [output-stem]}"
STEM="${2:-$(basename "${SRC%.md}")}"
OUT="docs"

pandoc "$SRC" -o "$OUT/$STEM.docx" \
  --resource-path=.:docs:charts/inv009 --standalone

pandoc "$SRC" -o "$OUT/$STEM.html" \
  --resource-path=.:docs:charts/inv009 --standalone
weasyprint "$OUT/$STEM.html" "$OUT/$STEM.pdf" --stylesheet docs/doc.css
rm -f "$OUT/$STEM.html"
echo "built $OUT/$STEM.docx and $OUT/$STEM.pdf"
