#!/bin/bash
# Export SimCRC slides to PDF without sidebar clipping
# Usage: ./export_pdf.sh [qmd_file]
#   Default: calibration_results_Adenoma_F_v0.13.0.qmd

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
QMD_FILE="${1:-calibration_results_Adenoma_F_v0.13.0.qmd}"
QMD_PATH="$SCRIPT_DIR/$QMD_FILE"
PDF_FILE="${QMD_FILE%.qmd}.pdf"
HTML_FILE="${QMD_FILE%.qmd}.html"

if [ ! -f "$QMD_PATH" ]; then
  echo "Error: $QMD_PATH not found"
  exit 1
fi

echo "==> Switching to PDF CSS..."
sed -i '' 's/css: simcrc\.css/css: [simcrc.css, simcrc-pdf.css]/' "$QMD_PATH"

echo "==> Rendering HTML..."
quarto render "$QMD_PATH"

echo "==> Exporting PDF with decktape..."
decktape reveal "$SCRIPT_DIR/$HTML_FILE" "$SCRIPT_DIR/$PDF_FILE" --size 2800x1800

echo "==> Restoring original CSS..."
sed -i '' 's/css: \[simcrc\.css, simcrc-pdf\.css\]/css: simcrc.css/' "$QMD_PATH"

echo "==> Re-rendering HTML with sidebar..."
quarto render "$QMD_PATH"

echo "==> Done! PDF saved to: $SCRIPT_DIR/$PDF_FILE"
