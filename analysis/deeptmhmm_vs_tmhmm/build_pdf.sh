#!/bin/bash
# Build REPORT.pdf from REPORT.md (pandoc + TeX Live xelatex on UCR HPCC).
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
module load pandoc
TL=/opt/linux/rocky/8.x/x86_64/pkgs/texlive/20220403/bin/x86_64-linux
pandoc REPORT.md -o REPORT.pdf --pdf-engine=$TL/xelatex -H pdf_header.tex \
  -V mainfont="DejaVu Sans" -V fontsize=10pt -V geometry:margin=2cm -V colorlinks=true
