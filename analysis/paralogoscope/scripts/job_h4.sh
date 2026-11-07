#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 16G --time 00:30:00 -p short -J pf_h4
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/h4.%j.log
set -e
E=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore
/usr/bin/python3.12 $E/scripts/make_h4_report.py
module load texlive/20220403
module load pandoc 2>/dev/null || true
cd $E/reports
pandoc H4_duplication_mode.md -o H4_duplication_mode.pdf --toc --pdf-engine=/opt/linux/rocky/8.x/x86_64/pkgs/texlive/20220403/bin/x86_64-linux/xelatex -V mainfont="DejaVu Sans" -V monofont="DejaVu Sans Mono" -V geometry:margin=2cm -V fontsize=9pt -V colorlinks=true --resource-path=$E/reports:$E
ls -la H4*
