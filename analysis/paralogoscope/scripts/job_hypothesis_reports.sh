#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 16G --time 01:00:00 -p short -J pf_hypreports
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/hypreports.%j.log
set -e
E=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore
cd $E
/usr/bin/python3.12 scripts/make_hypothesis_reports.py
module load texlive/20220403
module load pandoc 2>/dev/null || true
cd $E/reports
for r in H1_interior_peaks H2_young_paralogs H3_OrsD_case_study; do
  pandoc $r.md -o $r.pdf --toc --pdf-engine=/opt/linux/rocky/8.x/x86_64/pkgs/texlive/20220403/bin/x86_64-linux/xelatex -V mainfont="DejaVu Sans" -V monofont="DejaVu Sans Mono" -V geometry:margin=2cm -V fontsize=9pt -V colorlinks=true --resource-path=$E/reports:$E/reports/figures:$E
done
ls -la
