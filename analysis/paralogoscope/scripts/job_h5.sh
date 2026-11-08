#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 16G --time 00:30:00 -p short -J pf_h5
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore/data/h5.%j.log
set -e
E=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot/eval/explore
/usr/bin/python3.12 $E/scripts/make_h5_report.py
module load texlive/20220403
module load pandoc 2>/dev/null || true
cd $E/reports
pandoc H5_block_structure_window_codeml.md -o H5_block_structure_window_codeml.pdf --toc --pdf-engine=/opt/linux/rocky/8.x/x86_64/pkgs/texlive/20220403/bin/x86_64-linux/xelatex -V mainfont="DejaVu Sans" -V monofont="DejaVu Sans Mono" -V geometry:margin=2cm -V fontsize=9pt -V colorlinks=true --resource-path=$E/reports:$E
ls -la H5*
