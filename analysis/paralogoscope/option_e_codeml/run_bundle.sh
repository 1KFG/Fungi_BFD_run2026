#!/usr/bin/bash -l
# Option E: codeml second pass for genomes with a possible interior dS peak. Same settings as the yn00 run except --dnds-method codeml.
# Usage (sbatch): sbatch -c 16 --mem 32G --time 04:00:00 -p batch run_bundle.sh option_e_codeml/bundle_001.txt
# Does not touch any live Nextflow run: reads nf_full/input/cds and paralogoscope_full/wgd_dmd, writes only under option_e_codeml/.
#SBATCH -J optE_codeml
module load apptainer squashfuse fuse
P=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/paralogoscope_pilot
O=$P/option_e_codeml; mkdir -p $O/ks $O/logs $O/failed
BUNDLE=${1:?bundle file}; export TMPDIR=${SCRATCH:?}
SIF=/bigdata/stajichlab/shared/singularity_cache/wgd2_complete-2.0.38-subopt.sif
SCRIPT=$P/Fungi_BFD_wt/nextflow/bin/wgd_pair_ks.py
CPUS=${SLURM_CPUS_PER_TASK:-8}
BIND="--bind /bigdata:/bigdata,/rhome:/rhome,$TMPDIR:$TMPDIR"
while read -r id; do
  [ -z "$id" ] && continue
  out=$O/ks/${id}.cds-transcripts.fa.tsv.ks.tsv.gz
  [ -s "$out" ] && { echo "$id: done, skipped"; continue; }
  tab=$(awk -F'\t' -v g="$id" 'NR>1 && $1==g {print $5}' $O/genomes.tsv)
  tab=${tab:-1}
  case $tab in 1) table=1; icode=0;; 12) table=12; icode=8;; *) echo "$id: genetic code $tab not validated" >> $O/failed/$(basename $BUNDLE).txt; continue;; esac
  cdsgz=$P/nf_full/input/cds/${id}.cds-transcripts.fa.gz; cds=$TMPDIR/${id}.cds-transcripts.fa
  zcat $cdsgz > $cds
  fam=$P/paralogoscope_full/wgd_dmd/${id}.cds-transcripts.fa.tsv
  tmp=$TMPDIR/${id}.ks.tsv.gz
  echo "$(date +%T) $id start (table $table)"
  if apptainer exec $BIND --env PYTHONNOUSERSITE=1 $SIF python3 $SCRIPT --cds $cds --families $fam --out $tmp --threads $CPUS \
      --min-identity 0.4 --min-ng-codons 30 --dnds-method codeml --node-sample 5 --node-max-draws 20 --gate-max-ds 2.5 --gate-min-codons 30 \
      --table $table --icode $icode --tmpdir $TMPDIR < /dev/null > $O/logs/${id}.log 2>&1 && [ -s $tmp ]; then
     for f in $tmp $tmp.nodes.tsv.gz; do [ -s $f ] && cp $f $O/ks/$(basename $f | sed "s/^${id}.ks.tsv.gz/${id}.cds-transcripts.fa.tsv.ks.tsv.gz/"); done
     echo "$(date +%T) $id ok"
  else echo "$id: failed (see logs/${id}.log)" >> $O/failed/$(basename $BUNDLE).txt; echo "$id FAILED"; fi
  rm -f $cds $tmp $tmp.nodes.tsv.gz
done < $BUNDLE
