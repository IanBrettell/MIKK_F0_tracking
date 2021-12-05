#!/bin/bash

# NOTE: raw data locations:
# FTP: /nfs/ftp/private/indigene_ftp/upload/behaviour/transfer/20191111_panel_1/all_to_analyse

####################
# Codon cluster
####################

ssh codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -q datamover -M 20000 -Is bash
cd /hps/software/users/birney/ian/repos/MIKK_F0_tracking
conda activate snakemake_6.10.0
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -q {cluster.queue} -n {cluster.n} -M {cluster.memory} -o {cluster.outfile}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  -s workflow/Snakefile \
  -p

####################
# EBI cluster
####################

ssh codon
bsub -M 20000 -XF -Is bash
/nfs/software/birney/Fiji.app/ImageJ-linux64

####################
# Build custom containers
####################

RCONT=/hps/software/users/birney/ian/containers/MIKK_F0_tracking/R

OPENCVCONT=/hps/software/users/birney/ian/containers/opencv_4.5.1.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $OPENCVCONT \
    workflow/envs/opencv_4.5.1.def