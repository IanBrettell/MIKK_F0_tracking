#!/bin/bash

# NOTE: raw data locations:
# FTP: /nfs/ftp/private/indigene_ftp/upload/behaviour/transfer/20191111_panel_1/all_to_analyse

####################
# Codon cluster
####################

ssh codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -Is bash
cd /hps/software/users/birney/ian/repos/MIKK_F0_tracking
conda activate snakemake_6.12.1
#conda activate snakemake_6.13.1
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

# R
RCONT=/hps/nobackup/birney/users/ian/containers/MIKK_F0_tracking/R_4.1.2.sif
singularity build --remote \
    $RCONT \
    workflow/envs/R_4.1.2.def

# Open CV (python)
OPENCVCONT=/hps/nobackup/birney/users/ian/containers/MIKK_F0_tracking/opencv_4.5.1.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $OPENCVCONT \
    workflow/envs/opencv_4.5.1.def

####################
# Run RStudio Server
####################

ssh proxy-codon
bsub -q datamover -M 50000 -Is bash
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
RCONT=/hps/nobackup/birney/users/ian/containers/MIKK_F0_tracking/R_4.1.2.sif
singularity shell --bind /hps/software/users/birney/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/software/users/birney/ian/tmp:/tmp \
                  --bind /hps/software/users/birney/ian/run:/run \
                  $RCONT

rserver \
    --rsession-config-file /hps/software/users/birney/ian/repos/MIKK_F0_tracking/workflow/envs/rsession.conf \
    --server-user brettell

ssh -L 8787:hl-codon-37-04:8787 proxy-codon