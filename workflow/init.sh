#!/bin/bash

# NOTE: raw data locations:
# FTP: /nfs/ftp/private/indigene_ftp/upload/behaviour/transfer/20191111_panel_1/all_to_analyse

####################
# Codon cluster
####################

ssh codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -Is bash
# If needing to copy videos from FTP (rule copy_videos),
# Need to use the datamover queue so that it can see the FTP drive:
# bsub -M 20000 -q datamover -Is bash
cd /hps/software/users/birney/ian/repos/MIKK_F0_tracking
conda activate snakemake_6.12.1
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -q {cluster.queue} -n {cluster.n} -M {cluster.memory} -o {cluster.outfile}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  --restart-times 0 \
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
    workflow/envs/R_4.1.2/R_4.1.2.def

# Open CV (python)
OPENCVCONT=/hps/nobackup/birney/users/ian/containers/MIKK_F0_tracking/opencv_4.5.1.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $OPENCVCONT \
    workflow/envs/opencv_4.5.1.def

# idtrackerai
IDCONT=/hps/nobackup/birney/users/ian/containers/MIKK_F0_tracking/idtrackerai.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $IDCONT \
    workflow/envs/idtrackerai.def

####################
# Run RStudio Server
####################

ssh proxy-codon
bsub -M 50000 -Is bash
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
CONT=/hps/nobackup/birney/users/ian/containers/MIKK_F0_tracking/R_4.1.2.sif
singularity shell --bind /hps/nobackup/birney/users/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/nobackup/birney/users/ian/tmp:/tmp \
                  --bind /hps/nobackup/birney/users/ian/run:/run \
                  $CONT
rstudio-server kill-all
rserver \
    --rsession-config-file /hps/software/users/birney/ian/repos/MIKK_F0_tracking/workflow/envs/R_4.1.2/rsession.conf \
    --server-user brettell

ssh -L 8787:hl-codon-37-04:8787 proxy-codon

####################
# Run idtrackerai
####################

## GUI
ssh codon
bsub -M 10000 -q gui -XF -Is bash
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity shell docker://saulpierottiebi/idtrackerai_cpu_gui:latest
idtrackerai

## Terminal mode
singularity shell docker://saulpierottiebi/idtrackerai_cpu_gui:latest
INPUT_VIDEO=/nfs/research/birney/users/ian/MIKK_F0_tracking/split/open_field/20191121_1454_iCab_L_C_q4.mp4
VID_LENGTH=18178
idtrackerai terminal_mode \
            --_video $INPUT_VIDEO \
            --_bgsub 'True' \
            --_range [0,$VID_LENGTH] \
            --_session 20191121_1454_iCab_L_C_q4 \
            --exec track_video

# convert to .csv
python /hps/software/users/birney/ian/repos/MIKK_F0_tracking/workflow/scripts/trajectories_to_csv.py /nfs/research/birney/users/ian/MIKK_F0_tracking/split/open_field/session_20191121_1454_iCab_L_C_q4

idtrackerai terminal_mode \
            --_video $INPUT_VIDEO \
            --_bgsub 'True' \
            --_intensity [$int_floor,$int_ceiling] \
            --_area [$area_floor,$area_ceiling] \
            --_range [0,$vid_length] \
            --_nblobs 2 \
            --_session $in_sample \
            --exec track_video 


####################
# Copy videos from cluster to local
####################

# To set tracking parameters
rsync -aP brettell@codon:/hps/nobackup/birney/users/ian/MIKK_F0_tracking/recoded ~/Downloads/MIKK_F0
