# Pipeline for tracking the MIKK panel behavioural data

## Steps for running

1. Raw data can be found:

* On the FTP site: `indigene_ftp`

* Here on EBI's Codon cluster: `/nfs/ftp/private/indigene_ftp/upload/behaviour/transfer/20191111_panel_1`

2. 20191120_1217_106-2_R_C.avi comes in two parts due to accidental cessation of the recording between the open field and novel object assay. These need to get stitched together using this script: [workflow/scripts/join_20191120_1217_106-2_R_C.py](https://github.com/brettellebi/MIKK_F0_tracking/blob/master/workflow/scripts/join_20191120_1217_106-2_R_C.py) (instructions included in the header of the script; conda environment here: `workflow/envs/opencv_4.5.5.yaml`)

3. Create Snakemake environment

```bash
conda env create --file=workflow/envs/snakemake_6.12.1.yaml
```

4. They key file containing the parameters for each video is here: [config/samples.csv](https://github.com/brettellebi/MIKK_F0_tracking/blob/master/config/samples.csv).

5. (Split) videos with unresolved tracking errors are here: [config/tracking_fails.csv](https://github.com/brettellebi/MIKK_F0_tracking/blob/master/config/tracking_fails.csv). The workflow automaticaly excludes these. To include them – for debugging – hash out the desired row.

6. Run Snakemake by running the code set out here: [workflow/init.sh](https://github.com/brettellebi/MIKK_F0_tracking/blob/master/workflow/init.sh). Note this is not intended to be run as a script, but is rather just a file that contains useful bash code. It also shows how to build the containers used in the pipeline, and how to use the custom RStudio Server container for this analysis on the cluster.
