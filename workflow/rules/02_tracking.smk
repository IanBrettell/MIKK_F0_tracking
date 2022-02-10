# pull tracking parameters from config/samples.csv
def get_vid_length(wildcards):
    if wildcards.assay == "open_field":
        start = samples_df.loc[samples_df["sample"] == wildcards.sample, "of_start"]
        end = samples_df.loc[samples_df["sample"] == wildcards.sample, "of_end"]
    elif wildcards.assay == "novel_object":
        start = samples_df.loc[samples_df["sample"] == wildcards.sample, "no_start"]
        end = samples_df.loc[samples_df["sample"] == wildcards.sample, "no_end"]
    vid_length = int(end) - int(start)
    return(vid_length)

def get_bgsub(wildcards):
    if wildcards.assay == "open_field":
        bgsub = int(samples_df.loc[samples_df["sample"] == wildcards.sample, "bgsub_of"])
    elif wildcards.assay == "novel_object":
        bgsub = int(samples_df.loc[samples_df["sample"] == wildcards.sample, "bgsub_no"])
    return(bgsub)

def get_intensity_floor(wildcards):
    if wildcards.assay == "open_field":
        int_floor = int(samples_df.loc[samples_df["sample"] == wildcards.sample, "intensity_floor_of"])
    elif wildcards.assay == "novel_object":
        int_floor = int(samples_df.loc[samples_df["sample"] == wildcards.sample, "intensity_floor_no"])
    return(int_floor)

def get_intensity_ceiling(wildcards):
    if wildcards.assay == "open_field":
        int_ceiling = int(samples_df.loc[samples_df["sample"] == wildcards.sample, "intensity_ceiling_of"])
    elif wildcards.assay == "novel_object":
        int_ceiling = int(samples_df.loc[samples_df["sample"] == wildcards.sample, "intensity_ceiling_no"])
    return(int_ceiling)

def get_parameters(wildcards):
    if wildcards.assay == "open_field":
        start = samples_df.loc[samples_df["sample"] == wildcards.sample, "of_start"]
        end = samples_df.loc[samples_df["sample"] == wildcards.sample, "of_end"]

# adapt memory usage for tracking videos
def get_mem_mb(wildcards, attempt):
    return attempt * 10000

# Track videos with idtrackerai
## Note: `output` is set as `trajectories.npy` instead of `trajectories_wo_gaps.npy`, presumably because
## in videos where there are no crossovers, the latter file is not produced.
rule track_videos:
    input:
        rules.split_videos.output
    output:
        os.path.join(config["data_store_dir"], "split/{assay}/session_{sample}_{quadrant}/trajectories/trajectories.npy"),
    log:
        os.path.join(config["working_dir"], "logs/track_videos/{assay}/{sample}/{quadrant}.log"),
    params:
        vid_length = get_vid_length,
        vid_name = "{sample}_{quadrant}",
        intensity_floor = get_intensity_floor,
        intensity_ceiling = get_intensity_ceiling,
    resources:
        mem_mb = get_mem_mb
    container:
        config["idtrackerai"]
    shell:
        """
        idtrackerai terminal_mode \
            --_video {input} \
            --_bgsub 'True' \
            --_range [0,{params.vid_length}] \
            --_intensity [{params.intensity_floor},{params.intensity_ceiling}] \
            --_session {params.vid_name} \
            --exec track_video \
                2> {log}
        """

# Convert numpy arrays to .csv files
rule trajectories_to_csv:
    input:
        trajectories = rules.track_videos.output,
        script = "workflow/scripts/trajectories_to_csv.py"
    output:
        os.path.join(config["data_store_dir"], "split/{assay}/session_{sample}_{quadrant}/trajectories/trajectories.trajectories.csv")
    log:
        os.path.join(config["working_dir"], "logs/trajectories_to_csv/{assay}/{sample}/{quadrant}.log"),
    params:
        in_path = os.path.join(config["data_store_dir"], "split/{assay}/session_{sample}_{quadrant}")
    shell:
        """
        python {input.script} {params.in_path}
        """
