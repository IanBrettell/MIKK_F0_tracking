# pull tracking parameters from config/samples.csv
def get_vid_length(wildcards):
    if wildcards.assay == "open_field":
        start = samples_df.loc[samples_df["sample"] == wildcards.sample, "of_start"]
        end = samples_df.loc[samples_df["sample"] == wildcards.sample, "of_end"]
    elif wildcards.assay == "novel_object":
        start = samples_df.loc[samples_df["sample"] == wildcards.sample, "no_start"]
        end = samples_df.loc[samples_df["sample"] == wildcards.sample, "no_end"]
    vid_length = end[0] - start[0]
    return(vid_length)

rule track_videos:
    input:
        rules.split_videos.output
    output:
        os.path.join(config["data_store_dir"], "split/{assay}/session_{sample}_{quadrant}/trajectories/trajectories.npy")
    log:
        os.path.join(config["working_dir"], "logs/track_videos/{assay}/{sample}/{quadrant}.log"),
    params:
        vid_length = get_vid_length,
        vid_name = "{sample}_{quadrant}"
    container:
        config["idtrackerai"]
    shell:
        """
        idtrackerai terminal_mode \
            --_video {input} \
            --_bgsub 'True' \
            --_range [0,{params.vid_length}] \
            --_session {params.vid_name} \
            --exec track_video \
                2> {log}
        """

rule trajectories_to_csv:
    input:
        trajectories = rules.track_videos.output,
        script = "workflow/scripts/trajectories_to_csv.py"
    output:
        os.path.join(config["data_store_dir"], "split/{assay}/session_{sample}_{quadrant}/trajectories/trajectories.trajectories.csv")
    log:
        os.path.join(config["working_dir"], "logs/track_videos/{assay}/{sample}/{quadrant}.log"),
    params:
        in_path = lambda wildcards, output: os.path.split(output[0])[0]
    shell:
        """
        {input.script} {params.in_path}
        """
