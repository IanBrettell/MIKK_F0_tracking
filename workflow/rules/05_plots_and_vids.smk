rule get_line_ranks_and_colours:
    input:
        rules.merge_csvs.output,
    output:
        fig = "book/figs/line_mean_speed/line_mean_speed_{interval}.png",
        csv = "config/line_colours/line_colours_{interval}.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_line_ranks_and_colours/{interval}.log"
        ),
    params:
        interval = "{interval}"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 15000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/get_line_ranks_and_colours.R"

# Add coloured labels to original videos
rule coloured_trails_labels:
    input:
        video_object=rules.track_videos.output.video_obj,
        trajectories=get_trajectories_file,
        colours = rules.get_line_ranks_and_colours.output.csv,
    output:
        os.path.join(
            config["working_dir"],
            "tracked/{interval}/{assay}/{sample}/{quadrant}.avi",
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/coloured_trails_labels/{interval}/{assay}/{sample}/{quadrant}.log"
        ),
    params:
        sample = "{sample}",
        ref_loc = get_ref_loc,
    container:
        config["idtrackerai"]
    resources:
        mem_mb=5000,
    script:
        "../scripts/coloured_trails_labels.py"

rule stitch_tracked_vids:
    input:
        expand(os.path.join(
            config["working_dir"],
            "tracked/{{interval}}/{{assay}}/{{sample}}/{quadrant}.avi"
            ),
                quadrant = QUADRANTS_ALL
        ),
    output:
        os.path.join(
            config["working_dir"],
            "stitched/{interval}/{assay}/{sample}.avi",
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/stitch_tracked_vids/{interval}/{assay}/{sample}.log"
        ),
    container:
        config["opencv"]
    resources:
        mem_mb=5000,
    script:
        "../scripts/stitch_tracked_vids.py"

rule path_frames_to_vid:
    input:
        data = rules.merge_csvs.output,
        stitched_vid = rules.stitch_tracked_vids.output,
        colours = rules.get_line_ranks_and_colours.output.csv,
        dims = rules.get_split_video_dims.output,
    output:
        os.path.join(
            config["working_dir"],
            "path_vids/{interval}/{assay}/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/path_frames_to_vid/{interval}/{assay}/{sample}.log"
        ),
    params:
        assay = "{assay}",
        sample = "{sample}",
        tmpdir = config["tmpdir"]
    resources:
        mem_mb = 20000,
    container:
        config["opencv"]
    script:
        "../scripts/path_frames_to_vid.py"

rule compile_two_panel_vid:
    input:
        labels = rules.stitch_tracked_vids.output[0],
        paths = rules.path_frames_to_vid.output[0],
    output:
        os.path.join(
            config["working_dir"],
            "two_panel_vids/{interval}/{assay}/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/compile_two_panel_vid/{interval}/{assay}/{sample}.log"
        ),
    resources:
        mem_mb = 5000,
    container:
        config["opencv"]
    script:
        "../scripts/compile_two_panel_vid.py"

rule two_panel_short:
    input:
        rules.compile_two_panel_vid.output,
    output:
        avi = os.path.join(
            config["working_dir"],
            "two_panel_short_vids/{interval}/{assay}/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/two_panel_short/{interval}/{assay}/{sample}.log"
        ),
    params:
        tot_sec = 60,
    resources:
        mem_mb = 5000,
    container:
        config["opencv"]
    script:
        "../scripts/two_panel_short.py"         

# Creates point/path plots for each ref/test fish
# With each point coloured by HMM state
rule hmm_path_frames_to_vid:
    input:
        hmm = config["hmm_out"],
        dims = rules.get_split_video_dims.output,
    output:
        os.path.join(
            config["working_dir"],
            "hmm_path_vids/{interval}/dist_angle/15/{assay}/{sample}_{ref_test}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm_path_frames_to_vid/{interval}/dist_angle/15/{assay}/{sample}_{ref_test}.log"
        ),
    params:
        assay = "{assay}",
        sample = "{sample}",
        ref_test = "{ref_test}",
        fps = get_fps,
        tmpdir = os.path.join(
            config["working_dir"],
            "tmp_frames/0.05/dist_angle/15"
        ),
    resources:
        mem_mb = 20000,
    container:
        config["opencv"]
    script:
        "../scripts/hmm_path_frames_to_vid.py"

# Create four-panel videos with labelled videos (TL), path videos (TR),
# and HMM path videos for test (BL) and ref (BR) fishes.
rule compile_four_panel_vid:
    input:
        labels = rules.stitch_tracked_vids.output[0],
        paths = rules.path_frames_to_vid.output[0],
        hmm_ref = os.path.join(
            config["working_dir"],
            "hmm_path_vids/{interval}/dist_angle/15/{assay}/{sample}_ref.avi"
        ),
        hmm_test = os.path.join(
            config["working_dir"],
            "hmm_path_vids/{interval}/dist_angle/15/{assay}/{sample}_test.avi"
        ),
    output:
        os.path.join(
            config["working_dir"],
            "four_panel_vids/{interval}/dist_angle/15/{assay}/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/compile_four_panel_vid/{interval}/dist_angle/15/{assay}/{sample}.log"
        ),
    params:
        fps = get_fps,
    resources:
        mem_mb = 5000,
    container:
        config["opencv"]
    script:
        "../scripts/compile_four_panel_vid.py"

# Shorten to 
rule four_panel_short:
    input:
        rules.compile_four_panel_vid.output,
    output:
        avi = os.path.join(
            config["working_dir"],
            "four_panel_short_vids/{interval}/dist_angle/15/{assay}/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/four_panel_short/{interval}/dist_angle/15/{assay}/{sample}.log"
        ),
    params:
        tot_sec = 60,
    resources:
        mem_mb = 5000,
    container:
        config["opencv"]
    script:
        "../scripts/four_panel_short.py"    