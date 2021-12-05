rule copy_videos:
    input:


rule split_videos:
    input:
        join(config["input_dir"], "{sample}.avi")
    output:
        "split/{assay}/{sample}_{quadrant}_{assay}.mp4"
    params:
        samples_file = config["samples_file"]
    container:
        config["opencv"]
    script:
        config["split_script"]