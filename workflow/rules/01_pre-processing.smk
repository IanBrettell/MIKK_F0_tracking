rule copy_videos:
    input:
        os.path.join(config["raw_data_dir"], "{sample}.avi"),
    output:
        os.path.join(config["working_dir"], "raw_videos/{sample}.avi"),
    log:
        os.path.join(config["working_dir"], "logs/copy_videos/{sample}.log"),
    shell:
        """
        cp {input} {output}
        """

rule set_split_coords:
    input:
        video = os.path.join(config["working_dir"], "raw_videos/{sample}.avi"),
        samples_file = config["samples_file"]
    output:
        out = os.path.join(config["working_dir"], "results/split_coord_images/{sample}.jpg")
    log:
        os.path.join(config["working_dir"], "logs/set_split_coords/{sample}.log"),
    params:
        sample = "{sample}"
    container:
        config["opencv"]
    script:
        "/hps/software/users/birney/ian/repos/MIKK_F0_tracking/workflow/scripts/set_split_coords.py" 

rule split_videos:
    input:
        video = os.path.join(config["raw_data_dir"], "{sample}.avi"),
        samples_file = config["samples_file"]
    output:
        os.path.join(config["data_store_dir"], "split/{assay}/{sample}_{quadrant}_{assay}.mp4")
    params:
        assay = "{assay}"
    container:
        config["opencv"]
    script:
        "../scripts/split_videos.py"