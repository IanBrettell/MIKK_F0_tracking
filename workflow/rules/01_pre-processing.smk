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

#rule draw_split_lines:
#    input:
#        video = join(config["raw_data_dir"], "{sample}.avi"),
#        samples_file = config["samples_file"]
#
#rule split_videos:
#    input:
#        video = join(config["raw_data_dir"], "{sample}.avi"),
#        samples_file = config["samples_file"]
#    output:
#        "split/{assay}/{sample}_{quadrant}_{assay}.mp4"
#    params:
#        samples_file = config["samples_file"]
#    container:
#        config["opencv"]
#    script:
#        config["split_script"]