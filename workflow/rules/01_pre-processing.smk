# Copy videos from FTP to Codon (for ease downstream)
rule copy_videos:
    input:
        os.path.join(
            config["raw_data_dir"], 
            "{sample}.avi"
        ),
    output:
        os.path.join(
            config["working_dir"], 
            "raw_videos/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"], 
            "logs/copy_videos/{sample}.log"
        ),
    shell:
        """
        cp {input} {output} \
            2> {log}
        """

# Create copies recoded with opencv to create smaller videos
# That can be viewed with idtrackerai
rule recode_videos:
    input:
        rules.copy_videos.output,
    output:
        os.path.join(
            config["working_dir"],
            "recoded/{sample}.avi"
        ),        
    log:
        os.path.join(
            config["working_dir"],
            "logs/recode_videos/{sample}.log"
        ),
    container:
        config["opencv"]
    resources:
        mem_mb = 3000
    script:
        "../scripts/recode_videos.py"

# Generate single-frame grab showing coordinate of splits
rule set_split_coords:
    input:
        video = rules.recode_videos.output,
    output:
        fig = "results/split_coord_images/{assay}/{sample}.png",
    log:
        os.path.join(
            config["working_dir"], 
            "logs/set_split_coords/{assay}/{sample}.log"
        )
    params:
        assay = "{assay}",
        sample = "{sample}",
        samples_file = lambda wildcards: config["samples_long"]
    container:
        config["opencv"]
    resources:
        mem_mb = 500
    script:
        "../scripts/set_split_coords.py"


# Split videos into quadrants and assays (1 raw video * 4 quadrants * 2 assays = 8 output videos)
rule split_videos:
    input:
        rules.copy_videos.output,
    output:
        os.path.join(config["data_store_dir"], "split/{assay}/{sample}_{quadrant}.mp4"),
    log:
        os.path.join(config["working_dir"], "logs/split_videos/{assay}/{sample}/{quadrant}.log"),
    params:
        sample = "{sample}",
        assay = "{assay}",
        quadrant = "{quadrant}",
        samples_file = lambda wildcards: config["samples_file"]
    container:
        config["opencv"]
    script:
        "../scripts/split_videos.py"


