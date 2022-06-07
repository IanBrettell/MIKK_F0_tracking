# Choose input trajectories .csv file
def get_final_csvs(wildcards):
    # Get path of csv files
    ## Trajectories without gaps
    traj_wo_gaps_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        wildcards.sample,
        "session_" + wildcards.sample + "_" + wildcards.quadrant,
        "trajectories_wo_gaps",
        "trajectories_wo_gaps.trajectories.csv"
        )
    ## Trajectories (with gaps)
    traj_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        wildcards.sample,
        "session_" + wildcards.sample + "_" + wildcards.quadrant,
        "trajectories",
        "trajectories.trajectories.csv"
        )
    # If there is no "without gaps" file, return the "trajectories" file
    if os.path.exists(traj_wo_gaps_file):
        return(traj_wo_gaps_file)
    elif os.path.exists(traj_file):
        return(traj_file)


# Get frames-per-second
def get_fps(wildcards):
    fps = samples_df.loc[(samples_df['sample'] == wildcards.sample) & \
                         (samples_df['assay'] == wildcards.assay) & \
                         (samples_df['quadrant'] == wildcards.quadrant), \
                         'fps'].values[0]
    return(fps)

# Get relative location of reference iCab fish
def get_ref_loc(wildcards):
    cab_loc = samples_df.loc[(samples_df['sample'] == wildcards.sample) & \
                             (samples_df['assay'] == wildcards.assay) & \
                             (samples_df['quadrant'] == wildcards.quadrant), \
                             'cab_coords'].values[0]
    # `ref_loc` is nan if there are two iCabs, so convert to string
    if pd.isna(cab_loc):
        cab_loc = "NA"
    return(cab_loc)

# Assign reference and test fish IDs, and filter for frames up to 10 minutes
rule assign_ref_test:
    input:
        get_final_csvs,
    output:
        os.path.join(
            config["working_dir"],
            "final_tracks/{assay}/{sample}/{quadrant}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/assign_ref_test/{assay}/{sample}/{quadrant}.log"
        ),
    params:
        fps = get_fps,
        ref_loc = get_ref_loc,
    resources:
        mem_mb = 200
    script:
        "../scripts/assign_ref_test.py"

# Create .csv with proportion of frames tracked
rule tracking_success:
    input:
        expand(rules.assign_ref_test.output,
            zip,
            assay = ASSAYS_ZIP,
            sample = SAMPLES_ZIP,
            quadrant = QUADRANTS_ZIP                     
        ),
    output:
        "config/tracking_success.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/tracking_success/tracking_success.log"
        ),
    container:
        config["rocker_tidyverse"]
    resources:
        mem_mb = 1000
    script:
        "../scripts/tracking_success.R"
 
