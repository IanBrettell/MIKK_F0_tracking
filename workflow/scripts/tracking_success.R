# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
IN_FILES = list("/nfs/research/birney/users/ian/MIKK_F2_tracking/final_tracks/novel_object/20191115_1338_134-1_R_A/q1.csv",
                "/nfs/research/birney/users/ian/MIKK_F2_tracking/final_tracks/novel_object/20191115_1338_134-1_R_A/q2.csv")

## True
IN_FILES = snakemake@input
OUT_FILE = snakemake@output[[1]]

# Read files and process

final_df = purrr::map(IN_FILES, function(IN_FILE){
    # read in file
    df = readr::read_csv(IN_FILE, col_types = c("iddddd"))
    # create tibble
    out = tibble::tibble(
        "path" = IN_FILE,
        "total_seconds" = max(df$seconds),
        "success_count" = df %>% 
            dplyr::filter(complete.cases(.)) %>%
            nrow(.),
        "frame_count" = nrow(df)
        ) %>%
        # get `sample` and `assay`
        tidyr::separate(path,
                        into = c(rep(NA, 8), "assay", "sample", "quadrant"),
                        sep = "/") %>%
        # remove .csv extension from `quadrant`
        dplyr::mutate(quadrant = quadrant %>%
            stringr::str_remove(".csv")) %>%
        # get proportion of success
        dplyr::mutate("prop_success" = success_count / frame_count )

    return(out)

}) %>% 
    # bind into single DF
    dplyr::bind_rows() %>%
    # reorder columns
    dplyr::select(sample, assay, quadrant, everything()) %>%
    # order rows by `prop_success`
    dplyr::arrange(prop_success)

# Write to file

readr::write_csv(final_df, OUT_FILE)
