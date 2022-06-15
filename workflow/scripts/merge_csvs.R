# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
IN_FILES = list("/hps/nobackup/birney/users/ian/MIKK_F0_tracking/with_metrics/novel_object/20191121_1454_iCab_L_C/q1/0.08.csv",
                "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/with_metrics/novel_object/20191120_1026_84-2_R_C/q2/0.08.csv")

## True
IN_FILES = snakemake@input
OUT_FILE = snakemake@output[[1]]

# Get metadata from IN_FILES
names(IN_FILES) = purrr::map(IN_FILES, function(IN_FILE){
  elements = IN_FILE %>% 
    stringr::str_split(pattern = "/", simplify = T)
  n_el = length(elements)
  
  out = paste(elements[n_el - 3], elements[n_el-2], elements[n_el - 1], sep = "_")
  
  return(out)
}) 

# Read in file and process

out = purrr::map(IN_FILES, function(IN_FILE){
  readr::read_csv(IN_FILE)
}) %>% 
  # bind into single DF
  dplyr::bind_rows(.id = "video") %>% 
  # separate metadata
  tidyr::separate(col = video, into = c("assay1", "assay2", "date", "time", "test_fish", "tank_side", "block", "quadrant"),
                  sep = "_") %>% 
  # unite assay columns
  tidyr::unite(col = "assay",
               assay1, assay2,
               sep = "_",
               remove = T) %>% 
  # create "ref_fish", `pat_line` and `mat_line` columns,
  dplyr::mutate(ref_fish = "iCab",
                pat_line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                            fish == "test" ~ test_fish),
                mat_line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                            fish == "test" ~ test_fish)) %>%
  # remove unnecessary columns
  dplyr::select(-c(x_lag1, y_lag1, x_lag2, y_lag2, distance_b, block)) %>% 
  # order columns
  dplyr::select(assay, date, time, ref_fish, test_fish, dplyr::everything()) %>% 
  # drop NAs (they don't work with the HMM)
  tidyr::drop_na() 

# Write to file

readr::write_csv(out, OUT_FILE)
