# Get variables

IN_VIDEO = "/nfs/ftp/private/indigene_ftp/upload/behaviour/transfer/20191111_panel_1/all_to_analyse/20191120_1233_94-1_L_C.avi"
SAMPLES_FILE = "config/samples.csv"
TARGET_SAMPLE = "20191120_1233_94-1_L_C"

IN_VIDEO = snakemake@input[["video"]]
SAMPLES_FILE = snakemake@input[["samples_file"]]
OUT_FILE = snakemake@output[[1]]
LOG_FILE = snakemake@log[[1]]
TARGET_SAMPLE = snakemake@params[["target_sample"]]

# Send output to log

log <- file(LOG_FILE, open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(magick)

# Read in samples file

samples_df = readr::read_csv(SAMPLES_FILE)
sample_name = snakemake@params[["sample_name"]]

# Get lane coords

## Total lanes in video
lanes_n = samples_df %>%
  dplyr::filter(SAMPLE == sample_name) %>%
  dplyr::pull(TOTAL_LANES)

## Get lane coords
lane_coords = samples_df %>%
  dplyr::filter(SAMPLE == sample_name) %>%
  dplyr::select(starts_with("END_LANE")) %>%
  subset(select = 1:lanes_n -1) %>%
  unlist(use.names = F)

# read first frame of tiff
frame_1 = ijtiff::read_tif(snakemake@input[[1]], frames = 1)

# write to file (magick can't read directly...)
tmp_file_name = paste(sample_name, ".tmp.tif", sep = "")
ijtiff::write_tif(frame_1, tmp_file_name)

# read back in as magick image
frame_1_m = magick::image_read(tmp_file_name)

# add horizontal lines

img = magick::image_draw(frame_1_m) # make image object
abline(h = lane_coords, col = "white") # draw lines
lined_img = magick::image_scale(img, 1000) %>% # shrink
  magick::image_modulate(brightness = 700) # increase brightness
dev.off()

# write image

magick::image_write(lined_img, path = snakemake@output[[1]], format = "png")

# clean up

file.remove(tmp_file_name)