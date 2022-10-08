# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = "/hps/software/users/birney/ian/repos/MIKK_F0_tracking/book/figs/four_panel_frame_grabs/0.05/dist_angle/15/novel_object/novel_object_20191112_1236_18-2_L_A_300.png"

## True
IN = snakemake@input[[1]]
OUT = snakemake@output[["pdf"]]


# Get original dimensions
img = magick::image_read(IN)
dims = dim(magick::image_data(img))

WID = dims[2]
HEI = dims[3]

# Read in figure
fig = cowplot::ggdraw() + 
    cowplot::draw_image(IN)

# Save as pdf
ggsave(OUT, 
       fig,
       device = "pdf",
       width = WID, 
       hei = HEI, 
       units = "px")
