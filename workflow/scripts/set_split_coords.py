#!/usr/bin/env python3

# Send stdout and stderr to log file

import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

print("Script commenced")

## Import libraries
#import pandas as pd
#import numpy as np
#import cv2 as cv
#import os
#
#print("Libraries loaded")
#
## Get variables
#
### Debugging
##IN_FILE = "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/raw_videos/20191111_1527_5-1_L_A.avi"
##SAMPLES_FILE = "config/samples.csv"
##SAMPLE = "20191111_1527_5-1_L_A"
##OUT_FILE = "results/split_coord_images/20191111_1527_5-1_L_A.jpg"
#
### True
#IN_FILE = snakemake.input.video[0]
#SAMPLES_FILE = snakemake.input.samples_file[0]
#SAMPLE = snakemake.params.sample[0]
#OUT_FILE = snakemake.output[0]
#
## Read samples_file
#
#samples_df = pd.read_csv(SAMPLES_FILE, comment="#", skip_blank_lines=True, index_col=0)
#print("Samples df read")
### Get start frame for open field assay
#
#start = int(samples_df.loc[SAMPLE, "of_start"])
#
## Get crop adjustment values
### NOTE: Negative values for top/bottom shift boundary up
### NOTE: Negative values for left/right shift boundary left
#adj_top = int(samples_df.loc[SAMPLE, "adj_top"])
#adj_bottom = int(samples_df.loc[SAMPLE, "adj_bottom"])
#adj_left = int(samples_df.loc[SAMPLE, "adj_left"])
#adj_right = int(samples_df.loc[SAMPLE, "adj_right"])
#
## Read video from file
#cap = cv.VideoCapture(IN_FILE)
#print("Video captured")
## Frame width and height
#wid = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
#hei = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))
#
## Set adjusted midpoints
#mid_x = round(((wid - 1) / 2) + adj_right)
#mid_y = round(((hei - 1) / 2) + adj_top)
#
## Capture start frame
#cap.set(cv.CAP_PROP_POS_FRAMES, start)
#
## Read frame
#ret, frame = cap.read()
#
## Add vertical line 
#start_point = (mid_x, 0)
#end_point = (mid_x, hei)
#color = (255,0,0)
#thickness = 1
#frame = cv.line(frame, start_point, end_point, color, thickness)
#
## Add horizontal line
#start_point = (0, mid_y)
#end_point = (wid, mid_y)
#color = (255,0,0)
#thickness = 1
#frame = cv.line(frame, start_point, end_point, color, thickness)
#
## Write frame
#cv.imwrite(OUT_FILE, frame)
#
#cap.release()
#OUT_FILE.release()
