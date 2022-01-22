# Send stdout and stderr to log file
import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

# Import libraries
import pandas as pd
import numpy as np
import cv2 as cv
import os

# Get variables

## Debugging
IN_FILE = "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/raw_videos/20191111_1527_5-1_L_A.avi"
SAMPLES_FILE = "config/samples.csv"
SAMPLE = "20191111_1527_5-1_L_A"
OUT_FILE = "results/split_coord_images/20191111_1527_5-1_L_A.jpg"

## True
IN_FILE = snakemake.input.video
SAMPLES_FILE = snakemake.params.samples_file
SAMPLE = snakemake.params.sample
OUT_FILE = snakemake.output.fig

# Read samples_file
samples_df = pd.read_csv(SAMPLES_FILE, comment="#", skip_blank_lines=True, index_col=0)

# Get date
date = int(samples_df.loc[SAMPLE, "date"])

# Get start frame for open field assay
start = int(samples_df.loc[SAMPLE, "of_start"])

# Get crop adjustment values
## note: Negative values for top/bottom shift boundary up
## note: Negative values for left/right shift boundary left
adj_top = int(samples_df.loc[SAMPLE, "adj_top"])
adj_bottom = int(samples_df.loc[SAMPLE, "adj_bottom"])
adj_left = int(samples_df.loc[SAMPLE, "adj_left"])
adj_right = int(samples_df.loc[SAMPLE, "adj_right"])

# Read video from file
cap = cv.VideoCapture(IN_FILE)
print("Video captured")
# Frame width and height
wid = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
hei = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))

# Set adjusted midpoints
mid_x = round(((wid - 1) / 2) + adj_right)
mid_y = round(((hei - 1) / 2) + adj_top)

# Capture start frame
cap.set(cv.CAP_PROP_POS_FRAMES, start)

# Read frame
ret, frame = cap.read()

# Add vertical line 
start_point = (mid_x, 0)
end_point = (mid_x, hei)
color = (255,0,0)
thickness = 1
frame = cv.line(frame, start_point, end_point, color, thickness)

# Add horizontal line
start_point = (0, mid_y)
end_point = (wid, mid_y)
color = (255,0,0)
thickness = 1
frame = cv.line(frame, start_point, end_point, color, thickness)

# Add lines for outer boundaries for 20191111 videos (which have a black outer boundary for some reason)
left_side_width = 288
right_side_width = 290
if date == 20191111:
    left_start = (left_side_width, 0)
    left_end = (left_side_width, hei)
    right_start = (wid - right_side_width, 0)
    right_end = (wid - right_side_width, hei)
    # Add vertical line for left boundary
    frame = cv.line(frame, left_start, left_end, color, thickness)
    # Add vertical line for right boundary
    frame = cv.line(frame, right_start, right_end, color, thickness)

# Write frame
cv.imwrite(OUT_FILE, frame)

cap.release()
