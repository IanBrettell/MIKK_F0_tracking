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

# Import libraries

import numpy as np
import pandas as pd
import cv2
import os

# Get variables

## Debug
PATHS = "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/path_vids/0.08/open_field/20191113_1557_22-1_L_A.avi"
SECOND = 110
OUT = "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/tmp_out.png"

# True
PATHS = snakemake.input.paths
SECOND = int(snakemake.params.target_second)
OUT = snakemake.output[0]

########################
## Capture videos
########################

cap_path = cv2.VideoCapture(PATHS)

####################
# Get N frames of video, width, and height
####################

IN = {
    "path": {"path": PATHS},
    }

for key in IN:
    cap = cv2.VideoCapture(IN[key]['path'])
    IN[key]['n_frames'] = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    IN[key]['fps'] = int(cap.get(cv2.CAP_PROP_FPS))
    IN[key]['wid'] = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    IN[key]['hei'] = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

# Get FPS

FPS = IN['path']['fps']

# Get target frame

FRAME = FPS*SECOND

#######################
# Write frame
#######################

# Set frames
cap_path.set(cv2.CAP_PROP_POS_FRAMES, FRAME)

# Capture frames
ret, out = cap_path.read()

# Write to file
cv2.imwrite(OUT, out)
