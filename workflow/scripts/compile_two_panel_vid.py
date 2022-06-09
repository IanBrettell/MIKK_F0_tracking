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
LABELS = "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/stitched/novel_object/20190612_1636_icab_icab_R.avi"
PATHS = "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/path_vids/0.08/novel_object/20190612_1636_icab_icab_R.avi"
OUT = "/hps/nobackup/birney/users/ian/pilot/tmp_out.avi"

# True
LABELS = snakemake.input.labels
PATHS = snakemake.input.paths
OUT = snakemake.output[0]

########################
## Read in dims
########################
#
#dims = pd.read_csv(DIMS)
#dims = dims.loc[(dims['sample'] == SAMPLE) & (dims['assay'] == ASSAY)]
#N_FRAMES = max(dims['n_frames'])

# Capture videos

cap_vid = cv2.VideoCapture(LABELS)
cap_path = cv2.VideoCapture(PATHS)

####################
# Get N frames of video, width, and height
####################

IN = {
    "vid": {"path" : LABELS},
    "path": {"path": PATHS}
    }

for key in IN:
    cap = cv2.VideoCapture(IN[key]['path'])
    IN[key]['n_frames'] = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    IN[key]['fps'] = int(cap.get(cv2.CAP_PROP_FPS))
    IN[key]['wid'] = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    IN[key]['hei'] = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))


# Throw error if n_frames aren't equal for all videos

if not (IN['vid']['n_frames'] == IN['path']['n_frames']):
    raise Exception('Videos do not have the same number of frames')

# Get height and total width (by doubling original video) in pixels 

TOT_WID = IN['vid']['wid']*2
TOT_HEI = IN['vid']['hei']

# Get FPS

FPS = IN['vid']['fps']

####################
# Set up video writer
####################

fourcc = cv2.VideoWriter_fourcc('h', '2', '6', '4')

video_writer = cv2.VideoWriter(
    OUT,
    fourcc,
    FPS,
    (TOT_WID, TOT_HEI),
    isColor = True
)

#######################
# Write to video
#######################

start = 0
end = IN['vid']['n_frames']
i = start
while i in range(start,end):
    # Set frames
    cap_vid.set(cv2.CAP_PROP_POS_FRAMES, i)
    cap_path.set(cv2.CAP_PROP_POS_FRAMES, i)
    # Capture frames
    ret_vid, frame_vid = cap_vid.read()
    ret_path, frame_path = cap_path.read()
    # Stack frames horizontally to create top row
    top = np.hstack((frame_vid, frame_path))
    # Write to video
    video_writer.write(top)
    # Add to counter
    i += 1

video_writer.release()

