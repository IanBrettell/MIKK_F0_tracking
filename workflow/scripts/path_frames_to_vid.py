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
from plotnine import *
import os
import shutil

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/merged/0.08.csv"
STITCHED = "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/stitched/0.08/open_field/20191113_1557_22-1_L_A.avi"
COLOURS_CSV = "/hps/software/users/birney/ian/repos/MIKK_F0_tracking/config/line_colours/line_colours_0.08.csv"
DIMS = "/hps/software/users/birney/ian/repos/MIKK_F0_tracking/config/split_video_dims.csv"
ASSAY = "open_field"
SAMPLE = "20191113_1557_22-1_L_A"
TMP = "/hps/nobackup/birney/users/ian/tmp"

## True
IN = snakemake.input.data[0]
STITCHED = snakemake.input.stitched_vid[0]
COLOURS_CSV = snakemake.input.colours
DIMS = snakemake.input.dims[0]
OUT = snakemake.output[0]
ASSAY = snakemake.params.assay
SAMPLE = snakemake.params.sample
TMP = snakemake.params.tmpdir

#######################
# Get FPS of stitched vid
#######################

cap = cv2.VideoCapture(STITCHED)
FPS = int(cap.get(cv2.CAP_PROP_FPS))

#######################
# Set plotting parameters
#######################

# Get colours

col_df = pd.read_csv(COLOURS_CSV)
pal_dict = dict(zip(list(col_df.line), list(col_df.colour)))
# Add iCab ref/test to dictionary
pal_dict['iCab_ref'] = '#F1BB7B'
pal_dict['iCab_test'] = '#AB7535'

REF = "iCab"
TEST = SAMPLE.split('_')[2]
# If it is an outbred fish, e.g. outbred-A-2, make it just "outbred"
if 'outbred' in TEST:
    TEST = 'outbred'

if REF == TEST:
    line_dict = {
        "ref" : "iCab_ref",
        "test" : "iCab_test"
    }
else:
    line_dict = {
        "ref" : "iCab",
        "test" : TEST
    }

#######################
# Read in tracking data
#######################

# Get date and time
DATE = SAMPLE.split('_')[0]
TIME = SAMPLE.split('_')[1]

# Read in file

df = pd.read_csv(IN)
# Convert time to string
df['time'] = df['time'].astype('string')
# Add a 0 to the start of `time` if only 3 characters
def add_0(x):
    if len(x) == 4:
        return x
    elif len(x) == 3:
        return "0" + x

df['time'] = df['time'].apply(add_0)

# Filter for assay and sample
df_filt = df.loc[(df['assay'] == ASSAY) & (df['date'] == int(DATE)) & (df['time'] == TIME)].copy()

# Create `line` column and recode

df_filt['line'] = df_filt['fish'].map(line_dict)

#######################
# Read in dims
#######################

dims = pd.read_csv(DIMS)
dims = dims.loc[(dims['sample'] == SAMPLE) & (dims['assay'] == ASSAY)]
N_FRAMES = max(dims['n_frames'])

# Set order of quadrants
df_filt['quadrant'] = pd.Categorical(df_filt['quadrant'], categories = ["q2", "q1", "q3", "q4"])

# Get max width and height

wid = dims['wid'].max()
hei = dims['hei'].max()

# Get total height and width in pixels (to match with the raw videos)

TOT_WID = dims.loc[dims['quadrant'].isin(['q1', 'q2'])]['wid'].sum() 
TOT_HEI = dims.loc[dims['quadrant'].isin(['q1', 'q4'])]['hei'].sum() 

####################
# Set up video writer
####################

fourcc = cv2.VideoWriter_fourcc('h', '2', '6', '4')

## Debug
#video_writer = cv2.VideoWriter(
#    "/hps/nobackup/birney/users/ian/pilot/tmp_out.avi",
#    fourcc,
#    FPS,
#    (TOT_WID, TOT_HEI),
#    isColor = True
#)

video_writer = cv2.VideoWriter(
    OUT,
    fourcc,
    FPS,
    (TOT_WID, TOT_HEI),
    isColor = True
)

#######################
# Plot frames and write to video
#######################

# Get all available frames
all_frames = list(range(1, N_FRAMES + 1))

# Get all frames we have
rec_frames = df_filt['frame'].unique().tolist()
rec_frames.sort()

for i in all_frames:
    # If the frame is not included in the frames we have data for...
    if i not in rec_frames:
        # get the next frame we do have
        filt_frames = []
        for frame in rec_frames:
            if frame > i:
                filt_frames.append(frame)
        plot_frame = min(filt_frames)
    elif i in rec_frames:
        plot_frame = i
    print(plot_frame)
    # If file already exists, write directly to file
    out_path = os.path.join(
        TMP,
        "MIKK_F0_tracking",
        "path_frames",
        ASSAY,
        SAMPLE, 
        str(plot_frame) + ".png"
    )
    # Make directory
    os.makedirs(os.path.dirname(out_path), exist_ok = True)
    # If the plot .png is already there, read it in and write
    if os.path.exists(out_path):
        # read plot
        frame = cv2.imread(out_path)
        # resize
        frame_out = cv2.resize(frame, (TOT_WID, TOT_HEI))
        # write
        video_writer.write(frame_out)
    # otherwise create the plot
    else:
        ## Filter df
        dat = df_filt.loc[df_filt['frame'] <= plot_frame]
        ## Plot
        plot = (ggplot(dat) +
                geom_path(aes('x', 'y', colour = 'line')) +
                facet_wrap('quadrant', nrow = 2) +
                scale_color_manual(pal_dict) +
                scale_x_continuous(limits = [0,wid]) +
                scale_y_reverse(limits = [hei,0]) +
                guides(color = None) +
                theme(
                    aspect_ratio = 1,
                    strip_background = element_blank(),
                    strip_text = element_blank(),
                    axis_line = element_blank(),
                    axis_text = element_blank(),
                    axis_title = element_blank(),
                    axis_ticks = element_blank(),
                    panel_background = element_blank(),
                    plot_background = element_rect(fill = "white")
                )
        )
        # Save
        #out_path = os.path.join(TMP, SAMPLE, REF_TEST, str(plot_frame) + ".png")
        plot.save(out_path, width = 9, height = 9)
        # Write to video
        ## read image
        frame = cv2.imread(out_path)
        # resize
        frame_out = cv2.resize(frame, (TOT_WID, TOT_HEI))
        # write
        video_writer.write(frame_out)
        # delete file
        os.remove(out_path)

video_writer.release()

# remove temp folder

shutil.rmtree(os.path.dirname(out_path))
