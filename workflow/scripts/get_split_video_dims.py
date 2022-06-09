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

import cv2 as cv

# Get variables

## Debug
IN = ["/hps/nobackup/birney/users/ian/MIKK_F0_tracking/split/open_field/20191114_1129_32-2_R_A/q1.avi",
      "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/split/open_field/20191114_1129_32-2_R_A/q2.avi",
      "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/split/open_field/20191114_1129_32-2_R_A/q3.avi"]

## True
IN = snakemake.input
OUT = snakemake.output[0]

# Sort file names

IN = sorted(IN)

# Write lines to file and close

file = open(OUT, 'w')

# write header
file.writelines("assay,sample,quadrant,wid,hei,n_frames\n")

for VID in IN:
    cap = cv.VideoCapture(VID)
    wid = str(int(cap.get(cv.CAP_PROP_FRAME_WIDTH)))
    hei = str(int(cap.get(cv.CAP_PROP_FRAME_HEIGHT)))
    n_frames = str(int(cap.get(cv.CAP_PROP_FRAME_COUNT)))
    # get quadrant
    quadrant = os.path.basename(VID).replace('.avi', '')
    # get sample
    sample = VID.split('/')[-2]
    # get assay
    assay = VID.split('/')[-3]
    # compose full line
    line_to_write = ','.join([assay,sample,quadrant,wid,hei,n_frames]) + '\n'
    file.writelines(line_to_write)

file.close()
