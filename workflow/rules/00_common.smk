######################
# Libraries
######################

import pandas as pd
import numpy as np
import os
import itertools

######################
# Config file and sample sheets
######################

configfile: "config/config.yaml"

# Read in samples file
samples_df = pd.read_csv(config["samples_file"], comment = '#')

# Set variables

SAMPLES = samples_df["sample"]
ASSAYS = ["open_field", "novel_object"]
QUADRANTS = ["q1", "q2", "q3", "q4"]

# The novel_object assay is missing from 20191118_1311_80-2_L_B
# Therefore it needs to be removed from the combinations

## Create list of variable lists
full_list = [SAMPLES, ASSAYS, QUADRANTS]
## Create list of tuple combinations
combos = list(itertools.product(*full_list))
## Remove unavailable combinations
for QUADRANT in QUADRANTS:
    combos.remove(('20191118_1311_80-2_L_B', 'novel_object', QUADRANT))
## Create new lists of variables
SAMPLES = [i[0] for i in combos]
ASSAYS = [i[1] for i in combos]
QUADRANTS = [i[2] for i in combos]