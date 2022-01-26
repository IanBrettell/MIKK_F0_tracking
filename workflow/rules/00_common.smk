######################
# Libraries
######################

import pandas as pd
import numpy as np
import os

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
