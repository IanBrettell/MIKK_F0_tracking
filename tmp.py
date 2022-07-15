bsub -Is bash
conda activate snakemake_6.15.5
cd /hps/software/users/birney/ian/repos/MIKK_F0_tracking
python

###################

import pandas as pd
samples_df = pd.read_csv('config/samples_long.csv')

samples_df = pd.read_csv(config["samples_long"])

# Set variables

SAMPLES_ZIP = samples_df['sample']
ASSAYS_ZIP = samples_df['assay']
QUADRANTS_ZIP = samples_df['quadrant']

# Get lists of unique assay/sample (AS) combinations

sa_zip = sorted(list(set(zip(SAMPLES_ZIP, ASSAYS_ZIP))))
AS_SAMPLES = [i[0] for i in sa_zip]
AS_ASSAYS = [i[1] for i in sa_zip]

# All quadrants

QUADRANTS_ALL = ["q1", "q2", "q3", "q4"]

# Get samples with over 85% tracking success

## Create data frame with all
zip_df = pd.DataFrame(
    {
        'sample': SAMPLES_ZIP,
        'assay': ASSAYS_ZIP,
        'quadrant' : QUADRANTS_ZIP
    }
)
## Read in .csv
ts_df = pd.read_csv('config/tracking_success.csv')
## Filter for samples with < 85% tracking success
ts_df = ts_df.loc[ts_df['prop_success'] < 0.85]
## Remove samples with 
for i, row in ts_df.iterrows():
    target_assay = row['assay']
    target_sample = row['sample']
    target_quadrant = row['quadrant']
    zip_df.drop(
        zip_df[
            (zip_df['sample'] == target_sample) & \
            (zip_df['assay'] == target_assay) & \
            (zip_df['quadrant'] == target_quadrant)
            ].index,
            inplace = True
    )
## Get filtered samples and assays
SAMPLES_ZIP_TRK = zip_df['sample']
ASSAYS_ZIP_TRK = zip_df['assay']
QUADRANTS_ZIP_TRK = zip_df['quadrant']


#################################

assay='novel_object'
sample='20191118_1311_80-2_L_B'
quadrant='q1'

# get_vid_length
row = samples_df.loc[(samples_df['sample'] == sample) & \
                        (samples_df['assay'] == assay) & \
                        (samples_df['quadrant'] == quadrant)]
if assay == "open_field":
    start = int(row['of_start'])
    end = int(row['of_end'])
elif assay == "novel_object":
    start = int(row['no_start'])
    end = int(row['no_end'])
vid_length = int(end) - int(start)

#get_bgsub(wildcards):
bgsub = samples_df.loc[(samples_df['sample'] == sample) & \
                        (samples_df['assay'] == assay) & \
                        (samples_df['quadrant'] == quadrant), \
                        'bgsub'].values[0]

##def get_intensity_floor(wildcards):
int_floor = samples_df.loc[(samples_df['sample'] == sample) & \
                        (samples_df['assay'] == assay) & \
                        (samples_df['quadrant'] == quadrant), \
                        'intensity_floor'].values[0]

#def get_intensity_ceiling(wildcards):
int_ceiling = samples_df.loc[(samples_df['sample'] == sample) & \
                                (samples_df['assay'] == assay) & \
                                (samples_df['quadrant'] == quadrant), \
                                'intensity_ceiling'].values[0]

#def get_area_floor(wildcards):
area_floor = samples_df.loc[(samples_df['sample'] == sample) & \
                            (samples_df['assay'] == assay) & \
                            (samples_df['quadrant'] == quadrant), \
                            'area_floor'].values[0]

#def get_area_ceiling(wildcards):
area_ceiling = samples_df.loc[(samples_df['sample'] == sample) & \
                            (samples_df['assay'] == assay) & \
                            (samples_df['quadrant'] == quadrant), \
                            'area_ceiling'].values[0]
