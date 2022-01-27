# GUI defaults
NUMBER_OF_ANIMALS_DEFAULT = 2
MIN_THRESHOLD_DEFAULT, MAX_THRESHOLD_DEFAULT = 0, 200
AREA_LOWER, AREA_UPPER = 0, 60000
MIN_AREA_DEFAULT, MAX_AREA_DEFAULT = 150, 60000
RES_REDUCTION_DEFAULT = 1.0

# Parallel processes
NUMBER_OF_JOBS_FORBACKGROUND_SUBTRACTION = 1
NUMBER_OF_JOBS_FOR_SEGMENTATION = 1

# Knowledge transfer and identity transfer
KNOWLEDGE_TRANSFER_FOLDER_IDCNN = None
IDENTITY_TRANSFER = False
IDENTIFICATION_IMAGE_SIZE = None

# Memory management
SAVE_PIXELS = 'DISK' # 'RAM' or 'NOT'
SAVE_SEGMENTATION_IMAGE = 'DISK' # 'RAM' or 'NOT'

# Data management policy
DATA_POLICY = 'all' # 'trajectories', 'validation', 'knowledge_transfer', 'idmatcher.ai' or 'all'_
