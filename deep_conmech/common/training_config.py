from dataclasses import dataclass

import psutil
from conmech.helpers.config import Config


@dataclass
class TrainingConfig(Config):

    DEVICE: str = "cpu"
    # torch.autograd.set_detect_anomaly(True)
    # print(numba.cuda.gpus)

    DATALOADER_WORKERS = 4
    GENERATION_WORKERS = 2

    ############

    TOTAL_MEMORY_GB = psutil.virtual_memory().total / 1024 ** 3
    TOTAL_MEMORY_LIMIT_GB = round(TOTAL_MEMORY_GB * 0.9, 2)
    GENERATION_MEMORY_LIMIT_GB = round((TOTAL_MEMORY_GB * 0.8) / GENERATION_WORKERS, 2)

    #####################
    TEST: bool = False

    NORMALIZE_ROTATE: bool = True
    ############

    TRAIN_SCALE = 1.0
    VALIDATION_SCALE = 1.0
    PRINT_SCALE = 1.0

    FINAL_TIME = 2 if TEST else 4 #8

    MESH_DENSITY = 8 if TEST else 16
    ADAPTIVE_TRAINING_MESH = False #True #############

    ############

    FORCES_RANDOM_SCALE = 10 ################################1.3
    OBSTACLE_ORIGIN_SCALE = 2.0 * TRAIN_SCALE

    DATA_ZERO_FORCES = 0.5
    DATA_ROTATE_VELOCITY = 0.5

    
    U_RANDOM_SCALE = 0.2
    V_RANDOM_SCALE = 2.5

    U_NOISE_GAMMA = 0.1
    U_IN_RANDOM_FACTOR = 0.005 * U_RANDOM_SCALE
    V_IN_RANDOM_FACTOR = 0.005 * V_RANDOM_SCALE

    ############

    DATA_FOLDER = f"{MESH_DENSITY}"
    PRINT_DATA_CUTOFF = 0.1

    ############

    VALIDATE_AT_MINUTES = 1 if TEST else 5

    DATASET = "scenarios"  # synthetic # scenarios
    L2_LOSS = True  #!#
    BATCH_SIZE = 128  #!#
    VALID_BATCH_SIZE = 128  #!#
    SYNTHETIC_BATCHES_IN_EPOCH = 2 if TEST else 64  # 512  #!#
    SYNTHETIC_SOLVERS_COUNT = BATCH_SIZE * SYNTHETIC_BATCHES_IN_EPOCH

    ############

    USE_DATASET_STATS = False
    INPUT_BATCH_NORM = True
    INTERNAL_BATCH_NORM = False
    LAYER_NORM = True

    DROPOUT_RATE = None  # 0.0  # 0.1  # 0.2  0.05
    SKIP = True
    # GRADIENT_CLIP = 10.0

    ATTENTION_HEADS = None  # None 1 3 5

    INITIAL_LR = 1e-3  # 1e-4 # 1e-5
    LR_GAMMA = 1.0  # 0.999
    FINAL_LR = 1e-6

    LATENT_DIM = 128
    ENC_LAYER_COUNT = 2
    PROC_LAYER_COUNT = 0
    DEC_LAYER_COUNT = 2
    MESSAGE_PASSES = 8  # 5 # 10
