from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw_ct"
WORK_ROOT = PROJECT_ROOT / "results"

SANITY_CSV = WORK_ROOT / "sanity_summary.csv"
PREPROC_DIR = WORK_ROOT / "preprocessed_ct"
SV_DIR = WORK_ROOT / "supervoxels"
RF_DIR = WORK_ROOT / "supervoxel_features"
PREPROC_SUMMARY = WORK_ROOT / "ct_preprocessing_summary.csv"
ITED_CSV = WORK_ROOT / "iTED_ITH_features.csv"

TARGET_SPACING = (1.0, 1.0, 1.0)
HU_CLIP = (-200.0, 350.0)
AIR_HU = -1024.0
FEATURE_CLASSES = ["firstorder", "glcm", "glrlm", "glszm", "gldm", "ngtdm"]
MIN_SV_VOXELS = 10
RADIOMICS_SETTINGS = {
    "binWidth": 25,
    "normalize": False,
    "resampledPixelSpacing": None,
    "interpolator": "sitkBSpline",
    "label": 1,
}
K_RANGE = range(1, 6)
RANDOM_STATE = 42
NON_FEATURE_COLS = ["Patient_ID", "SV_ID", "SV_Label", "Voxel_Count"]
