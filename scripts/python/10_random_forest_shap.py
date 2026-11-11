"""Script 10B | Training-cohort random forest and validation-cohort SHAP.

Run from the repository root: python scripts/python/10_random_forest_shap.py
"""
from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score

try:
    import shap
except ImportError as exc:
    raise SystemExit("Install SHAP first: pip install shap") from exc

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "clinical"
RESULTS_DIR = ROOT / "results" / "clinical_association"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

TRAINING_FILE = DATA_DIR / "ITH_cohort1.csv"
VALIDATION_FILES = [DATA_DIR / "ITH_cohort2.csv", DATA_DIR / "ITH_cohort3.csv"]
OUTCOME_COLUMN = "CR"  
PATIENT_ID_COLUMN = "Patient_ID"  
FEATURES = [
    "iTED_original_firstorder_Minimum",
    "iTED_original_firstorder_RobustMeanAbsoluteDeviation",
    "iTED_original_glcm_ClusterProminence",
    "iTED_original_glcm_Contrast",
    "iTED_original_glszm_SmallAreaLowGrayLevelEmphasis",
]
RANDOM_STATE = 42
N_SPLITS = 5


def read_csv_robust(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required cohort file not found: {path}")
    last_error: Exception | None = None
    for encoding in ("utf-8-sig", "gb18030"):
        try:
            return pd.read_csv(path, encoding=encoding, low_memory=False)
        except (UnicodeDecodeError, LookupError) as exc:
            last_error = exc
    raise RuntimeError(f"Could not read {path} as UTF-8 or GB18030") from last_error


def encode_cr(series: pd.Series, source: str) -> pd.Series:
    values = series.astype("string").str.strip().str.lower()
    mapped = values.map({
        "cr": 1, "1": 1, "yes": 1, "true": 1,
        "non-cr": 0, "noncr": 0, "0": 0, "no": 0, "false": 0,
    })
    if mapped.isna().any():
        bad = sorted(values[mapped.isna()].dropna().unique().tolist())
        raise ValueError(f"Unrecognized {OUTCOME_COLUMN} values in {source}: {bad}")
    return mapped.astype(int)


def shap_values_for_positive_class(explainer, X: pd.DataFrame, positive_index: int) -> np.ndarray:
    """Normalize SHAP's binary-class output across supported SHAP versions."""
    raw = explainer.shap_values(X)
    if isinstance(raw, list):
        values = np.asarray(raw[positive_index])
    else:
        values = np.asarray(raw)
        if values.ndim == 3:
            
            values = values[:, :, positive_index]
        elif values.ndim != 2:
            raise ValueError(f"Unexpected SHAP output dimensions: {values.shape}")
    if values.shape != X.shape:
        raise ValueError(f"SHAP matrix shape {values.shape} does not match X {X.shape}")
    return values


train = read_csv_robust(TRAINING_FILE)
validation_parts = [read_csv_robust(p) for p in VALIDATION_FILES]
validation = pd.concat(validation_parts, ignore_index=True)

for name, frame in (("training", train), ("validation", validation)):
    needed = [OUTCOME_COLUMN, *FEATURES]
    missing = [col for col in needed if col not in frame.columns]
    if missing:
        raise KeyError(f"Missing required columns in {name} cohort: {missing}")
    frame[OUTCOME_COLUMN] = encode_cr(frame[OUTCOME_COLUMN], name)
    for feature in FEATURES:
        frame[feature] = pd.to_numeric(frame[feature], errors="coerce")
    if frame[FEATURES].isna().any().any():
        missing_counts = frame[FEATURES].isna().sum()
        raise ValueError(
            f"Missing/non-numeric iTED values in {name} cohort:\n"
            f"{missing_counts[missing_counts > 0].to_string()}\n"
            "Resolve missingness using the prespecified analysis rule; this script does not impute."
        )
    if frame[OUTCOME_COLUMN].nunique() != 2:
        raise ValueError(f"Both CR and Non-CR outcomes must be present in the {name} cohort.")

X_train = train[FEATURES].astype(float)
y_train = train[OUTCOME_COLUMN].astype(int)
X_valid = validation[FEATURES].astype(float)
y_valid = validation[OUTCOME_COLUMN].astype(int)


model = RandomForestClassifier(
    n_estimators=1000,
    max_depth=3,
    min_samples_split=2,
    min_samples_leaf=4,
    max_features=None,          
    criterion="gini",
    bootstrap=True,
    class_weight="balanced",
    random_state=RANDOM_STATE,
    n_jobs=-1,
)
cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
cv_auc = cross_val_score(model, X_train, y_train, cv=cv, scoring="roc_auc", n_jobs=-1)
cv_summary = pd.DataFrame({
    "fold": [str(i) for i in range(1, N_SPLITS + 1)] + ["mean"],
    "roc_auc": list(cv_auc) + [float(np.mean(cv_auc))],
})
cv_summary.to_csv(RESULTS_DIR / "random_forest_cv_auc.csv", index=False)

model.fit(X_train, y_train)
classes = list(model.classes_)
if 1 not in classes:
    raise RuntimeError("CR (coded as 1) is not present as the positive class.")
pos_idx = classes.index(1)
valid_prob_cr = model.predict_proba(X_valid)[:, pos_idx]
valid_auc = roc_auc_score(y_valid, valid_prob_cr)
pd.DataFrame([{
    "training_n": len(train),
    "validation_n": len(validation),
    "cv_mean_auc": float(np.mean(cv_auc)),
    "cv_sd_auc": float(np.std(cv_auc, ddof=1)),
    "validation_auc": float(valid_auc),
}]).to_csv(RESULTS_DIR / "random_forest_auc_summary.csv", index=False)

explainer = shap.TreeExplainer(model)
shap_matrix = shap_values_for_positive_class(explainer, X_valid, pos_idx)
shap_df = pd.DataFrame(shap_matrix, columns=[f"SHAP_CR_{f}" for f in FEATURES])
if PATIENT_ID_COLUMN in validation.columns:
    shap_df.insert(0, PATIENT_ID_COLUMN, validation[PATIENT_ID_COLUMN].astype(str).values)
else:
    warnings.warn(f"'{PATIENT_ID_COLUMN}' not found; SHAP rows will use validation row numbers.")
    shap_df.insert(0, "validation_row", np.arange(1, len(validation) + 1))
shap_df.to_csv(RESULTS_DIR / "validation_shap_values_cr.csv", index=False)

importance = pd.DataFrame({
    "feature": FEATURES,
    "mean_absolute_shap_cr": np.mean(np.abs(shap_matrix), axis=0),
}).sort_values("mean_absolute_shap_cr", ascending=False)
importance.to_csv(RESULTS_DIR / "validation_shap_global_importance.csv", index=False)

print(f"Training cohort: n={len(train)}; external validation cohort: n={len(validation)}")
print(f"5-fold stratified CV AUC: {np.mean(cv_auc):.3f} ± {np.std(cv_auc, ddof=1):.3f}")
print(f"External validation AUC: {valid_auc:.3f}")
print(f"Results saved under: {RESULTS_DIR}")
