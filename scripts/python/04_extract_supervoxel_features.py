import time
import traceback
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor
from tqdm import tqdm
import config as C


def extract_patient_SV(patient_id):
    results, log_lines = [], []
    image_path = C.PREPROC_DIR / f"{patient_id}_preprocessed.nii.gz"
    mask_path  = C.SV_DIR / f"{patient_id}.SVmask.nii.gz"
    if not image_path.exists() or not mask_path.exists():
        return f"{patient_id} missing files"

    try:
        mask_array = sitk.GetArrayFromImage(sitk.ReadImage(str(mask_path)))
        labels = [int(l) for l in np.unique(mask_array) if l != 0]
        if not labels:
            return f"{patient_id} no valid SV"

        print(f"{patient_id} SV count: {len(labels)}")

        extractor = featureextractor.RadiomicsFeatureExtractor()
        for fc in C.FEATURE_CLASSES:
            extractor.enableFeatureClassByName(fc)
        extractor.settings.update(C.RADIOMICS_SETTINGS)

        image_sitk = sitk.ReadImage(str(image_path))
        for label in labels:
            nv = int(np.sum(mask_array == label))
            if nv < C.MIN_SV_VOXELS:
                continue
            sm_sitk = sitk.GetImageFromArray((mask_array == label).astype(np.uint16))
            sm_sitk.CopyInformation(image_sitk)
            try:
                fd = extractor.execute(image_sitk, sm_sitk)
                fd = {k: v for k, v in fd.items()
                      if any(fc in k for fc in C.FEATURE_CLASSES)}
                fd['Patient_ID']  = patient_id
                fd['SV_Label']    = label
                fd['Voxel_Count'] = nv
                results.append(fd)
            except Exception as e:
                log_lines.append(f"SV{label}: {e}\n{traceback.format_exc()}")
    except Exception as e:
        log_lines.append(f"{patient_id}: {e}\n{traceback.format_exc()}")

    if results:
        pd.DataFrame(results).to_csv(C.RF_DIR / f"{patient_id}_RF_AllSV.csv", index=False)
    if log_lines:
        with open(C.RF_DIR / f"{patient_id}_RF_log.txt", 'w', encoding='utf-8') as f:
            f.write("\n".join(log_lines))
    return f"{patient_id} done, SV: {len(results)}"


def main():
    C.RF_DIR.mkdir(parents=True, exist_ok=True)
    patient_ids = [f.name.replace('.SVmask.nii.gz', '')
                   for f in C.SV_DIR.glob("*.SVmask.nii.gz")]
    print(f"Patients: {len(patient_ids)}, extracting (binWidth=25, HU)...\n")

    max_workers = min(4, len(patient_ids)) if patient_ids else 1
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(extract_patient_SV, pid): pid for pid in patient_ids}
        for i, fut in enumerate(tqdm(as_completed(futures), total=len(futures),
                                     desc="Extracting"), 1):
            print(f"({i}/{len(futures)}) {fut.result()}")
    print(f"\nDone in {(time.time()-t0)/60:.2f} min")

    all_csv = list(C.RF_DIR.glob("*_RF_AllSV.csv"))
    if all_csv:
        df_all = pd.concat([pd.read_csv(f) for f in all_csv], ignore_index=True)
        out_path = C.RF_DIR / "RF_AllPatients_merged.csv"
        df_all.to_csv(out_path, index=False)
        print(f"Merged -> {out_path}")
    else:
        print("No feature CSV found.")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
