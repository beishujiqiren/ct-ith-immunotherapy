import numpy as np
import pandas as pd
import SimpleITK as sitk
from pathlib import Path
from tqdm import tqdm
import config as C


def resample_to_isotropic(img, new_spacing=C.TARGET_SPACING, is_mask=False, default_value=0.0):
    orig_spacing = np.array(img.GetSpacing(), dtype=float)
    orig_size    = np.array(img.GetSize(),   dtype=float)
    new_size = np.round(orig_size * orig_spacing / np.array(new_spacing, dtype=float)).astype(int)
    new_size = [int(max(1, s)) for s in new_size]

    rs = sitk.ResampleImageFilter()
    rs.SetOutputSpacing(tuple(float(s) for s in new_spacing))
    rs.SetSize(new_size)
    rs.SetOutputDirection(img.GetDirection())
    rs.SetOutputOrigin(img.GetOrigin())
    rs.SetTransform(sitk.Transform())
    rs.SetDefaultPixelValue(float(default_value))
    rs.SetInterpolator(sitk.sitkNearestNeighbor if is_mask else sitk.sitkBSpline)
    return rs.Execute(img)


def main():
    C.PREPROC_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(C.SANITY_CSV)
    df_valid = df[df["Mask_Path"].astype(str).str.len() > 0].copy()
    print(f"Candidates: {len(df_valid)}")

    records = []
    for row in tqdm(df_valid.itertuples(index=False), total=len(df_valid), desc="Preprocessing"):
        pid       = row.Patient_ID
        img_path  = Path(str(row.Image_Path))
        mask_path = Path(str(row.Mask_Path))
        if not img_path.exists() or not mask_path.exists():
            tqdm.write(f"[skip] {pid} image or mask not found")
            continue

        try:
            img  = sitk.ReadImage(str(img_path), sitk.sitkFloat32)
            mask = sitk.Cast(sitk.ReadImage(str(mask_path)) > 0, sitk.sitkUInt8)

            img_rs  = resample_to_isotropic(img,  is_mask=False, default_value=C.AIR_HU)
            mask_rs = resample_to_isotropic(mask, is_mask=True,  default_value=0.0)

            if C.HU_CLIP is not None:
                img_rs = sitk.Clamp(img_rs, lowerBound=float(C.HU_CLIP[0]),
                                    upperBound=float(C.HU_CLIP[1]))

            img_arr  = sitk.GetArrayFromImage(img_rs)
            mask_arr = sitk.GetArrayFromImage(mask_rs)
            roi = img_arr[mask_arr > 0]
            voxels = int(roi.size)
            if voxels == 0:
                tqdm.write(f"[skip] {pid} empty ROI after resampling")
                continue

            sitk.WriteImage(img_rs,  str(C.PREPROC_DIR / f"{pid}_preprocessed.nii.gz"))
            sitk.WriteImage(mask_rs, str(C.PREPROC_DIR / f"{pid}_mask_preprocessed.nii.gz"))

            records.append({
                "Patient_ID":  pid,
                "Orig_Spacing": tuple(round(s, 3) for s in img.GetSpacing()),
                "New_Size":     tuple(int(s) for s in img_rs.GetSize()),
                "ROI_Voxels":   voxels,
                "ROI_Mean_HU":  round(float(roi.mean()), 2),
                "ROI_Std_HU":   round(float(roi.std()), 2),
            })
        except Exception as e:
            tqdm.write(f"[error] {pid}: {e}")

    pd.DataFrame(records).to_csv(C.PREPROC_SUMMARY, index=False, encoding="utf-8-sig")
    print(f"\nProcessed {len(records)} cases")
    print(f"Images/masks -> {C.PREPROC_DIR}")
    print(f"Summary -> {C.PREPROC_SUMMARY}")


if __name__ == "__main__":
    main()
