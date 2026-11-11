import logging
import numpy as np
import SimpleITK as sitk
from skimage.segmentation import slic
from skimage.morphology import ball, opening
from tqdm import tqdm
import config as C


def main():
    C.SV_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=str(C.SV_DIR.parent / "slic.log"), level=logging.INFO,
                        format="%(asctime)s | %(levelname)s | %(message)s")
    logging.getLogger().addHandler(logging.StreamHandler())

    vols = sorted([f for f in C.PREPROC_DIR.glob("*_preprocessed.nii.gz")
                   if "_mask_" not in f.name])
    print(f"To segment: {len(vols)}")

    for vol_path in tqdm(vols, desc="SLIC"):
        pid = vol_path.name.replace("_preprocessed.nii.gz", "")
        mask_path = C.PREPROC_DIR / f"{pid}_mask_preprocessed.nii.gz"
        out_path  = C.SV_DIR / f"{pid}.SVmask.nii.gz"
        if not mask_path.exists():
            logging.warning(f"{pid} missing mask")
            continue

        try:
            img_sitk = sitk.ReadImage(str(vol_path))
            img = sitk.GetArrayFromImage(img_sitk)
            denom = img.max() - img.min()
            img_norm = (img - img.min()) / denom if denom > 0 else img.astype(np.float32)
            img_norm = np.sqrt(img_norm).astype(np.float32)

            mask = sitk.GetArrayFromImage(sitk.ReadImage(str(mask_path)))
            mask_clean = opening(mask, ball(1))

            roi_voxels = img_norm[mask_clean > 0]
            if roi_voxels.size == 0:
                logging.warning(f"{pid} empty ROI")
                continue
            roi_size = roi_voxels.size
            roi_std  = float(np.std(roi_voxels))

            if roi_size < 500:
                n_segments = 20
            elif roi_size < 2000:
                n_segments = 40
            else:
                n_segments = 80
            if roi_std < 0.05:
                compactness = 0.1
            elif roi_std < 0.1:
                compactness = 0.3
            else:
                compactness = 1.0

            logging.info(f"{pid} | ROI={roi_size}, std={roi_std:.4f}, "
                         f"n_segments={n_segments}, compactness={compactness}")

            slic_mask = slic(img_norm, n_segments=n_segments, compactness=compactness,
                             mask=mask_clean, start_label=1, channel_axis=None)
            out_sitk = sitk.GetImageFromArray(np.uint16(slic_mask))
            out_sitk.CopyInformation(img_sitk)
            sitk.WriteImage(out_sitk, str(out_path))
        except Exception as e:
            logging.error(f"{pid} SLIC failed: {e}")

    print(f"\nSV masks -> {C.SV_DIR}")


if __name__ == "__main__":
    main()
