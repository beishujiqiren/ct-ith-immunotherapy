import pandas as pd
import config as C

records = []
for img_file in C.RAW_DIR.glob("*.nii*"):
    if "mask" in img_file.name.lower():
        continue
    pid = img_file.name.replace(".nii.gz", "").replace(".nii", "")
    mask_file = C.RAW_DIR / f"{pid}mask.nii.gz"
    if not mask_file.exists():
        mask_file = C.RAW_DIR / f"{pid}mask.nii"
    records.append({
        "Patient_ID": pid,
        "Status": "Complete" if mask_file.exists() else "Missing mask",
        "Image_Path": str(img_file),
        "Mask_Path": str(mask_file) if mask_file.exists() else "",
    })

summary = pd.DataFrame(records)
C.SANITY_CSV.parent.mkdir(parents=True, exist_ok=True)
summary.to_csv(C.SANITY_CSV, index=False, encoding="utf-8-sig")

n_ok = (summary["Status"] == "Complete").sum()
print(summary.to_string(index=False))
print(f"\nTotal {len(summary)}, complete {n_ok}, missing mask {len(summary)-n_ok}")
print(f"Saved: {C.SANITY_CSV}")
