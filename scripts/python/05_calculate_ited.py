import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from scipy.stats import entropy
from tqdm import tqdm
import config as C


def main():
    df = pd.read_csv(C.RF_DIR / "RF_AllPatients_merged.csv")
    print(f"Read: {len(df)} SVs, {df['Patient_ID'].nunique()} patients, {df.shape[1]} cols")

    feature_cols = [c for c in df.columns if c not in C.NON_FEATURE_COLS]

    rows = []
    for pid, group in tqdm(df.groupby("Patient_ID"), desc="iTED"):
        row = {"Patient_ID": pid}
        for feature in feature_cols:
            X = group[[feature]].replace([np.inf, -np.inf], np.nan).dropna().values
            if X.size == 0 or len(np.unique(X)) < 2:
                row[f"iTED_{feature}"] = 0.0
                continue

            best_bic, best_gmm = np.inf, None
            for K in C.K_RANGE:
                try:
                    gmm = GaussianMixture(n_components=K, random_state=C.RANDOM_STATE)
                    gmm.fit(X)
                    b = gmm.bic(X)
                    if b < best_bic:
                        best_bic, best_gmm = b, gmm
                except Exception:
                    continue
            if best_gmm is None:
                row[f"iTED_{feature}"] = 0.0
                continue

            labels = best_gmm.predict(X)
            probs = np.bincount(labels) / len(labels)
            row[f"iTED_{feature}"] = float(entropy(probs, base=2))
        rows.append(row)

    out = pd.DataFrame(rows)
    C.ITED_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(C.ITED_CSV, index=False)
    print(f"\niTED ({out.shape[0]} patients x {out.shape[1]-1} features) -> {C.ITED_CSV}")


if __name__ == "__main__":
    main()
