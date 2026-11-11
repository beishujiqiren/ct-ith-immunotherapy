# CT Habitat Radiomics Analysis
## Repository layout

```text
scripts/
  python/
    config.py
    01_check_ct_inputs.py
    02_preprocess_ct.py
    03_generate_supervoxels.py
    04_extract_supervoxel_features.py
    05_calculate_ited.py
  R/
    06_build_ith_phenotypes.R
    07_validate_and_visualize_phenotypes.R
    08_single_patient_supervoxel_heatmap.R
    09_baseline_table.R
    10_clinical_models_forest_and_subgroups.R
    11_sensitivity_and_model_diagnostics.R
    archive/earlier_clinical_model_version.R
```


## Inputs

Keep all patient-level files local. They are excluded from Git by `.gitignore`; do not commit identifiable imaging, clinical data, accession keys, or patient-level outputs.

- CT volumes and matching masks: `data/raw_ct/`. The Python check expects image files ending in `.nii` or `.nii.gz`, and a same-ID mask named `<ID>mask.nii` or `<ID>mask.nii.gz`.
- Phenotype-development clinical/features CSVs: `data/clinical/ITH_cohort1.csv`, `ITH_cohort2.csv`, and `ITH_cohort3.csv`. Script 06 treats cohort 1 as training and cohorts 2–3 as external validation.
- Clinical analysis files: `data/clinical/clinical_data_training.csv` and, for script 09, `data/clinical/clinical_data_validation_imputed.xlsx`.

The clinical scripts retain the source column names and coding. Check those names, category labels, and CR coding against your de-identified files before running them. R scripts can be run from the repository root. To run them from elsewhere, set `ITH_PROJECT_ROOT` to the repository root first.

## Run order

### CT habitat features

From the repository root, run the Python scripts in order:

```bash
python scripts/python/01_check_ct_inputs.py
python scripts/python/02_preprocess_ct.py
python scripts/python/03_generate_supervoxels.py
python scripts/python/04_extract_supervoxel_features.py
python scripts/python/05_calculate_ited.py
```

Outputs are written beneath `results/`: preprocessed CT and masks, supervoxel masks, per-supervoxel radiomics tables, the merged feature table, and `iTED_ITH_features.csv`. Python dependencies are listed in `requirements-python.txt`; the manuscript reports Python 3.10 and PyRadiomics 3.1.0.

