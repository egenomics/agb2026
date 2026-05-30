# metadata_obtain_clean

This folder contains the scripts used to **select, cleaan and label** the metadata from American Gut Project (AGP) before it enters the pipeline.

The process runs in two sequential setps, both implemented in `cleaning_metadata.py`.

---

## Files

| File | Description |
|------|-------------|
| `cleaning_metadata.py` | Main script: column subsetting + metadata cleaning + sample splitting (healthy/unhealthy) |
| `sample_information_from_prep_1834.tsv` | Raw meadata export from Qiita (Study ID: 10317, Analysis ID: 1834) |
| `library.xslx` | Variable dictionary: store all types of variables of metadata and its information | 

---

## Requirements

**Python ≥ 3.10** 

### Dependencies

```bash
conda install pandas openpyxl  
```
> `pandas` is required to read `sample_information_from_prem_1834.tsv`and its management
> `openpyxl`is required to real `library.xlsx`, the AGP data dictionary. 

---

## Usage 

Make sure both `sample_information_from_prep_1834q.tsv`and `library.slsx`are in the same directory as `cleaning_metadata.py`, then run:

```bash
python cleaning_metadata.py
```

---

## Steps

### 1) Column subsetting

Reads the raw Qiita export and retains only the columns could realistically be collected in a clinical setting (e.g., age, sex, BMI, dietary habits, antibiotic history and disgnosed medical contitions). THe selection is guided by `library.xlsx`information.

**Output:** `raw_metadata_1834.tsv`

### 2) Cleaning

Takes `raw_metadata_1834.tsv`and applies the following:

- Drops columns that exceed the missing data threshold (> 20% NaN), are entirely invalid (e.g. "not providede2, "not applicable") or are constant across all samples. Then save a summary of all dropped columns and the reason for removal

**Output:** `dropped_columns_summary_1834.tsv`

### 3) Sample splitting

Split samples into two groups:
 
- **Unhealthy**: individuals with an autoinmune disease diagnosed by a mediacl professional, without gut-related comorbidities (cancer, C.diff, diabetes, fungal overgrowth, IBD, IBS, kidney disease, liver diseade, SIBO) and without andibiotic use in the past 6 months.
- **Healthy**: individuals with a normal BMI, no diagnosed condition across all disease columns, and no antibiotic use in the past year.

The adds a `healthy`column (`yes`/`no`) according to the previous criteira. And rename `sample_name`column to `sample_name_id` which is necessay to have no confict to retrieve sra in next steps.

**Output:** `sample_information_cleaned_1834.tsv`

## Configuration

There are some key parameters that can be adjusted:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nan_threshold` | `0.20` | Maximum fraction of NaN values allowed per column |
| `invalid_threshold` | `None` | If `None`, only drops columns where **all** values are invalid |
| `invalid_values` | `["not provided", "not applicable", "n/a", "not collected"]` | Strings treated as missing |
| `required_columns` | `["sample_name", "sex", "age_cat", "bmi_cat", "autoimmune"]` | Columns that must be present in the input |


## Data source

- **Platform:** Qiita - Stydy ID:10317, Analysis ID: 1834
- **Reference:** McDonald, D. et al. (2018). American Gut: an Open Platform for Citizen Science Microbiome Research. mSystems, 3(3). https://doi.org/10.1128/mSystems.00031-18

## Next step

The output file `sample_information_cleaned_1834.tsv`is passed to `retrieve_sra.sh` to resolva NCBI SRA accession numbers for each sample.