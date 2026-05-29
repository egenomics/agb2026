# %% [markdown]
# # Step 1: Subset Metadata
# 
# In this first part of the notebook we will **build a subset of the metadata** only using a selected group of columns that we consider that can be obtained as clinical data in a hospital.
# 
# * The **output file** of this section is the supposed ***raw_metadata.tsv*** that we would be collected in the hospital by the doctor, for example.
# 
# * To filter the columns we have used the *description_file* of the American Gut Project.
# 
# 
# For this pipeline, we are using a curated subset of the American Gut Project (AGP) dataset.
# 
#     Source Platform: Qiita (Study ID: 10317)
# 
#     Target Subset: 649 samples, retrieved via Qiita. Analysis ID: 1834.
# 
# For **executing** this pipeline, **we need** to have in the same folder as this script the file: ***sample_information_1834.tsv***, that can be obtained from Qiita.
# 
# 
# **Reference Publication:**
# McDonald, D. et al. (2018). "American Gut: an Open Platform for Citizen Science Microbiome Research". mSystems, 3(3). DOI: https://doi.org/10.1128/mSystems.00031-18

# %%
#Load libraries
import pandas as pd
import logging #for logging messages about the cleaning process (e.g., which columns were dropped and why)
from pathlib import Path
from dataclasses import dataclass, field #for creating configuration classes with default values
from typing import Optional #for type hinting of optional parameters
import re #for regular expressions

# Structured logging (shows timestamp, level and message)
logging.basicConfig(
    level=logging.INFO, #minimum messages showed (INFO, WARNING, ERROR, CRITICAL)
    format="%(asctime)s [%(levelname)s] %(message)s", #timestamp, log level and message
    datefmt="%Y-%m-%d %H:%M:%S", #reading format for the timestamp
)
logger = logging.getLogger(__name__) #the name of the logger is the name of the module (file) where it is used


# %% [markdown]
# 

# %% [markdown]
# ## Configuration
# 
# Here we specify the tunable parameters.

# %%
@dataclass
class RawMetadataConfig:
    """Central configuration for the metadata cleaning pipeline."""

    # File paths
    AGP_metadata_file: Path = Path("/sample_information_from_prep_1834.tsv") #input file with the metadata to be cleaned
    raw_metadata_file: Path = Path("raw_metadata_1834.tsv") #output file with the raw metadata (after initial loading and basic cleaning)

    #Columns to mantain in the metadata file
    columns_to_keep_not_medical_condition: list = field(default_factory=lambda: ["exercise_frequency", "probiotic_frequency", "env_material", "host_taxid", "taxon_id", "age_cat", "antibiotic_history", "appendix_removed", "bmi_cat", "cancer_treatment", "chickenpox", "contraceptive", "country_of_birth", "country_residence", "csection", "diabetes_type", "diet_type", "fed_as_infant", "flu_vaccine_date", "gluten", "height_cm", "ibd_diagnosis", "lactose", "pregnant", "race", "sample_name", "sex", "tonsils_removed", "weight_kg"])

    columns_to_keep_medical_condition: list = field(default_factory=lambda: ["acid_reflux", "add_adhd", "alzheimers", "asd", "autoimmune", "cancer", "cardiovascular_disease", "cdiff", "clinical_condition", "depression_bipolar_schizophrenia", "diabetes", "epilepsy_or_seizure_disorder", "fungal_overgrowth", "ibd",  "ibs", "kidney_disease", "liver_disease", "lung_disease", "mental_illness", "mental_illness_type_anorexia_nervosa", "mental_illness_type_bipolar_disorder", "mental_illness_type_bulimia_nervosa", "mental_illness_type_depression", "mental_illness_type_schizophrenia",  "mental_illness_type_substance_abuse", "migraine", "pku", "sibo", "skin_condition", "thyroid"])

    # Post-cleaning expectations
    #expected_max_columns = int(sum(len(columns_to_keep_not_medical_condition), len(columns_to_keep_medical_condition))) #maximum number of columns expected after cleaning


# %% [markdown]
# ## Removing columns not required

# %%
def good_metadata(cfg: RawMetadataConfig):
    """Cleans the metadata file according to the specified configuration."""

    # Load the raw metadata file
    df = pd.read_csv(cfg.AGP_metadata_file, sep="\t")
    logger.info("Metadata loaded from %s with %d samples and %d columns", cfg.AGP_metadata_file, df.shape[0], df.shape[1])

    # Filter the dataframe to keep only the specified columns
    columns_to_keep = cfg.columns_to_keep_not_medical_condition + cfg.columns_to_keep_medical_condition
    df_clean = df[columns_to_keep]
    logger.info("Filtered metadata to keep specified columns. Remaining columns: %d", len(df_clean.columns))

    # Save the raw metadata file (after initial loading and basic cleaning)
    df_clean.to_csv(cfg.raw_metadata_file, sep="\t", index=False)
    logger.info("Raw metadata saved to %s", cfg.raw_metadata_file)

# %%
# Pipeline execution
if __name__ == "__main__":
    cfg = RawMetadataConfig()
    good_metadata(cfg)


# %% [markdown]
# # Step 2: Metadata Cleaning
# 
# **Outputs**:
# * *sample_information_cleaned_1834.tsv*: file that contains the clean information with just the samples that passed the filters.
# * *dropped_columns_summary_1834.tsv*: file that contains the removed columns from the *raw_metadata.tsv* file.

# %%
#Load libraries
import pandas as pd
import logging #for logging messages about the cleaning process (e.g., which columns were dropped and why)
from pathlib import Path
from dataclasses import dataclass, field #for creating configuration classes with default values
from typing import Optional #for type hinting of optional parameters
import re #for regular expressions

# Structured logging (shows timestamp, level and message)
logging.basicConfig(
    level=logging.INFO, #minimum messages showed (INFO, WARNING, ERROR, CRITICAL)
    format="%(asctime)s [%(levelname)s] %(message)s", #timestamp, log level and message
    datefmt="%Y-%m-%d %H:%M:%S", #reading format for the timestamp
)
logger = logging.getLogger(__name__) #the name of the logger is the name of the module (file) where it is used

# %% [markdown]
# ## Configuration
# 
# Here we can specify the tunable parameters.
# 
# **For invalid values, columns are only removed if all entries are invalid. This is done to avoid discarding variables that may distinguish controls from patients or capture relevant conditions important for downstream classification.**
# 
# For other datasets, these parameters can be adjusted and refined in the next step of the pipeline.

# %%
@dataclass
class MetadataConfig:
    """Central configuration for the metadata cleaning pipeline."""

    # File paths
    dictionary: Path = Path("library.xlsx")
    input_file: Path = Path("raw_metadata_1834.tsv")
    clean_file: Path = Path("sample_information_cleaned_1834.tsv")
    summary_file: Path = Path("dropped_columns_summary_1834.tsv")
    #prep_info_file: Path = Path("10317_prep_1834_20210907-154515.txt")

    # Dropping thresholds
    nan_threshold: float = 0.20          # Drop columns with > this fraction of NaN
    invalid_threshold: Optional[float] = None  # None = only drop if ALL values invalid

    # Domain-specific invalid values
    # This creates a new list for each instance of MetadataConfig with the specified default values
    invalid_values: list = field(default_factory=lambda: [
        "not provided",
        "not applicable",
        "n/a",
        "not collected"
    ])

    # Pre-cleaning expectations (set to None to skip)
    expected_min_rows: Optional[int] = None     # Dataset must have at least this many rows (are there enough samples?)
    expected_min_cols: Optional[int] = 5       # Dataset must have at least this many cols (are the required ones)
    required_columns: list = field(default_factory=lambda: [
        # Written in normalized snake_case (post-normalization names)
         "sample_name", "sex","age_cat", "bmi_cat", "autoimmune"
    ])

    # Disease columns to keep
    columns_to_keep_medical_condition: list = field(default_factory=lambda: ["acid_reflux", "add_adhd", "alzheimers", "asd", "autoimmune", "cancer", "cardiovascular_disease", "cdiff", "clinical_condition", "depression_bipolar_schizophrenia", "diabetes", "epilepsy_or_seizure_disorder", "fungal_overgrowth", "ibd",  "ibs", "kidney_disease", "liver_disease", "lung_disease", "mental_illness", "mental_illness_type_anorexia_nervosa", "mental_illness_type_bipolar_disorder", "mental_illness_type_bulimia_nervosa", "mental_illness_type_depression", "mental_illness_type_schizophrenia",  "mental_illness_type_substance_abuse", "migraine", "pku", "sibo", "skin_condition", "thyroid"])

    # Post-cleaning expectations
    max_allowed_nan_frac: float = 0.20        # No column in output should exceed this
    min_retained_cols: Optional[int] = None    # Output must retain at least this many cols (I don't know what a good number is, so I'm leaving it as None for now)

cfg = MetadataConfig()
logger.info("Configuration loaded.")

# %% [markdown]
# ## Column name normalization
# 
# Auto-rule applied to every column.
# 
# In this case all columns are in snake_case lowercase but we take into account all possible cases.

# %%
def normalize_column_names(df: pd.DataFrame, cfg: MetadataConfig) -> pd.DataFrame:
    """
    Homogenize column names to snake_case lowercase.
    """

    def auto_snake(name: str) -> str:
        name = name.strip()
        name = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", name)  # SRAStudy -> SRA_Study
        name = re.sub(r"([a-z])([A-Z])", r"\1_\2", name)        # camelCase -> camel_Case
        name = re.sub(r"[\s\-]+", "_", name)                    # spaces or - -> _
        name = re.sub(r"[^\w]", "", name)                       # remove special chars
        name = re.sub(r"_+", "_", name)                         # collapse multiple _
        return name.lower()


    logger.info("Column names normalized.")
    return df

# %% [markdown]
# ## Validation helpers
# 
# We have different validators:
# * One for the raw input
# * One for the cleaned output.
# 
# Both of them return a list of issues. An empty list means that all checks are passed.
# 
# (Validators recieve the dataframe after normalization, so column names are already in snake_case).

# %%
def validate_input(df: pd.DataFrame, cfg: MetadataConfig) -> list[str]:
    """
    Pre-cleaning checks on the raw dataframe.
    Returns a list of issue strings (empty = all OK).
    """
    issues = []

    #Check if the dataset meets the minimum row and column requirements (if specified in the config)
    if cfg.expected_min_rows and len(df) < cfg.expected_min_rows:
        issues.append(
            f"Too few rows: got {len(df)}, expected >= {cfg.expected_min_rows}"
        )

    if cfg.expected_min_cols and df.shape[1] < cfg.expected_min_cols:
        issues.append(
            f"Too few columns: got {df.shape[1]}, expected >= {cfg.expected_min_cols}"
        )

    # Check for required columns (after normalization)
    missing_required = [c for c in cfg.required_columns if c not in df.columns]
    if missing_required:
        issues.append(f"Required columns missing: {missing_required}")

    #Reporting duplicate rows
        n_dup = df.duplicated().sum()
        issues.append(f"Duplicate rows detected: {n_dup}")

    return issues


def validate_output(df_clean: pd.DataFrame, df_raw: pd.DataFrame, cfg: MetadataConfig) -> list[str]:
    """
    Post-cleaning checks on the cleaned dataframe.
    Returns a list of issue strings (empty = all OK).
    """
    issues = []

    # Row count must be preserved
    if len(df_clean) != len(df_raw):
        issues.append(
            f"Row count changed during cleaning: {len(df_raw)} -> {len(df_clean)}"
        )


    # No column should exceed the NaN threshold
    nan_fracs = df_clean.isna().mean()
    bad_cols = nan_fracs[nan_fracs > cfg.max_allowed_nan_frac].index.tolist()
    if bad_cols:
        issues.append(
            f"Columns exceeding NaN threshold in output: {bad_cols}"
        )

    # Minimum number of columns retained
    if cfg.min_retained_cols and df_clean.shape[1] < cfg.min_retained_cols:
        issues.append(
            f"Too few columns retained: {df_clean.shape[1]} < {cfg.min_retained_cols}"
        )

    return issues


def report_validation(issues: list[str], stage: str) -> None:
    """Log validation results. Raises if issues found."""
    if not issues:
        logger.info("[%s] All checks passed.", stage)
    else:
        for issue in issues:
            logger.warning("[%s] %s", stage, issue)
        raise ValueError(
            f"{stage} validation failed with {len(issues)} issue(s). See warnings above."
        )

# %% [markdown]
# ## Cleaning function
# 

# %%
def clean_metadata(df: pd.DataFrame, cfg: MetadataConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Apply general cleaning criteria to a metadata dataframe.

    Drops columns that are:
      - Entirely NaN
      - Above the NaN threshold
      - Entirely invalid values (or above invalid threshold if set)
      - Constant (single unique non-NaN value)

    Returns
    -------
    df_clean : cleaned dataframe
    summary  : dataframe with one row per dropped column and reason
    """
    invalid_set = {v.lower() for v in cfg.invalid_values}
    nan_limit = int(cfg.nan_threshold * len(df))   # max allowed NaN count

    dropped: dict[str, str] = {}

    for col in df.columns:
        series = df[col]
        n_nan = series.isna().sum()

        # All NaN
        if n_nan == len(df):
            dropped[col] = "All values are NaN"
            continue

        # Above NaN threshold
        if n_nan > nan_limit:
            dropped[col] = f">{cfg.nan_threshold:.0%} NaN ({n_nan}/{len(df)})"
            continue

        # Invalid values (normalise to lowercase for comparison)
        is_invalid = series.dropna().astype(str).str.lower().isin(invalid_set)
        n_invalid = is_invalid.sum()
        n_non_nan = len(series.dropna())

        if cfg.invalid_threshold is None:
            # Drop only if ALL non-NaN values are invalid
            if n_non_nan > 0 and n_invalid == n_non_nan:
                dropped[col] = f"All non-NaN values are invalid ({n_invalid} values)"
                continue
        else:
            if n_non_nan > 0 and (n_invalid / n_non_nan) > cfg.invalid_threshold:
                dropped[col] = (
                    f">{cfg.invalid_threshold:.0%} invalid values "
                    f"({n_invalid}/{n_non_nan})"
                )
                continue

        # 4. Constant column
        if series.dropna().nunique() == 1:
            dropped[col] = f"Constant value: '{series.dropna().iloc[0]}'"
            continue

    df_clean = df.drop(columns=list(dropped.keys()))

    # summary
    summary = pd.DataFrame(
        list(dropped.items()), columns=["Column", "Reason"]
    )

    logger.info(
        "Cleaning complete: %d -> %d columns (%d dropped).",
        df.shape[1], df_clean.shape[1], len(dropped),
    )
    return df_clean, summary

# %% [markdown]
# ## Run pipeline
# 

# %%
# Load files
df_raw = pd.read_csv(cfg.input_file, sep=None, engine="python") # automatically detects separator
logger.info("Raw metadata loaded: %s", df_raw.shape)

# Normalize column names
df_raw = normalize_column_names(df_raw, cfg)

# Pre-cleaning validation
pre_issues = validate_input(df_raw, cfg)
report_validation(pre_issues, stage="PRE-CLEANING")

# Clean
df_clean, summary = clean_metadata(df_raw, cfg)

# Post-cleaning validation
post_issues = validate_output(df_clean, df_raw, cfg)
report_validation(post_issues, stage="POST-CLEANING")

#Save
df_clean.to_csv(cfg.clean_file, sep="\t", index=False)
logger.info("Cleaned metadata saved to: %s", cfg.clean_file)

summary.to_csv(cfg.summary_file, sep="\t", index=False)
logger.info("Dropped columns summary saved to: %s", cfg.summary_file)

logger.info("Pipeline finished. Final shape: %s", df_clean.shape)

# %% [markdown]
# ## Split into two datasets
# 
# We are dividing samples between controls ("healthy people") and cases (others).
# 
# Selection of **controls** with these parameters:
# - Using the study variable: *bmi_cat* = 18.5- 30 kg/m2
# - Removing patients with a set of diseases and those who had taken antibiotics in the past year (defined in df_healthy).
# 
# 
# 
# Selection of **cases**:
# - Patients with an **autoimmune disease** ("Diagnosed by a medical professional (doctor, physician assistant)")
# - We have removed autoimmune patients with diseases highly related to gut microbiome and patients who had taken antibiotics in the past 6 months (defined in df_disease).

# %%
# Constants
diagnosed_by_doctor = "Diagnosed by a medical professional (doctor, physician assistant)"
not_taken_antibiotics = "I have not taken antibiotics in the past year."

def no_disease_mask(df: pd.DataFrame, disease_cols: list[str], exclude: list[str] | None = None) -> pd.Series:
    cols = [c for c in disease_cols if c not in (exclude or [])]
    return ~df[cols].eq(diagnosed_by_doctor).any(axis=1)

# Subgroups
disease_cols = cfg.columns_to_keep_medical_condition

mask_no_antibiotics = df_clean["antibiotic_history"] == not_taken_antibiotics
mask_autoimmune     = df_clean["autoimmune"] == diagnosed_by_doctor
mask_healthy_bmi    = df_clean["bmi_cat"] == "Normal"

logger.info("Samples with only autoimmune diagnosis: %d", mask_autoimmune.sum())
logger.info("Samples with normal BMI: %d", mask_healthy_bmi.sum())

df_disease = df_clean[
    mask_autoimmune &
    no_disease_mask(df_clean, disease_cols, exclude=["autoimmune"]) &
    mask_no_antibiotics
]

df_healthy = df_clean[
    mask_healthy_bmi &
    mask_no_antibiotics &
    no_disease_mask(df_clean, disease_cols)
]

# logger.info("Autoimmune group: %d samples", len(df_disease))
logger.info("Healthy controls: %d samples", len(df_healthy))
# print(f"Autoimmune group : {len(df_disease):>6} samples")
# print(f"Healthy controls : {len(df_healthy):>6} samples")


# %%
# Reference dataframe that contains the filters based on the information provided in the paper
df_clean = pd.read_csv(cfg.clean_file, sep="\t")

df_disease = df_clean[
        #(df_clean["acid_reflux"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["add_adhd"]!= "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["alzheimers"] != "Diagnosed   by a medical professional (doctor, physician assistant)") &
        #(df_clean["asd"] != "Diagnosed by a medical professional (doctor, physician assistant)")&
        (df_clean["autoimmune"] == "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["cancer"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["cardiovascular_disease"]!= "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["cdiff"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["clinical_condition"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
       #(df_clean["depression_bipolar_schizophrenia"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["diabetes"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["epilepsy_or_seizure_disorder"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["fungal_overgrowth"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["ibd"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["ibs"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["kidney_disease"] != "Diagnosed by a medical professional (doctor, physician assistant)")  &
        (df_clean["liver_disease"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["lung_disease"] != "Diagnosed by a medical professional (doctor, physician assistant)")  &
        #(df_clean["mental_illness"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["mental_illness_type_anorexia_nervosa"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["mental_illness_type_bipolar_disorder"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["mental_illness_type_bulimia_nervosa"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["mental_illness_type_depression"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["mental_illness_type_schizophrenia"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["mental_illness_type_substance_abuse"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["migraine"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["pku"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["sibo"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["skin_condition"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        #(df_clean["thyroid"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        df_clean["antibiotic_history"].isin(["I have not taken antibiotics in the past year.", "6 months"]) &
        (df_clean["env_material"] == "feces")
      ]

print('Filtered autoinmune cases:', len(df_disease))

# Healthy subset
df_healthy = df_clean[
        #(df_clean["age_cat"].isin(["20s", "30s", "40s", "50s", "60s"])) &
        (df_clean["bmi_cat"].isin(["Normal"])) &
        (df_clean["acid_reflux"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["add_adhd"]!= "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["alzheimers"] != "Diagnosed   by a medical professional (doctor, physician assistant)") &
        (df_clean["asd"] != "Diagnosed by a medical professional (doctor, physician assistant)")&
        (df_clean["autoimmune"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["cancer"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["cardiovascular_disease"]!= "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["cdiff"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["clinical_condition"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["depression_bipolar_schizophrenia"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["diabetes"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["epilepsy_or_seizure_disorder"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["fungal_overgrowth"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["ibd"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["ibs"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["kidney_disease"] != "Diagnosed by a medical professional (doctor, physician assistant)")  &
        (df_clean["liver_disease"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["lung_disease"] != "Diagnosed by a medical professional (doctor, physician assistant)")  &
        (df_clean["mental_illness"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["mental_illness_type_anorexia_nervosa"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["mental_illness_type_bipolar_disorder"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["mental_illness_type_bulimia_nervosa"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["mental_illness_type_depression"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["mental_illness_type_schizophrenia"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["mental_illness_type_substance_abuse"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["migraine"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["pku"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["sibo"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["skin_condition"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["thyroid"] != "Diagnosed by a medical professional (doctor, physician assistant)") &
        (df_clean["antibiotic_history"] == "I have not taken antibiotics in the past year.")]

print('Healthy controls:',len(df_healthy))


# Add column indicating level of healthy
df_disease["healthy"] = "no"
df_healthy["healthy"] = "yes"

# Combine results in the final data frame
df_output = pd.concat([df_disease, df_healthy], ignore_index=True)

#Rename sample_name to sample_name_id
df_output.rename(columns={'sample_name': 'sample_name_id'}, inplace=True)

# Save result
df_output.to_csv(cfg.clean_file, sep="\t", index=False)
logger.info("Combined file saved with %d samples", len(df_output))




