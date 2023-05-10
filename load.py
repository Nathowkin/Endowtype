# Taken from Meow

import numpy as np
import pandas
import pathlib
from sklearn import preprocessing
from paths import get_tcga_paths
from typing import Callable, List, Union

#
# Constants
#
COHORTS = {
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
    "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
    "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
    "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS",
    "UVM"
}

TASKS = {
    "OS", "PFS"
}

FEATURES = {
    "moco",
}

NORMALIZATIONS = {
    "raw", "norm", "rpkm", "tpm",
}

#
# Load functions
#
def load_metadata(cohort: Union[str, List[str]]) -> pandas.DataFrame:
    """Load metadata for a TCGA cohort.

    Parameters
    ----------
    cohort : str, List[str]
        Name of the TCGA cohort, or list of cohorts
    Returns
    -------
    pd.DataFrame
    """
    data_paths = get_tcga_paths().data_paths

    if isinstance(cohort, List):
        return pandas.concat([load_metadata(cc) for cc in cohort], axis='index').dropna(axis=1)

    # Check input cohort
    assert cohort in COHORTS, f"Please select a cohort in {COHORTS}"

    # Load dataframe
    if isinstance(data_paths["path_to_labels"], Callable):
        path = data_paths["path_to_labels"](cohort)
    else:
        path = data_paths["path_to_labels"]
    df_metadata = pandas.read_excel(path)

    # Select cohort
    df_metadata = df_metadata[(df_metadata.type == cohort)]
    df_metadata = df_metadata.rename(
        columns={
            "bcr_patient_barcode": "patient_id",
            "OS": "OS_event",
            "OS.time": "OS_time",
            "PFI": "PFS_event",
            "PFI.time": "PFS_time",
        }
    )
    df_metadata.set_index("patient_id", inplace=True)
    return df_metadata


def load_labels(cohort: Union[str, List[str]], task: str = "OS") -> pandas.DataFrame:
    """Load survival labels for a TCGA cohort.

    Parameters
    ----------
    cohort : str, List[str]
        Name of the TCGA cohort, or list of cohorts
    task : str
        Survival task (OS or PFS).
    Returns
    -------
    pd.DataFrame
    """
    data_paths = get_tcga_paths().data_paths
    if isinstance(data_paths["path_to_labels"], Callable):
        path = data_paths["path_to_labels"](cohort)
    else:
        path = data_paths["path_to_labels"]

    if isinstance(cohort, List):
        return pandas.concat([load_labels(cc, task=task) for cc in cohort], axis='index').dropna(axis=1)

    # Check input cohort and task
    assert cohort in COHORTS, f"Please select a cohort in {COHORTS}"
    assert task in TASKS, f"Please select a task in {TASKS}"

    # Load dataframe
    df_labels = pandas.read_excel(path)

    # Select cohort
    df_labels = df_labels[(df_labels.type == cohort)]
    df_labels = df_labels.rename(
        columns={
            "bcr_patient_barcode": "patient_id",
            "OS": "OS_event",
            "OS.time": "OS_time",
            "PFI": "PFS_event",
            "PFI.time": "PFS_time",
        }
    )

    if task == "OS":
        df_labels.dropna(subset=["OS_event", "OS_time"], inplace=True)
        df_labels["label"] = df_labels.apply(
            lambda row: row.OS_time if row.OS_event == 1 else -row.OS_time, axis=1
        )
    elif task == "PFS":
        df_labels.dropna(subset=["PFS_event", "PFS_time"], inplace=True)
        df_labels["label"] = df_labels.apply(
            lambda row: row.PFS_time if row.PFS_event == 1 else -row.PFS_time, axis=1
        )

    df_labels = df_labels[["patient_id", "label"]]
    df_labels.set_index("patient_id", inplace=True)
    return df_labels


def load_clinical(cohort: Union[str, List[str]]) -> pandas.DataFrame:
    """Load clinical variables for a TCGA cohort.

    Parameters
    ----------
    cohort : str, List[str]
        Name of the TCGA cohort, or list of cohorts
    Returns
    -------
    pd.DataFrame
    """
    data_paths = get_tcga_paths().data_paths

    if isinstance(cohort, List):
        return pandas.concat([load_clinical(cc) for cc in cohort], axis='index').dropna(axis=1)

    # Check input cohort
    assert cohort in COHORTS, f"Please select a cohort in {COHORTS}"

    # Load dataframe
    if isinstance(data_paths["path_to_clinical"], Callable):
        path = data_paths["path_to_clinical"](cohort)
    else:
        path = pathlib.Path(data_paths["path_to_clinical"]) / f"TCGA-{cohort}_clinical.tsv.gz"

    df_clin = pandas.read_csv(path, sep="\t")

    # Rename columns
    df_clin = df_clin.rename(columns={"patient": "patient_id"})

    # Drop columns
    to_keep = [
        "patient_id",
        "age_at_diagnosis",
        "gender",
        "ethnicity",
        "race",
        "ajcc_pathologic_stage",
        "high_grade",
        "prior_malignancy",
        "definition",
        "synchronous_malignancy",
        "Purety_CCF",
        "treatments.treatment_or_therapy",
        "treatments.treatment_type",
    ]
    to_keep = list(set(to_keep).intersection(set(df_clin.columns)))
    df_clin = df_clin[to_keep]

    # Remove columns with only NaN values
    for col in df_clin.columns:
        if df_clin[col].isna().sum() == len(df_clin):
            df_clin = df_clin.drop(col, 1)

    # Fill NaN values
    for col in df_clin.columns:
        try:
            df_clin[col] = df_clin[col].fillna(df_clin[col].median())
        except:
            df_clin[col] = df_clin[col].fillna("Unknown")

    # Encode categorical variables
    categorical_vars = [
        "gender",
        "ethnicity",
        "race",
        "high_grade",
        "prior_malignancy",
        "synchronous_malignancy",
        "treatments.treatment_type",
        "treatments.treatment_or_therapy",
        "definition",
        "ajcc_pathologic_stage",
        "ajcc_pathologic_t",
        "ajcc_pathologic_n",
        "ajcc_pathologic_m",
        "paper_MSI_status",
        "prior_malignancy",
        "prior_treatment",
        "synchronous_malignancy",
        "morphology",
        "paper_methylation_subtype",
        "paper_histological_type",
        "paper_history_of_colon_polyps",
        "paper_lymphatic_invasion_present",
        "paper_lymphnode_pathologic_spread",
        "paper_primary_tumor_pathologic_spread",
    ]
    categorical_vars = list(
        set(categorical_vars).intersection(set(df_clin.columns)))
    for col in categorical_vars:
        df_clin[col] = preprocessing.LabelEncoder().fit_transform(df_clin[col]).astype(int)

    # Scale continuous variables
    continuous_vars = [
        "age_at_diagnosis",
        "age_at_index",
        "initial_weight",
        "Purety_CCF",
    ]
    continuous_vars = list(
        set(continuous_vars).intersection(set(df_clin.columns)))
    for col in continuous_vars:
        df_clin[col] = (df_clin[col] - np.min(df_clin[col])) / (
            np.max(df_clin[col]) - np.min(df_clin[col])
        )

    # df_clin.columns == ["patient"] + clinical_variables
    df_clin.set_index("patient_id", inplace=True)
    return df_clin


def load_project_meta(project : str):
    data_paths = get_tcga_paths().data_paths

    if isinstance(data_paths["path_to_rnaseq"], Callable):
        path = pathlib.Path(data_paths["path_to_rnaseq"](project)) / "raw" / "metadata.tsv.gz"
    else:
        path = f'{data_paths["path_to_rnaseq"]}/{project}/Data/metadata.tsv.gz'

    meta = pandas.read_csv(path, sep='\t', index_col='external_id')

    assert(meta['tcga.tcga_barcode'].isna().sum() == 0)
    #assert(meta['tcga.cgc_file_last_modified_date'].isna().sum() == 0)

    meta['tcga.cgc_file_last_modified_date'] = pandas.to_datetime(meta['tcga.cgc_file_last_modified_date'], utc=True)
    meta['patient_id'] = meta['tcga.tcga_barcode'].map(lambda x : x[:12])
    meta['sample_id'] = meta['tcga.tcga_barcode'].map(lambda x : x[:15])
    meta['barcode'] = meta['tcga.tcga_barcode']
    meta['is_normal_tissue'] = meta['sample_id'].map(lambda x : 10 <= int(x[-2:]) <= 29)

    # Select latest sample per barcode (discard multiple uploads)
    meta = meta.sort_values('tcga.cgc_file_last_modified_date').groupby(['tcga.tcga_barcode']).tail(1)
    return meta


def load_rnaseq(cohort: Union[str, List[str]], normalization: str = "norm", index_by_barcode: bool = False) -> pandas.DataFrame:
    """Load RNASeq data for a TCGA cohort.

    Parameters
    ----------
    cohort : str, List[str]
        Name of the TCGA cohort, or list of cohorts
    normalization : {"raw", "norm", "rpkm", "tpm"}
        Normalization type

    Returns
    -------
    pandas.DataFrame
    """
    data_paths = get_tcga_paths().data_paths

    if isinstance(data_paths["path_to_rnaseq"], Callable):
        path = pathlib.Path(data_paths["path_to_rnaseq"](cohort)) / "processed" / f"Counts_{normalization}_float32.parquet"
    else:
        path = f'{data_paths["path_to_rnaseq"]}/{cohort}/Data/Counts_{normalization}.tsv.gz'

    if isinstance(cohort, List):
        return pandas.concat([load_rnaseq(cc, normalization=normalization, index_by_barcode=index_by_barcode) for cc in cohort], axis='index').dropna(axis=1)

    # Check input cohort and normalization method
    assert cohort in COHORTS, f"Please select a cohort in {COHORTS}"
    assert normalization in NORMALIZATIONS, (
        f"Please select a normalization method in {NORMALIZATIONS}"
    )

    meta = load_project_meta(cohort)
    if str(path).endswith(".tsv.gz"):
        counts = pandas.read_csv(path, engine='c', sep='\t', index_col='Hugo').transpose()
    else:
        counts = pandas.read_parquet(path).set_index("Hugo").transpose()

    counts = counts.loc[meta.query("is_normal_tissue == False").index]
    if index_by_barcode:
        counts = counts.rename(index=meta['barcode'].to_dict())
        counts.index.rename("barcode", inplace=True)
    else:
        counts = counts.rename(index=meta['patient_id'].to_dict())
        counts.index.rename("patient_id", inplace=True)
    return counts


def load_histo(cohort: Union[str, List[str]], features: str = "moco") -> pandas.DataFrame:
    """Load histo data for a TCGA cohort.
    Parameters
    ----------
    cohort : str, List[str]
        Name of the TCGA cohort, or list of cohorts
    features: {"imagenet", "moco"}
        Which histo features to use.

    Returns
    -------
    Dataframe containing histo paths
    """
    tcga_paths = get_tcga_paths()

    if isinstance(cohort, List):
        return pandas.concat([load_histo(cc, features=features) for cc in cohort], axis='index').dropna(axis=1)

    # Check input cohort and features
    assert cohort in COHORTS, f"Please select a cohort in {COHORTS}"
    assert features in FEATURES, f"Please select features in {FEATURES}"

    if tcga_paths.hostname != "TDE_ABSTRA":
        if features == "moco":
            path = tcga_paths.moco_path_primary(cohort)
            if not path.exists():
                path = tcga_paths.moco_path_secondary(cohort)
            if not path.exists():
                path = tcga_paths.moco_path_tertiary(cohort)
            assert path.exists()
        else:
            raise ValueError(f"{features} path not defined.")
    else:
        path = pathlib.Path(tcga_paths.data_paths["path_to_histo"](cohort, features))

    hpaths = list(path.glob("*/features.npy"))

    df_histo = pandas.DataFrame({"histo_path": hpaths})
    df_histo["slide_id"] = df_histo.histo_path.apply(lambda x: x.parent.name[:-4])
    df_histo["patient_id"] = df_histo.histo_path.apply(lambda x: x.parent.name[:12])

    df_histo.set_index("patient_id", inplace=True)
    return df_histo


if __name__ == "__main__":
    cohort = "LUAD"
    df_meta = load_metadata(cohort)
    df_labels = load_labels(cohort)
    df_clin = load_clinical(cohort)
    df_rnaseq = load_rnaseq(cohort)
    df_histo = load_histo(cohort, "moco")
    print(df_histo)
