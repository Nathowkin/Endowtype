# Taken from Meow

"""Define paths to load TCGA data in different hosts."""

import os
import pathlib
from typing import Callable, Dict, NamedTuple, Union

from loguru import logger


class TCGAPaths(NamedTuple):
    data_paths: Dict[str, Union[str, Callable[[], str]]]
    hostname: str
    moco_path_primary: Callable[[str], str]
    moco_path_secondary: Callable[[str], str]
    moco_path_tertiary: Callable[[str], str]

POSSIBLE_HOSTNAMES =  {"OWKINLAB_GCP", "TARGET_ENGINE_GCP", "MIDDLE_EARTH", "TDE_ABSTRA"}


def get_tcga_paths():
    # Get hostname from environment variable
    
    hostname = 'TDE_ABSTRA'

    def _path_to_labels(cohort):
        del cohort
        return "/home/owkin/project/dataset/tcga/labels/1-s2.0-S0092867418302290-mmc1.xlsx"

    def _path_to_clinical(cohort):
        return f"../../data/TCGA_{cohort}/clinical/raw/TCGA-{cohort}_clinical.tsv.gz"

    def _path_to_rnaseq(cohort):
        return f"../../data/TCGA_{cohort}/omics"

    def _path_to_histo(cohort, features):
        return  f"../../data/TCGA_{cohort}/histology/processed/MoCoWideResNetCOAD"
    
    moco_path_primary = None
    moco_path_secondary = None
    moco_path_tertiary = None

    data_paths = {
        "path_to_labels": _path_to_labels,
        "path_to_clinical": _path_to_clinical,
        "path_to_rnaseq": _path_to_rnaseq,
        "path_to_histo": _path_to_histo,
    }

    return TCGAPaths(
        data_paths=data_paths,
        moco_path_primary=moco_path_primary,
        moco_path_secondary=moco_path_secondary,
        moco_path_tertiary=moco_path_tertiary,
        hostname=hostname,
    )