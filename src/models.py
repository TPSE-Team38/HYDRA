from dataclasses import dataclass
import numpy as np

@dataclass
class AnalysisConfig:

    ms1_path: str
    protein_name: str
    protein_mz: float
    mz_window: float
    charge_state: int
    charge_range: int

    temperature: float
    viscosity: float
    capillary_radius: float
    capillary_length: float
    flow_rate: float


@dataclass
class EICResult:
    # change indicator
    changed:bool
    # metadata
    protein_name: str
    protein_mz: float
    mz_window: float
    charge_state: int
    charge_range: int

    seconds: np.ndarray
    final_intensities: np.ndarray
    removed_dip: np.ndarray
    removed_dip_fitted: np.ndarray


    r2: float
    tR: float
    sigma: float
    D: float
    Rh: float
    t: float
    p: float

    original_removed_dip:np.ndarray
    original_removed_dip_fitted: np.ndarray

    original_r2: float
    original_tR: float
    original_sigma: float
    original_D: float
    original_Rh: float
    original_t: float
    original_p: float
