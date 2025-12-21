from dataclasses import dataclass
import numpy as np

@dataclass
class AnalysisConfig:
    ms1_path: str
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
