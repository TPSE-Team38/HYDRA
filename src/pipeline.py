from .EIC_extraction import load_ms1, get_final_eic_intensities
from .eic import extract_eic
from .Fitting_and_masking import get_peaks, gaussian_fit, gaussian_fit_no_mask
from .Calculations import *
from .models import EICResult, AnalysisConfig
from sklearn.metrics import r2_score

def run_analysis(spectra,config:AnalysisConfig):

    seconds, final_intensities = extract_eic(
        spectra,
        config.protein_mz,
        config.mz_window,
        config.charge_state,
        config.charge_range,
        config.measurement_start,
        config.measurement_end
    )

    if not config.disable_masking:
        removed_dip, xc_guess = get_peaks(final_intensities)
        mask = ~np.isnan(removed_dip)
        removed_dip_fitted, sigma = gaussian_fit(removed_dip, seconds, xc_guess)
        r2 = r2_score(removed_dip[mask], removed_dip_fitted[mask])
    else:
        removed_dip=final_intensities
        xc_guess=np.argmax(removed_dip)
        removed_dip_fitted,sigma=gaussian_fit_no_mask(removed_dip, seconds, xc_guess)
        r2=r2_score(removed_dip,removed_dip_fitted)

    tR = np.argmax(removed_dip_fitted) + 1
    D = diffusion_coefficient(config.capillary_radius, sigma, tR)
    Rh = hydrodynamic_radius(config.temperature, config.viscosity, D)
    t = tau(config.temperature, config.capillary_length,
            config.viscosity, config.flow_rate, Rh)
    p = peclet(Rh, config.temperature, config.viscosity,
               config.capillary_radius, config.flow_rate)

    return EICResult(
        False,
        config.protein_name,
        config.protein_mz,
        config.mz_window,
        config.charge_state,
        config.charge_range,
        seconds,
        final_intensities,
        removed_dip,
        removed_dip_fitted,
        r2,
        tR,
        sigma,
        D,
        Rh,
        t,
        p,
        removed_dip,
        removed_dip_fitted,
        r2,
        tR,
        sigma,
        D,
        Rh,
        t,
        p,
        config.measurement_start,
        config.measurement_end,
        config.disable_masking
        )
