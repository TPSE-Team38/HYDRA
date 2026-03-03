from .EIC_extraction import load_ms1, get_final_eic_intensities
from .eic import extract_eic
from .Fitting_and_masking import get_peaks, gaussian_fit,different_approach_gaus_jonathan
from .Calculations import *
from .models import EICResult
from sklearn.metrics import r2_score

def run_analysis(spectra,config):

    seconds, final_intensities = extract_eic(
        spectra,
        config.protein_mz,
        config.mz_window,
        config.charge_state,
        config.charge_range
    )

    removed_dip, xc_guess = get_peaks(final_intensities)
    mask = ~np.isnan(removed_dip)
    removed_dip_fitted, sigma = gaussian_fit(removed_dip, seconds, xc_guess)
    r2 = r2_score(removed_dip[mask], removed_dip_fitted[mask])

    #refitting with different gaussian if r2 bad
    if r2<0.9:
        temp_fit,temp_sigma=different_approach_gaus_jonathan(removed_dip, seconds, xc_guess)
        temp_r2=r2_score(removed_dip[mask], temp_fit[mask])
        if temp_r2>r2:
            r2=temp_r2
            removed_dip_fitted=temp_fit
            sigma=temp_sigma


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
        p
        )
