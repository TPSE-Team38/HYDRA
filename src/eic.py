import numpy as np
from .Calculations import get_z_vals, mz_to_mz

def extract_eic(spectra, protein_mz, mz_window, charge_state, charge_range):
    z_vals = get_z_vals(charge_state, charge_range)

    all_eics = []

    for z in z_vals:
        mz_target = mz_to_mz(protein_mz, charge_state, z)
        intensities = []

        for spectrum in spectra:
            mz = spectrum["m/z array"]
            inten = spectrum["intensity array"]

            mask = (mz >= mz_target - mz_window) & (mz <= mz_target + mz_window)
            intensities.append(inten[mask].sum())

        all_eics.append(intensities)

    summed = np.sum(all_eics, axis=0)
    x = np.arange(1, len(summed) + 1)

    return x, summed
