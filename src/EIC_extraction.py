import random

import scipy.optimize
from fontTools.subset import intersect
from numpy.ma.core import left_shift, masked
from pyteomics import ms1
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import scipy.stats as stats
import scipy.signal as sig
import scipy.optimize as opt
import argparse
import pathlib
from scipy.optimize import curve_fit,least_squares
import lmfit
from scipy.signal import find_peaks
from sklearn.metrics import r2_score
from .Fitting_and_masking import *
from .plotting import ResultPlot
from .Calculations import *
def load_ms1(path):
    """Read all MS1 spectra from a file into a list."""
    spectra = []
    with open(path, "r") as f:
        for spect in ms1.read(f):
            spectra.append(spect)
    return spectra

def is_in_region(spectrum, min_mz, max_mz):
    """ return true or false , if m/z inside region. select only m/z , that inside the region [min_mz, max_mz]"""
    mz_array = spectrum['m/z array']
    return (mz_array>= min_mz) & (mz_array <= max_mz)

def get_intensities_of_region(spectrum, min_mz, max_mz):
    """Return all intensities of only m/z , that inside the region [min_mz, max_mz]."""
    intensity_array = spectrum["intensity array"]
    mask = is_in_region(spectrum, min_mz, max_mz)
    return intensity_array[mask]

def get_final_eic_intensities(spectra, protein_mz, protein_sampling_range)->np.ndarray:
    """Return final intensities array of our EIC by summing intensity in the m/z window for each spectrum."""

    #  array of sum of intensities in the region per spectrum
    final_intensities = []

    for spect in spectra:

        # array of intensities in one region
        intensities_per_region = get_intensities_of_region(spect, protein_mz-protein_sampling_range, protein_mz+protein_sampling_range)

        # sum of intensities in each region
        sum_of_intensities_per_region = np.sum(intensities_per_region)

        # final EIC Intensities
        final_intensities.append(sum_of_intensities_per_region)

    return np.array(final_intensities)

def get_all_intensity(spectra,protein_mz,protein_sampling_range,original_charge_state,charge_list):
    """sums intensities for each second including only m/z extracted using charge_list"""
    final_intensities_arr=[]
    for z in charge_list:
        final_intensities_arr.append(
            get_final_eic_intensities(spectra, mz_to_mz(protein_mz, original_charge_state, z), protein_sampling_range))

    # summation of intensities of mz(s)
    final_intensities = [0] * len(final_intensities_arr[0])
    for final_intensity in final_intensities_arr:
        for i in range(len(final_intensity)):
            final_intensities[i] += final_intensity[i]

    return final_intensities

def recalculate(peaks,y,x,params):
    masked_y=np.concatenate((y[:peaks[0][0]],[np.nan]*(peaks[1][0]-peaks[0][0]),y[peaks[1][0]:]))
    mask=~np.isnan(masked_y)
    fitted,sigma=gaussian_fit(masked_y,x,x[-1]/2)
    r2 = r2_score(masked_y[mask], fitted[mask])
    t_R = np.argmax(y) + 1
    D = diffusion_coefficient(float(params[2]), sigma, t_R)
    R_h = hydrodynamic_radius(float(params[0]), float(params[1]), D)
    return masked_y,fitted,r2,t_R,D,R_h

results=[]

def main():

    parser = argparse.ArgumentParser(
        description="Extract EICs from an MS1 file for one or many m/z regions."
    )

    parser.add_argument(
        "path",
        type=str,
        help="Path to the .ms1 file"
    )

    parser.add_argument(
        "--smooth","-S",
        nargs=2,
        action="store",
        metavar=('Window Length', 'Poly order'),
        required=True,
        help="Savitzky-Golay Filter Parameters, ex. --smooth 10 3,! 1st param>=2nd param"
    )

    parser.add_argument(
        "--fit","-F",
        required=False,
        action="store_true",
        default=False,
        help="set this flag if a gaussian fitting should be done"
    )

    parser.add_argument(
        "--region","-R",
        nargs=5,
        action="append",
        metavar=('Protein_MZ', 'Protein_Smapling_Range','Z','Z_range','name'),
        required=True,
        help="-R <m/z value> <region size>, <charge_state>, <other_charge_states>(region),<name> , example: --region 1689 4 7 2 FKBP12 (can be repeated)"
    )

    parser.add_argument(
        "--parameters","-P",
        nargs=5,
        action="store",
        metavar=('temperature','viscosity','capillary_radius','capillary_length','flow_rate'),
        required=True,
        help="-P <temperature> <viscosity> <capillary_radius> <capillary_length> <flow_rate> \n example -P 22 0.0009544 127 114.5 20"
    )
    args = parser.parse_args()
    params=np.asarray([float(x) for x in args.parameters],dtype=float)

    path = pathlib.Path(args.path)

    if not path.exists() or path.suffix.lower() != ".ms1":
        print("ERROR: The input file must exist and must be .ms1 format.")
        return

    # Convert region strings to floats
    regions = [(float(a), float(b),int(c),int(d),str(n)) for a, b, c, d,n in args.region]
    # print(params)
    # Load spectra
    print(f"Loading MS1 file: {path}")
    spectra = load_ms1(path)
    print(f"Loaded {len(spectra)} spectra.")
    plt.style.use('ggplot')
    # Extract EICs for each region
    for (i,(protein_mz, protein_sampling_range,charge_state,charge_state_range,protein_name)) in enumerate(regions):

        fig,ax=plt.subplots()
        final_intensities_arr=[]
        z_vals=get_z_vals(charge_state,charge_state_range)
        #
        # for z in z_vals:
        #     final_intensities_arr.append(get_final_eic_intensities(spectra, mz_to_mz(protein_mz,charge_state,z), protein_sampling_range))
        #
        # #summation of intensities of mz(s)
        # final_intensities=[0]*len(final_intensities_arr[0])
        # for final_intensity in final_intensities_arr:
        #     for i in range(len(final_intensity)):
        #         final_intensities[i]+=final_intensity[i]
        final_intensities=get_all_intensity(spectra,protein_mz,charge_state,protein_sampling_range,z_vals)
        seconds = np.arange(1, len(final_intensities) + 1)
        ax.set_xlim(0, seconds[-1])
        ax.set_ylim(min(final_intensities), max(final_intensities))
        #smoothing
        smoothed_intensities=sig.savgol_filter(final_intensities, int(args.smooth[0]), int(args.smooth[1]))
        # ax.plot(seconds,smoothed_intensities, label=f"smoothed EIC intensity of Protein {protein_mz}")

        if args.fit:
            print("Fitting EIC")
            #masking
            removed_dip,xc_guess=get_peaks(final_intensities)
            # removed_dip,xc_guess=new_masking(final_intensities,seconds,params)
            print("done masking")
            mask=~np.isnan(removed_dip)
            ax.scatter(seconds,removed_dip,label="EIC after Masking")
            ax.scatter(seconds, [v if v not in removed_dip else np.nan for v in final_intensities ],label=f"EIC of Protein {[float(mz_to_mz(protein_mz, charge_state, z)) for z in z_vals]} with range {protein_sampling_range}")

            #fitting and r2 score

            # removed_dip_fitted,sigma = different_approach_gaus_jonathan(removed_dip, seconds, xc_guess)
            # r2=r2_score(removed_dip[mask],removed_dip_fitted[mask])
            # if r2<0.75:
            # print(f"using normal gaus for EIC {protein_mz}")
            # removed_dip_fitted,sigma= different_approach_gaus_jonathan(removed_dip, seconds, xc_guess)
            removed_dip_fitted,sigma= gaussian_fit(removed_dip, seconds, xc_guess)
            r2=r2_score(removed_dip[mask],removed_dip_fitted[mask])
            t_R=np.argmax(removed_dip_fitted)+1
            D=diffusion_coefficient(float(params[2]),sigma,t_R)
            R_h=hydrodynamic_radius(float(params[0]),float(params[1]),D)
            t=tau(params[0],params[3],params[1],params[4],R_h)
            p=peclet(R_h,params[0],params[1],params[2],params[4])
            ax.set_ylim(min(*removed_dip_fitted,*final_intensities), max(*removed_dip_fitted,*final_intensities))
            ax.plot(seconds,removed_dip_fitted,'--',label=f"EIC with R^2 value of {r2} \n and R_h of {R_h}  \n D: {D} \n sigma: {sigma} \n tau: {t} \n peclet: {p}")

        ax.set_xlabel("Seconds")
        ax.set_ylabel("Total intensity")
        ax.set_title(f"Extracted Ion Chromatograms (EIC) of {protein_name} with {protein_mz}m/z; sampling range:{protein_sampling_range}")
        ax.legend()
        plt.grid(True)
        # ax.tight_layout()
        results.append(ResultPlot(final_intensities, seconds, params, fig, ax, recalculate))

    plt.show()


if __name__ == "__main__":
    main()
