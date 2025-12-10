import scipy.optimize
from fontTools.subset import intersect
from numpy.ma.core import left_shift
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
import os
from BaselineRemoval import BaselineRemoval
from scipy.constants import k as boltzmann_c



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

def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / ( 2*sigma**2))



def get_peaks(y):
    #try getting peaks beginning from left then beginning from right (sigmoidal fit)

    maxDiff=0
    for i in range(0,len(y)-2):
        maxDtemp=max(abs(y[i]-y[i+1]),abs(y[i]-y[i+2]))
        if maxDtemp>maxDiff or maxDtemp>maxDiff:
            maxDiff=maxDtemp
    '''
    maxDiff=0
    for i in range(0,len(y)-1):
        maxDtemp=abs(y[i]-y[i+1])
        if maxDtemp>maxDiff or maxDtemp>maxDiff:
            maxDiff=maxDtemp
    '''
    leftP=0
    for i,x in enumerate(y):
        if x>y[leftP]:
            leftP=i
        elif abs(x-y[leftP])>maxDiff:
            break
    rightP=0
    for i in reversed(range(1,len(y))):
        if y[i]>y[rightP]:
            rightP=i
        elif abs(y[i]-y[rightP])>maxDiff:
            break
    # y_new=np.concatenate((y[:leftP],[np.nan]*len(y[leftP:rightP]),y[rightP:]))
    actual_right = -1
    actual_left = -1
    right_side = y[rightP:]
    left_side = y[:leftP+1]
    done = False
    boundL, boundR = 0.94, 0.96
    while not done:
        for (i, x) in enumerate(right_side):
            if boundL * y[rightP] <= x <= boundR * y[rightP]:
                actual_right = rightP+i
                break
        if actual_right == -1:
            boundL -= 0.0001
        else:
            done = True
    done = False
    boundL, boundR = 0.94, 0.96
    while not done:
        for (i, x) in enumerate(left_side):
            if boundL * left_side[leftP] <= x <= boundR * left_side[leftP]:
                actual_left = i
                break
        if actual_left == -1:
            boundL -= 0.0001
        else:
            done = True
    y_new=np.concatenate((y[:actual_left+1],[np.nan]*len(y[actual_left+1:actual_right]),y[actual_right:]))
    return y_new,leftP+((rightP-leftP)/2)

def get_2_peaks(y:np.ndarray[float],y_nonsmoothed)->np.ndarray[float]:
    """
    Returns y-values with the found noise created by ion suppression removed(replaced by np.nan)
    """

    y_smoothed=sig.savgol_filter(y_nonsmoothed,window_length=11,polyorder=2,mode='mirror')
    # plt.figure()
    # plt.plot(np.arange(1,len(y_smoothed)+1),y_smoothed)
    # plt.savefig("fig.png")
    peaks,_=sig.find_peaks(y_smoothed,height=y_smoothed.mean())
    x_c=peaks[0]+np.argmin(y[peaks[0]:peaks[-1]])

    right_peak=np.argmax(y_nonsmoothed[x_c+1:])
    left_peak=np.argmax(y_nonsmoothed[:x_c])
    actual_right=-1
    actual_left=-1
    right_side=y_nonsmoothed[right_peak+x_c:]
    left_side=y_nonsmoothed[:left_peak]
    done=False
    boundL,boundR=0.94,0.96
    while not done:
        for (i,x) in enumerate(right_side):
            if boundL*right_side[right_peak]<= x <=boundR*right_side[right_peak]:
                actual_right=x_c+i
                break
        if actual_right==-1:
            boundL-=0.01
        else:
            done=True
    done=False
    boundL,boundR=0.94,0.96
    while not done:
        for (i,x) in enumerate(left_side):
            if boundL*left_side[left_peak-1]<= x <=boundR*left_side[left_peak-1]:
                actual_left=i
                break
        if actual_left==-1:
            boundL-=0.01
        else:
            done=True
    maxPeak=max(y_nonsmoothed[actual_left],y_nonsmoothed[actual_right])
    if (maxPeak-y_nonsmoothed[x_c])/maxPeak>0.3:
        new_y=np.concatenate((y_nonsmoothed[:actual_left],[np.nan]*len(y[actual_left:actual_right]),y_nonsmoothed[actual_right:]))
        # print((max(y_nonsmoothed[actual_left],y_nonsmoothed[actual_right])-y_nonsmoothed[x_c])/max(y_nonsmoothed[actual_left],y_nonsmoothed[actual_right]))
    else:
        new_y=y_nonsmoothed.copy()
        fit=different_approach_gaus(new_y,np.arange(1,len(new_y)+1))
        inter=np.intersect1d(y,fit)
        # print(inter)

    return new_y,x_c

def different_approach_gaus_jonathan(y:np.ndarray[float],x:np.ndarray[float],xc):
    """
    fits the curve of the tails created by ion suppression removal
    """
    mask=~np.isnan(y)
    x_fit=x[mask]
    y_fit=y[mask]
    # print(np.median(y_fit,axis=0))
    p0 = [np.median(y_fit,axis=0), xc, np.log1p(y_fit.max())*y_fit.std(), (y_fit.max() - y_fit.min())**2]
    params,_=curve_fit(new_gauss_from_jonathan,x_fit,y_fit,p0=p0,maxfev=10000)
    sigma = abs(params[2] / 2)
    return new_gauss_from_jonathan(x,*params),sigma

def different_approach_gaus(y:np.ndarray[float],x:np.ndarray[float],xc):
    """
    fits the curve of the tails created by ion suppression removal
    """
    mask=~np.isnan(y)
    x_fit=x[mask]
    y_fit=y[mask]
    p0=[(y_fit.max() - y_fit.min())**2,xc,np.log1p(y_fit.max())*y_fit.std()]
    params,_=curve_fit(gaus,x_fit,y_fit,maxfev=10000)
    sigma=abs(params[2])
    return gaus(x,*params),sigma



def new_gauss_from_jonathan(y:np.ndarray[float],y_0,xc,w,A):
    return y_0 +((A/(w*np.sqrt(np.pi/2)))*np.exp((-2)*(((y-xc)/w)**2)))

def get_z_vals(charge_state,charge_state_range):
    offset=np.floor(charge_state_range / 2)
    if charge_state_range % 2 == 0:
        return np.arange(charge_state - offset,
                           charge_state + offset + 1)
    else:
        return np.arange(charge_state - offset,
                           charge_state + offset + 2)

def mz_to_mz(original_mz,original_charge_state,new_charge_state):
    return original_mz*original_charge_state/new_charge_state

def diffusion_coefficient(capillary_radius,standard_deviation : float,t_R:float):
    return(((capillary_radius**2)*t_R)/(24*(standard_deviation**2)))

def hydrodynamic_radius(temp:float,viscosity:float,diffusion_coefficient:float):
    from scipy.constants import k as boltzmann_c
    return (boltzmann_c*temp)/(6*np.pi*viscosity*diffusion_coefficient)

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
        nargs=4,
        action="append",
        metavar=('MIN_MZ', 'Protein_Smapling_Range','Z','Z_range'),
        required=True,
        help="-R <m/z value> <region size> , example: --region 1689 4 7 2 (can be repeated)"
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=22,

        help=" <temperature in Kelvin>  (default: 22 C°)"
    )
    parser.add_argument(
        "--viscosity",
        type=float,
        default=0.9544e-3,

        help="<viscosity> (default: 0.0009544 kg m^-1 s^-1)"
    )

    parser.add_argument(
        "--parameters","-P",
        nargs=2,
        action="store",
        metavar=('capillary_radius','capillary_length'),
        required=True,
        help="-P <capillary_radius> <capillary_length>\n example -P 0.000127 1.1455 "
    )
    args = parser.parse_args()

    path = pathlib.Path(args.path)

    if not path.exists() or path.suffix.lower() != ".ms1":
        print("ERROR: The input file must exist and must be .ms1 format.")
        return

    # Convert region strings to floats
    regions = [(float(a), float(b),int(c),int(d)) for a, b, c, d in args.region]
    # Load spectra
    print(f"Loading MS1 file: {path}")
    spectra = load_ms1(path)
    print(f"Loaded {len(spectra)} spectra.")
    plt.style.use('ggplot')

    # Extract EICs for each region
    for (i,(protein_mz, protein_sampling_range,charge_state,charge_state_range)) in enumerate(regions):

        plt.figure(figsize=(10, 6))
        final_intensities_arr=[]
        z_vals=get_z_vals(charge_state,charge_state_range)

        for z in z_vals:
            final_intensities_arr.append(get_final_eic_intensities(spectra, mz_to_mz(protein_mz,charge_state,z), protein_sampling_range))

        #summation of intensities of mz(s)
        final_intensities=[0]*len(final_intensities_arr[0])
        for final_intensity in final_intensities_arr:
            for i in range(len(final_intensity)):
                final_intensities[i]+=final_intensity[i]

        seconds = np.arange(1, len(final_intensities) + 1)

        plt.scatter(seconds, final_intensities, label=f"EIC of Protein {protein_mz} with range {protein_sampling_range}")

        #smoothing
        smoothed_intensities=sig.savgol_filter(final_intensities, int(args.smooth[0]), int(args.smooth[1]))
        plt.plot(seconds,smoothed_intensities, label=f"smoothed EIC intensity of Protein {protein_mz}")

        if args.fit:
            print("Fitting EIC")
            #masking
            removed_dip,xc_guess=get_peaks(final_intensities)
            mask=~np.isnan(removed_dip)
            plt.plot(seconds,removed_dip,label="EIC after Masking",color='yellow')

            #fitting and r2 score
            removed_dip_fitted,sigma = different_approach_gaus_jonathan(removed_dip, seconds, xc_guess)
            r2=r2_score(removed_dip[mask],removed_dip_fitted[mask])
            if r2<0.5:
                removed_dip_fitted,sigma = different_approach_gaus(removed_dip, seconds, xc_guess)
                r2=r2_score(removed_dip[mask],removed_dip_fitted[mask])
            t_R=np.argmax(removed_dip_fitted)+1
            #sigma=np.std(removed_dip_fitted)
            D=diffusion_coefficient(float(args.parameters[0]),sigma,t_R)
            R_h=hydrodynamic_radius(float(args.temperature+273.15),float(args.viscosity),D)
            plt.plot(seconds,removed_dip_fitted,'--',label=f"EIC with R^2 value of {r2} \n \n and R_h of {R_h}")

            print("diffusion_coefficient =" ,D)
            print("standard deviation of the Gaussian fit =",sigma)
            print("retention time=",t_R,"seconds")
            print("R_h =", R_h)
        plt.xlabel("Seconds")
        plt.ylabel("Total intensity")
        plt.title(f"Extracted Ion Chromatograms (EIC) of {protein_mz}m/z of range {protein_sampling_range}")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
# --smooth 10 3 --fit --region 1689 4 7 2 -P 0.000127 1.145
