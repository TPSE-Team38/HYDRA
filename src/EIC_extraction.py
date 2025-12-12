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

def gaus(x, a, x0, sigma,c):
    return c + a * np.exp(-(x - x0)**2 / (2* sigma**2))

#not being used
def fit_curve(y):

    x=np.arange(1,len(y)+1)
    mu=np.mean(x)
    sigma=np.sqrt(5)

    initial_guess = [np.max(y), mu, sigma]

    params, error = curve_fit(gaus, x, y,p0=initial_guess)
    fit=gaus(x,*params)

    return fit

#not being used
#finds the peaks of the smoothed graph
#finds the smallest value between first and last peak: min
#then sets the mid point of interpolation to double of: biggest peak - min
#interpolates between first peak, mid point, last peak quadratically
def dip_detect_correct(y):
    """
    finds the peaks of the smoothed graph
    finds the smallest value between first and last peak: min
    then sets the mid point of interpolation to double of: biggest peak - min
    interpolates between first peak, mid point, last peak quadratically
    y: smoothed values of the y-axis
    Returns graph with dip replaced with interpolated values"""

    #get middle point by adding the highest peak to the difference between the highest peek and the lowest point
    peak, _ = sig.find_peaks(y, height=np.mean(y))
    max_peak=max(y[peak])
    minPoint = max_peak + (max_peak - min(y[peak[0]:peak[-1]]))

    pointsX=np.array([peak[0] - 1,
                      peak[-1] - ((peak[-1] - peak[0]) / 2),
                      peak[-1] + 1])
    pointsY=np.array([y[peak[0] - 1]
                         , minPoint
                         ,y[peak[-1] + 1]])
    #interpolate using first and last peak in x-axis and middle-point
    predict = interpolate.interp1d(pointsX,pointsY,kind='quadratic')

    #connect the interpolated values to the tails on either side of the left and right peak
    interpolated = np.concatenate((y[:(peak[0] - 1)], predict(np.arange(1,len(y))[peak[0] - 1:peak[len(peak) - 1] + 1]),
                                           y[peak[len(peak) - 1] + 1:]))

    return interpolated

def mask_region(y,left,right):
    return np.concatenate((y[:left],[np.nan]*(right-left),y[right:]))

def new_masking(y:np.ndarray[float],x:np.ndarray[float]):
    points=np.arange(1,len(y)-1)
    rng=np.random.default_rng()
    mid=rng.choice(points)
    np.delete(points,mid)
    leftP,rightP=mid-1,mid+1
    highest_r2=-1
    highest_r2_point=(-1,-1)
    highest_midpoint=-1
    i=0
    while len(points)>0:
        print(i)
        i+=1
        while leftP>0 and rightP<len(y):
            print(i)
            masked_y=mask_region(y,leftP,rightP)
            mask=~np.isnan(masked_y)
            fitted_y,sigma=different_approach_gaus(masked_y,x,len(y)/2)
            r_2=r2_score(masked_y[mask],fitted_y[mask])
            if r_2>highest_r2:
                highest_r2=r_2
                highest_r2_point=(leftP,rightP)
                highest_midpoint=mid
        mid=rng.choice(points)
        np.delete(points,mid)

    return mask_region(y,*highest_r2_point),mid

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
    boundL, boundR = 0.94, 0.95
    while not done:
        for (i, x) in enumerate(right_side):
            if boundL * y[rightP] <= x <= boundR * y[rightP] :
                actual_right = rightP+i
                break
        if actual_right == -1:
            boundL -= 0.0001
        else:
            done = True
    done = False
    boundL, boundR = 0.94, 0.95
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

#not being used
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

    return new_gauss_from_jonathan(x,*params),params[2]

def different_approach_gaus(y:np.ndarray[float],x:np.ndarray[float],xc):
    """
    fits the curve of the tails created by ion suppression removal
    """
    mask=~np.isnan(y)
    x_fit=x[mask]
    y_fit=y[mask]
    p0=[(y_fit.max() - y_fit.min())**2,xc,np.log1p(y_fit.max())*y_fit.std(),np.median(y_fit,axis=0)]
    params,_=curve_fit(gaus,x_fit,y_fit,maxfev=10000)

    return gaus(x,*params),params[2]

#not being used
def gauss_constant_fit(y:np.ndarray[float]):
    mod_gaus = lmfit.models.GaussianModel(nan_policy='propagate')
    mod_const=lmfit.models.ConstantModel(nan_policy='propagate')
    # model=mod_gaus+mod_const
    mask=~np.isnan(y)
    xdat=np.arange(1,len(y)+1)
    # pars = mod.make_params(c=np.mean(y_work),
    #                        center=xdat.mean(),
    #                        sigma=xdat.std(),
    #                        amplitude=xdat.std() * np.ptp(y_work)
    #                        )
    # pars=mod_gaus.guess(y[mask],x=xdat[mask])
    mod=mod_gaus+mod_const
    params=mod_gaus.guess(y[mask],x=xdat[mask])
    out=mod.fit(y[mask],x=xdat[mask],maxfev=10000)
    return out.eval(x=xdat)
    print(params)
    return mod.eval(params,x=xdat)
    # out = mod.fit(x=xdat[mask],params=params)
    # return []

#not being used
def gauss_from_jonathan_lmfit(y:np.ndarray[float],y_original:np.ndarray[float]):
    mod = lmfit.models.Model(new_gauss_from_jonathan,nan_policy='propagate')
    # mod = lmfit.models.GaussianModel(nan_policy='omit')
    y_work=y[~np.isnan(y)]
    xdat=np.arange(1,len(y)+1)
    pars=mod.make_params()
    # pars = mod.make_params(c=np.mean(y_work),
    #                        center=xdat.mean(),
    #                        sigma=xdat.std(),
    #                        amplitude=xdat.std() * np.ptp(y_work)
    #                        )
    out = mod.fit(y, pars, x=xdat)
    plt.figure(5, figsize=(8, 8))
    out.plot_fit()


def new_gauss_from_jonathan(y:np.ndarray[float],y_0,xc,w,A):
    return y_0 +((A/(w*np.sqrt(np.pi/2)))*np.exp((-2)*(((y-xc)/w)**2)))

#not being used
def gaussian_fit_with_removed_dip(y,y_original,x_c):
    obj=BaselineRemoval(y_original)

    y_mask=~np.isnan(y)
    xdat=np.arange(1,len(y)+1)
    y_0=y_original[obj.ZhangFit().argmin()]
    w=np.sqrt(2)*np.std(y_original)
    A=xdat.std() * np.ptp(y_original)

    initial_p=[np.float64(y_0),np.float64(x_c),w,A]

    params,_=curve_fit(new_gauss_from_jonathan,xdat[y_mask],y[y_mask],p0=initial_p)

    return new_gauss_from_jonathan(xdat,*params)

def get_z_vals(charge_state,charge_state_range):
    offset=np.floor(charge_state_range / 2)
    if charge_state_range % 2 == 0:
        return np.arange(charge_state - offset,
                           charge_state + offset + 1)
    else:
        return np.arange(charge_state - offset,
                           charge_state + offset + 2)

def mz_to_mz(original_mz,original_charge_state,new_charge_state):
    #mzOld*zOld / zNew
    return original_mz*original_charge_state/new_charge_state

def diffusion_coefficient(capillary_radius,standard_deviation : float,t_R:float):
    return((capillary_radius**2)*t_R)/(24*(standard_deviation**2))

def hydrodynamic_radius(temp:float,viscosity:float,diffusion_coefficient:float):
    from scipy.constants import k as boltzmann_c
    return (boltzmann_c*temp)/(6*np.pi*viscosity*diffusion_coefficient)

def recalculate(peaks,y,x,params):
    masked_y=np.concatenate((y[:peaks[0][0]],[np.nan]*(peaks[1][0]-peaks[0][0]),y[peaks[1][0]:]))
    mask=~np.isnan(masked_y)
    fitted,sigma=different_approach_gaus(masked_y,x,x[-1]/2)
    r2 = r2_score(masked_y[mask], fitted[mask])
    t_R = np.argmax(x) + 1
    D = diffusion_coefficient(float(params[2]), sigma, t_R)
    R_h = hydrodynamic_radius(float(params[0]), float(params[1]), D)
    return masked_y,fitted,r2,t_R,D,R_h

class result():
    def __init__(self,y,x,params,fig,ax):
        self.peaks=[]
        self.y=y
        self.x=x
        self.changed=False
        self.params=params
        self.fig,self.ax=fig,ax
        self.cid=self.fig.canvas.mpl_connect("button_press_event",self.on_click)
        self.recalculated_fit = None
        # self.recalculated_mask = None
        # self.recalculated_scatter = None

    def on_click(self,event):
        if not event.inaxes or len(self.peaks)>2:
            return
        c, v = int(event.xdata), event.ydata
        self.peaks.append((c, v))
        self.ax.plot(c, v, 'ro')
        self.fig.canvas.draw()
        if len(self.peaks) == 2:
            if self.changed:
                self.recalculated_fit.remove()
                # self.recalculated_mask.remove()
                # self.recalculated_scatter.remove()
            if self.peaks[0]>self.peaks[1]:
                temp=self.peaks[0]
                self.peaks[0]=self.peaks[1]
                self.peaks[1]=temp
                del temp
            self.changed=True
            masked_y,fitted_y,r2,t_R,D,R_h=recalculate(self.peaks, self.y, self.x, self.params)
            self.recalculated_fit=self.ax.plot(self.x, fitted_y,"--",label=f"recalculated_fit with r_2 score of {r2}\n t_R of {t_R}\n and diffusion coefficient of {D}\n and R_h of {R_h}")[0]
            # self.recalculated_mask=self.ax.plot(self.x, masked_y,"--",label="recalculated_mask")[0]
            # self.recalculated_scatter=self.ax.scatter(self.x, self.y,label="original")
            self.ax.set_ylim(min(*fitted_y,*self.y), max(*fitted_y,*self.y))
            self.fig.canvas.draw()
            self.ax.legend()
            self.peaks=[]
            plt.show()

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

    path = pathlib.Path(args.path)

    if not path.exists() or path.suffix.lower() != ".ms1":
        print("ERROR: The input file must exist and must be .ms1 format.")
        return

    # Convert region strings to floats
    regions = [(float(a), float(b),int(c),int(d),str(n)) for a, b, c, d,n in args.region]
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

        for z in z_vals:
            final_intensities_arr.append(get_final_eic_intensities(spectra, mz_to_mz(protein_mz,charge_state,z), protein_sampling_range))

        #summation of intensities of mz(s)
        final_intensities=[0]*len(final_intensities_arr[0])
        for final_intensity in final_intensities_arr:
            for i in range(len(final_intensity)):
                final_intensities[i]+=final_intensity[i]

        seconds = np.arange(1, len(final_intensities) + 1)
        ax.scatter(seconds, final_intensities, label=f"EIC of Protein {[float(mz_to_mz(protein_mz,charge_state,z)) for z in z_vals ]} with range {protein_sampling_range}")
        ax.set_xlim(0, seconds[-1])
        ax.set_ylim(min(final_intensities), max(final_intensities))
        #smoothing
        smoothed_intensities=sig.savgol_filter(final_intensities, int(args.smooth[0]), int(args.smooth[1]))
        # ax.plot(seconds,smoothed_intensities, label=f"smoothed EIC intensity of Protein {protein_mz}")

        if args.fit:
            print("Fitting EIC")
            #masking
            removed_dip,xc_guess=get_peaks(final_intensities)
            print(removed_dip)
            removed_dip,xc_guess=new_masking(final_intensities,seconds)
            mask=~np.isnan(removed_dip)
            print(removed_dip)
            ax.plot(seconds,removed_dip,label="EIC after Masking")

            #fitting and r2 score

            removed_dip_fitted,sigma = different_approach_gaus_jonathan(removed_dip, seconds, xc_guess)
            r2=r2_score(removed_dip[mask],removed_dip_fitted[mask])
            if r2<0.75:
                print(f"using normal gaus for EIC {protein_mz}")
                removed_dip_fitted,sigma= different_approach_gaus(removed_dip, seconds, xc_guess)
                r2=r2_score(removed_dip[mask],removed_dip_fitted[mask])
            t_R=np.argmax(removed_dip_fitted)+1
            D=diffusion_coefficient(float(args.parameters[2]),sigma,t_R)
            R_h=hydrodynamic_radius(float(args.parameters[0]),float(args.parameters[1]),D)
            ax.set_ylim(min(*removed_dip_fitted,*final_intensities), max(*removed_dip_fitted,*final_intensities))
            ax.plot(seconds,removed_dip_fitted,'--',label=f"EIC with R^2 value of {r2} \n and R_h of {R_h}  \n D: {D} \n sigma: {sigma}")

        ax.set_xlabel("Seconds")
        ax.set_ylabel("Total intensity")
        ax.set_title(f"Extracted Ion Chromatograms (EIC) of {protein_name} with {protein_mz}m/z; sampling range:{protein_sampling_range}")
        ax.legend()
        plt.grid(True)
        # ax.tight_layout()
        results.append(result(final_intensities,seconds,args.parameters,fig,ax))

    plt.show()


if __name__ == "__main__":
    main()
