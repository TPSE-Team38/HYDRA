import numpy as np
from scipy.stats import boltzmann


def diffusion_coefficient(capillary_radius,standard_deviation : float,t_R:float):
    if standard_deviation <= 0:
        return np.nan
    return((capillary_radius**2)*t_R)/(24*(standard_deviation**2))

def hydrodynamic_radius(temp:float,viscosity:float,diffusion_coefficient:float):
    from scipy.constants import k as boltzmann_c
    return (boltzmann_c*temp)/(6*np.pi*viscosity*diffusion_coefficient)

def peclet(R_h,temperature,viscosity,capillary_radius,flow_rate):
    if boltzmann == 0 or capillary_radius <= 0 or temperature <= 0:
        return np.nan
    from scipy.constants import k as boltzmann_c
    return 6 * viscosity * (flow_rate*(10**-9) / 60) * R_h / (boltzmann_c * (temperature) * capillary_radius)

def tau(T,L,viscosity,Q,R_h):
    if viscosity <= 0 or Q <= 0 or R_h <= 0 :
        return np.nan
    from scipy.constants import k as k_b
    return (k_b*T*L)/(6*viscosity*Q*R_h)

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

''''''
def gaus(x, a, x0, sigma,c):
    if sigma == 0:
        return np.nan
    return c + a * np.exp(-((x - x0)**2) / (2* sigma**2))

def new_gauss_from_jonathan(y:np.ndarray[float],y_0,xc,w,A):
    return y_0 +((A/(w*np.sqrt(np.pi/2)))*np.exp((-2)*(((y-xc)/w)**2)))

''''''