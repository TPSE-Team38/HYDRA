import numpy as np

def diffusion_coefficient(capillary_radius,standard_deviation : float,t_R:float):
    return(((capillary_radius*(10**-6))**2)*t_R)/(24*(standard_deviation**2))

def hydrodynamic_radius(temperature:float,viscosity:float,diffusion_coefficient:float):
    from scipy.constants import k as boltzmann_c
    return (boltzmann_c*(temperature+273.15))/(6*np.pi*viscosity*diffusion_coefficient)

def peclet(R_h,temperature,viscosity,capillary_radius,flow_rate):
    from scipy.constants import k as boltzmann_c
    return 6 * viscosity * (flow_rate*(10**-9)/60) * R_h / (boltzmann_c * (temperature+273.15) * capillary_radius*(10**-6))

def tau(temperature,capillary_length,viscosity,flow_rate,estimated_R_h):
    from scipy.constants import k as boltzmann_c
    return boltzmann_c*(temperature+273.15)*(capillary_length* 10**-2)/(6*viscosity*(flow_rate*(10**-9)/60)*(estimated_R_h))

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
    return c + a * np.exp(-((x - x0)**2) / (2* sigma**2))

def new_gauss_from_jonathan(y:np.ndarray[float],y_0,xc,w,A):
    return y_0 +((A/(w*np.sqrt(np.pi/2)))*np.exp((-2)*(((y-xc)/w)**2)))

''''''