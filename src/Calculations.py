import numpy as np

from scipy.constants import k as boltzmann_c


def diffusion_coefficient(capillary_radius,standard_deviation : float,t_R:float):

    radius_in_meter = capillary_radius * (10 ** -6)

    #if standard_deviation <= 0:
      #  return np.nan
    return((radius_in_meter**2)*t_R)/(24*(standard_deviation**2))

def hydrodynamic_radius(temp:float,viscosity:float,diffusion_coefficient:float):
    T_in_kelvin = temp + 273.15
    return (boltzmann_c*T_in_kelvin)/(6*np.pi*viscosity*diffusion_coefficient)

def peclet(R_h,temperature,viscosity,capillary_radius,flow_rate):

    T_in_kelvin = temperature + 273.15
    radius_in_meter = capillary_radius * (10 ** -6)
    Q_in_meter_cube_per_sec = (flow_rate * (10 ** -9)) / 60

    ##if boltzmann_c == 0 or radius_in_meter <= 0 or T_in_kelvin <= 0:
      #  return np.nan
    return (6 * viscosity * Q_in_meter_cube_per_sec * R_h) / (boltzmann_c * T_in_kelvin * radius_in_meter)

def tau(T,L,viscosity,Q,R_h):
    #if viscosity <= 0 or Q <= 0 or R_h <= 0 :
     #   return np.nan

    T_in_kelvin = T + 273.15
    L_in_meter = L * (10**-2)
    Q_in_meter_cube_per_sec = (Q * (10**-9)) / 60

    return (boltzmann_c * T_in_kelvin * L_in_meter) / (6 * viscosity * Q_in_meter_cube_per_sec * R_h)

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