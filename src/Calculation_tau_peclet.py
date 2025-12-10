import argparse
from scipy.constants import k as boltzmann_c

parser = argparse.ArgumentParser(
        description="Calculation of tau and Peclet based on theoretical Rh"
    )

parser.add_argument("--Rh", type=float,required=True,
                    help="Rh")

parser.add_argument("--temperature", type=float, required=True,
                    help="temperature in C°")

parser.add_argument("--viscosity", type=float, required=True,
                    help="viscosity")

parser.add_argument("--capillary_radius", type=float, required=True,
                    help="Capillary radius in meters")

parser.add_argument("--capillary_length", type=float, required=True,
                    help="Capillary length in meters")

parser.add_argument("--flow_rate", type=float, required=True,
                    help="Flow rate in µL/min")

def calculate_Tau_and_Peclet(estimated_R_h,temperature ,viscosity,capillary_radius,capillary_length,flow_rate):
    tau= boltzmann_c*(temperature+273.15)*capillary_length/(6*viscosity*(flow_rate*(10**-9)/60)*estimated_R_h)
    peclet=6*viscosity*(flow_rate*(10**-9)/60)*estimated_R_h/(boltzmann_c*(temperature+273.15)*capillary_radius)
    print(f"tau={tau} peclet={peclet}")

if __name__ == "__main__":
    args = parser.parse_args()

    calculate_Tau_and_Peclet(
        args.Rh,
        args.temperature,
        args.viscosity,
        args.capillary_radius,
        args.capillary_length,
        args.flow_rate
    )

# --Rh 2.216366203600859e-09 --temperature 22 --viscosity 0.9544e-3 --capillary_radius 0.000127 --capillary_length 1.145 --flow_rate 20


