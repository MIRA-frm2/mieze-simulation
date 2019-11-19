import numpy as np
from utils.physics_constants import unit


# Variables are stored as tuples with numerical value and units
I_sf1 = 1.6 * unit.A
I_hsf1 = 1.6 * unit.A
lambda_n = 4.3 / unit.angstrom

I_real_coil = 10 * unit.A

startpoint = 0.001 * unit.m
endpoint = 0.10 * unit.m  # Positions.get_position_coilA()
absolute_x_position = np.linspace(startpoint.magnitude, endpoint.magnitude, num=100)
step = (absolute_x_position[1] - absolute_x_position[0]) * unit.m

polarizer_position = 0 * unit.m
HelmholtzSpinFlipper_position_HSF1 = 0 * unit.m

# SpinFlipper positions
SpinFlipper_position1 = 0 * unit.m
SpinFlipper_position2 = 0 * unit.m

# MIEZE coil set values
L_IN = 0.086 * unit.m
L_OUT = 0.05 * unit.m
N_IN = 168 * unit.m
N_OUT = 48 * unit.m
R_IN = 0.177 / 2. * unit.m
R_OUT = 0.13 * unit.m

DISTANCE_2ND_COIL = 0.073 * unit.m
DISTANCE_3RD_COIL = 0.187 * unit.m
DISTANCE_4TH_COIL = 0.260 * unit.m

# Square coil Parameters
# ToDo: What are these values?
SQUARE_COIL_POSITION_1ST = 0.1 * unit.m
SQUARE_COIL_POSITION_2ND = 1 * unit.m

R_HSF = 5.38e-2 * unit.m  # radius of helmholtz coils at the spin flippers; achieved by fitting real magnetic field
WIDTH_CBOX = 2 * (50 + 86) * 1e-3 * unit.m  # width of one group of coils
L1 = 0.53 * unit.m  # coil group A to B
L2 = 2.58 * unit.m  # coil group B to detector
Ls = L2 - 0.62 * unit.m  # sample to detector
POLARISER_HSF1 = 10e-2 * unit.m  # distance between polariser and first coil of hsf1

POLARISATOR = 0.0 * unit.m
HSF1 = POLARISATOR + POLARISER_HSF1 + R_HSF / 2.0
SF1 = HSF1
COIL_A = HSF1 + R_HSF / 2.0
COIL_B = COIL_A + L1
HSF2 = COIL_B + R_HSF / 2.0
SF2 = HSF2
