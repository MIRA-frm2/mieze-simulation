import numpy as np

I_sf1 = 1.6  # A
I_hsf1 = 1.6  # A
lambda_n = 4.3  # Ångström^{-1}

I_real_coil = 10  # A

endpoint = 0.10  # Positions.get_position_coilA() #m
absolute_x_position = np.linspace(0.001, endpoint, num=1000)
step = 0.001  # absolute_y_position[1] - absolute_y_position[0] #m

polarizer_position = 0
HelmholtzSpinFlipper_position_HSF1 = 0

# SpinFlipper positions
SpinFlipper_position1 = 0
SpinFlipper_position2 = 0

# MIEZE coil set values
L_IN = 0.086
L_OUT = 0.05
N_IN = 168
N_OUT = 48
R_IN = 0.177 / 2.
R_OUT = 0.13

DISTANCE_2ND_COIL = 0.073
DISTANCE_3RD_COIL = 0.187
DISTANCE_4TH_COIL = 0.260

# Square coil Parameters
# ToDo: What are these values
SQUARE_COIL_POSITION_1ST = 0.1
SQUARE_COIL_POSITION_2ND = 1

R_HSF = 5.38e-2  # [m], radius of helmholtz coils at the spin flippers; achieved by fitting real magnetic field
WIDTH_CBOX = 2 * (50 + 86) * 1e-3  # [m], width of one group of coils
L1 = 0.53  # [m], coil group A to B
L2 = 2.58  # [m], coil group B to detector
Ls = L2 - 0.62  # [m], sample to detector
POLARISER_HSF1 = 10e-2  # [m], distance between polariser and first coil of hsf1

POLARISATOR = 0.0
HSF1 = POLARISATOR + POLARISER_HSF1 + R_HSF / 2.0
SF1 = HSF1
COIL_A = HSF1 + R_HSF / 2.0
COIL_B = COIL_A + L1
HSF2 = COIL_B + R_HSF / 2.0
SF2 = HSF2