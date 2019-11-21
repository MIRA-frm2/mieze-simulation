import numpy as np


# Variables are stored as tuples with numerical value and units
I_sf1 = 1.6  # [A]
I_hsf1 = 1.6  # [A]
lambda_n = 4.3  # 1/[Angstrom]

I_real_coil = 10  # [A]

startpoint = 0.001  # [A]
endpoint = 0.10  # [A]  # Positions.get_position_coilA()
npoints = 100
absolute_x_position = np.linspace(startpoint, endpoint, num=npoints)
step = (absolute_x_position[1] - absolute_x_position[0])  # [A]

polarizer_position = 0  # [A]
HelmholtzSpinFlipper_position_HSF1 = 0  # [A]

# SpinFlipper positions
SpinFlipper_position1 = 0  # [A]
SpinFlipper_position2 = 0  # [A]

# MIEZE coil set values
L_IN = 0.086  # [A]
L_OUT = 0.05   # [A]
N_IN = 168   # [A]
N_OUT = 48   # [A]
R_IN = 0.177 / 2.   # [A]
R_OUT = 0.13   # [A]

DISTANCE_2ND_COIL = 0.073   # [A]
DISTANCE_3RD_COIL = 0.187   # [A]
DISTANCE_4TH_COIL = 0.260   # [A]

# Square coil Parameters
# ToDo: What are these values?
SQUARE_COIL_POSITION_1ST = 0.1   # [A]
SQUARE_COIL_POSITION_2ND = 1   # [A]

R_HSF = 5.38e-2   # [A]  # radius of helmholtz coils at the spin flippers; achieved by fitting real magnetic field
WIDTH_CBOX = 2 * (50 + 86) * 1e-3   # [A]  # width of one group of coils
L1 = 0.53   # [A]  # coil group A to B
L2 = 2.58   # [A]  # coil group B to detector
Ls = L2 - 0.62   # [A]  # sample to detector
POLARISER_HSF1 = 10e-2   # [A]  # distance between polariser and first coil of hsf1

POLARISATOR = 0.0   # [A]
HSF1 = POLARISATOR + POLARISER_HSF1 + R_HSF / 2.0
SF1 = HSF1
COIL_A = HSF1 + R_HSF / 2.0
COIL_B = COIL_A + L1
HSF2 = COIL_B + R_HSF / 2.0
SF2 = HSF2
