import numpy as np


# Variables are stored as tuples with numerical value and units
I_sf1 = 1.6  # [A]
I_hsf1 = 1.6  # [A]
lambda_n = 4.3  # 1/[Angstrom]

I_real_coil = 10  # [A]

startpoint = 0.001  # [A]
endpoint = 10  # [A]  # Positions.get_position_coilA()
npoints = 100
absolute_x_position = np.linspace(startpoint, endpoint, num=npoints)
step = (absolute_x_position[1] - absolute_x_position[0])  # [A]

# MIEZE coil set values
L_IN = 0.086  # [m]
L_OUT = 0.05   # [m]
N_IN = 168   # [m]
N_OUT = 48   # [m]
R_IN = 0.177 / 2.   # [m]
R_OUT = 0.13   # [m]

DISTANCE_2ND_COIL = 0.073   # [m]
DISTANCE_3RD_COIL = 0.187   # [m]
DISTANCE_4TH_COIL = 0.260   # [m]

# Square coil Parameters
SQUARE_COIL_POSITION_1ST = 0.1   # [A]
SQUARE_COIL_POSITION_2ND = 1   # [A]

# Rectangular coils
RECTANGULAR_COIL_LENGTH = 1 * 1e-2  # [m]
RECTANGULAR_COIL_WIDTH = 10 * 1e-2  # [m]
RECTANGULAR_COIL_HEIGHT = 13 * 1e-2  # [m]
WIRE_D = 1 * 1e-2  # [m]
WINDINGS = 1

# Coil group
R_HSF = 5.38e-2   # [A]  # radius of helmholtz coils at the spin flippers; achieved by fitting real magnetic field
WIDTH_CBOX = 2 * (50 + 86) * 1e-3   # [A]  # width of one group of coils
L1 = 0.53   # [A]  # coil group A to B
L2 = 2.58   # [A]  # coil group B to detector
Ls = L2 - 0.62   # [A]  # sample to detector
POLARISER_HSF1 = 10e-2   # [A]  # distance between polariser and first coil of hsf1

POLARISATOR = 0.0   # [A]
HelmholtzSpinFlipper_position_HSF1 = 2  # POLARISATOR + POLARISER_HSF1 + R_HSF / 2.0
SpinFlipper_position1 = HelmholtzSpinFlipper_position_HSF1
COIL_A = HelmholtzSpinFlipper_position_HSF1 + R_HSF / 2.0
COIL_B = COIL_A + L1
HelmholtzSpinFlipper_position_HSF2 = COIL_B + R_HSF / 2.0
SpinFlipper_position2 = HelmholtzSpinFlipper_position_HSF2
