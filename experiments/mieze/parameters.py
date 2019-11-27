import numpy as np


# Variables are stored as tuples with numerical value and units
I_sf1 = 1.6  # [A]
I_hsf1 = 1.6  # [A]
lambda_n = 4.3  # 1/[Angstrom]

I_real_coil = 10  # [A]

startpoint = 0.001  # [m]
endpoint = 10  # [m]  # Positions.get_position_coilA()
npoints = 100
absolute_x_position = np.linspace(startpoint, endpoint, num=npoints)
step = (absolute_x_position[1] - absolute_x_position[0])  # [A]

# MIEZE coil set values
LENGTH_COIL_INNER = 0.086  # [m]
LENGTH_COIL_OUTER = 0.050   # [m]
N_WINDINGS_COIL_INNER = 168   # [m]
N_WINDINGS_COIL_OUTER = 48   # [m]
RADIUS_COIL_INNER_EFFECTIVE = 0.177 / 2.   # [m]
RADIUS_COIL_OUTER_EFFECTIVE = 0.13   # [m]

RADIUS_COIL_INNER_MIN = 50 * 1e-3  # [m]
RADIUS_COIL_INNER_MAX = 126 * 1e-3  # [m]

DISTANCE_2ND_COIL = 73 * 1e-3   # [m]
DISTANCE_3RD_COIL = 187 * 1e-3   # [m]
DISTANCE_4TH_COIL = 260 * 1e-3   # [m]

WIRE_D = 5 * 1e-3  # [m]
COIL_SET_CURRENT = 5  # [A]

# Square coil Parameters
SQUARE_COIL_POSITION_1ST = 0.1   # [A]
SQUARE_COIL_POSITION_2ND = 1   # [A]

# Rectangular coils
RECTANGULAR_COIL_LENGTH = 1 * 1e-2  # [m]
RECTANGULAR_COIL_WIDTH = 10 * 1e-2  # [m]
RECTANGULAR_COIL_HEIGHT = 13 * 1e-2  # [m]
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
