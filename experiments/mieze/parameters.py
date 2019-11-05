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
SQUARE_COIL_POSITION_1ST = 0
SQUARE_COIL_POSITION_2ND = 0
