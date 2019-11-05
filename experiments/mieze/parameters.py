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

