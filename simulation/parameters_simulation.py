import numpy as np

from simulation.beamline.beamline_properties import neutron_speed

startpoint = 0.000  # [m]
beamend = 0.25  # [m]

npoints = 250

absolute_x_position = np.linspace(startpoint, beamend, num=npoints)

step_x = (absolute_x_position[1] - absolute_x_position[0])

default_beam_grid = {'x_start': startpoint, 'x_end': beamend, 'x_step': step_x,
                     'y_start': -0.0, 'y_end': 0.0, 'z_start': -0.0, 'z_end': 0.0, 'yz_step': 0.1}

total_simulation_time = beamend / neutron_speed
