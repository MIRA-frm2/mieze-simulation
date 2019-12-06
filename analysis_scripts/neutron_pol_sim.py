# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the flow of the particles through the magnetic field."""

# import numpy as np

from particles.beam import NeutronBeam

from experiments.mieze.parameters import I_hsf1, step_x, startpoint, beamend, beamsize, step_x, velocity

from utils.physics_constants import earth_field


# def b_function(vec):
#     x = vec[0]
#     y = vec[1]
#     z = vec[2]
#
#     flipper = HelmholtzSpinFlipper()
#     flipper.create_element(position=0.1, current=I_hsf1)
#
#     return flipper.hsf(x, y, z) + Polariser.B_field(x, y, z) + earth_field
#

def main():
    simulation = NeutronBeam(beamsize=beamsize,
                             incrementsize=step_x,
                             number_of_neutrons=1,
                             velocity=velocity,
                             totalflightlength=beamend)

    grid_size = {'x_start': startpoint, 'x_end': beamend, 'x_step': step_x,
                 'y_start': -1.0, 'y_end': 1.0, 'z_start': -1.0, 'z_end': 1.0,
                 'yz_step':  (1.0 - -1.0) / 20}

    simulation.initialize_computational_space(**grid_size)

    simulation.create_neutrons()
    simulation.reset_pol()

    print(simulation.get_pol())
    print(simulation.get_neutron_position())

    # simulation.create_b_map(b_function, (0.001, 0.001))
    simulation.simulate_neutron_trajectory()
    print(simulation.get_pol())
    print(simulation.get_neutron_position())

    # 
    # B_polariser = polariser.pol(absolute_y_position)
    # B_extra = RealCoil.Bz(z=absolute_y_position, coil_mid_pos=0.05, rho=0.01, l=0.1, N=100, I=I_real_coil, R=0.05)
    # 
    # Bx = B_polariser + Spin_Flipper.sf('sf1', I_sf1, absolute_y_position) #G
    # By = Helmholtz_Spin_Flipper.hsf('hsf1', I_hsf1, absolute_y_position) + B_extra #G
    # Long_Coil.hsf('long', 20.0, 0.05, absolute_y_position)+ Single_Coil.hsf('extra', 10.0, 0.02, absolute_y_position) + Single_Coil.hsf('extra', 10.0, 0.03, absolute_y_position)+ Single_Coil.hsf('extra', 10.0, 0.04, absolute_y_position) ++ Single_Coil.hsf('extra', 10.0, 0.05, absolute_y_position)+ Single_Coil.hsf('extra', 10.0, 0.06, absolute_y_position)
    # B = np.sqrt(Bx**2+By**2)
    # theta = np.degrees(np.arctan(By/Bx))
    # dtheta_dy = gradient(theta, step)


if __name__ == "__main__":
    main()
