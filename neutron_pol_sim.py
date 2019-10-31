#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 12:37:37 2018

@author: henrik
"""

import numpy as np

from elements.polariser import polariser
from elements.helmholtz_spin_flipper import HelmholtzSpinFlipper
from particles.neutron import Neutrons
#from extra_helmholtz import Extra_Helmholtz
#from single_coil import Single_Coil
#from long_coil import Long_Coil
from elements.coils import RealCoil

from experiments.mieze.parameters import I_hsf1, I_real_coil, step, endpoint

from utils.physics_constants import earth_field


def b_function(vec):
    x = vec[0]
    y = vec[1]
    z = vec[2]
    return HelmholtzSpinFlipper.hsf(x, y, z, current=I_hsf1, mid_pos=0.1) + polariser.B_field(x,y,z) + earth_field


def main():

    # Bx = lambda x,y,z: RealCoil.Bx(x, coil_mid_pos=0.05, rho=np.sqrt(z**2+y**2), l=0.1, N=100, I=I_real_coil, R=0.05) + RealCoil.Bx(x, coil_mid_pos=0.05, rho=np.sqrt(z**2+y**2), l=0.1, N=100, I=I_real_coil, R=0.05)
    # By = lambda x,y,z: RealCoil.pol(x) + RealCoil.Brho(x, coil_mid_pos=0.05, rho=y, l=0.1, N=100, I=I_real_coil, R=0.05)
    # Bz = lambda x,y,z: RealCoil.Brho(x, coil_mid_pos=0.05, rho=z, l=0.1, N=100, I=I_real_coil, R=0.05)

    simulation = Neutrons(10000, incrementsize=step, totalflightlength=endpoint)
    # simulation.reset_pol()

    print(simulation.get_pol())
    # simulation.create_B_map(B_function, (0.001, 0.001))
    simulation.simulate_neutrons()
    print(simulation.get_pol())

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
