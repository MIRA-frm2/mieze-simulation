# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the MIEZE setup."""

from experiments.setup import Setup
from experiments.mieze.parameters import (DISTANCE_2ND_COIL, DISTANCE_3RD_COIL, DISTANCE_4TH_COIL,
                                          I_real_coil, I_hsf1, I_sf1,
                                          L_IN, L_OUT, N_IN, N_OUT, R_IN, R_OUT,
                                          SpinFlipper_position1, HelmholtzSpinFlipper_position_HSF1,
                                          SQUARE_COIL_POSITION_1ST, SQUARE_COIL_POSITION_2ND)

from elements.coils import SquareCoil, RealCoil
from elements.spin_flipper import SpinFlipper
from elements.helmholtz_spin_flipper import HelmholtzSpinFlipper


class Mieze(Setup):
    def __init__(self, coil_type, sample_distance, coil_distance, detector_distance, increment=0.001):

        super(Mieze, self).__init__(increment=increment)

        self.sample_distance = sample_distance
        self.detector_distance = detector_distance
        self.coil_distance = coil_distance

        self.coil_type = coil_type

    # def _create_coil_set(self, first_coil_pos=0, current=5):
    #     """Crate the magnetic coil sets."""
    #     self.create_element(element_class=self.coil_type, coil_mid_pos=first_coil_pos, length=L_OUT,
    #                         windings=N_OUT, current=-current, r=R_OUT)
    #     self.create_element(element_class=self.coil_type, coil_mid_pos=first_coil_pos + DISTANCE_2ND_COIL,
    #                         length=L_IN, windings=N_IN, current=current, r=R_IN)
    #     self.create_element(element_class=self.coil_type, coil_mid_pos=first_coil_pos + DISTANCE_3RD_COIL,
    #                         length=L_IN, windings=N_IN, current=current, r=R_IN)
    #     self.create_element(element_class=self.coil_type, coil_mid_pos=first_coil_pos + DISTANCE_4TH_COIL,
    #                         length=L_OUT, windings=N_OUT, current=-current, r=R_OUT)

    def set_plot_ticks(self, set_ticks=False):
        if set_ticks:
            self.x_ticks = [0, self.coil_distance, SQUARE_COIL_POSITION_1ST, SQUARE_COIL_POSITION_2ND]
            self.x_ticks_labels = ['1st Coil Set', '2nd Coil Set', '1st Square Coil ', '2nd Square Coil']

    def create_setup(self, current):

        current1 = current
        current2 = current * self.detector_distance / (self.detector_distance - self.coil_distance)

        self.create_element(element_class=RealCoil, position=(0.05, 0, 0), length=0.1, windings=100,
                            current=I_real_coil, r=0.05)
        self.create_element(element_class=SquareCoil, position=(SpinFlipper_position1, 0, 0), current=I_sf1)
        self.create_element(element_class=HelmholtzSpinFlipper, position=(HelmholtzSpinFlipper_position_HSF1, 0, 0), current=I_hsf1)

        # self.create_element(element_class=SquareCoil, coil_mid_pos=SQUARE_COIL_POSITION_1ST)
        #
        # self._create_coil_set(current=current1)
        # self._create_coil_set(current=current2, first_coil_pos=self.coil_distance)

        # self.create_element(element_class=SquareCoil, coil_mid_pos=SQUARE_COIL_POSITION_2ND)

        self.set_plot_ticks()
