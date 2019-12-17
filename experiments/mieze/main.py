# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the MIEZE setup."""

from experiments.setup import Setup
from experiments.mieze.parameters import (
    I_hsf1, I_sf1, HelmholtzSpinFlipper_position_HSF1, SQUARE_COIL_POSITION_1ST, POLARISATOR, R_HSF)

from elements.coils import Coil
from elements.coil_set import CoilSet
from elements.helmholtz_pair import HelmholtzPair
from elements.polariser import Polariser
# from elements.spin_flipper import SpinFlipper


class Mieze(Setup):
    def __init__(self, coil_type, sample_distance, coil_distance, detector_distance):

        super(Mieze, self).__init__()

        self.sample_distance = sample_distance
        self.detector_distance = detector_distance
        self.coil_distance = coil_distance

        self.coil_type = coil_type

    def create_setup(self, current):
        self.create_element(element_class=Polariser,
                            position=(POLARISATOR, 0, 0))

        self.create_element(coil_type=Coil,
                            current=I_hsf1,
                            element_class=HelmholtzPair,
                            position=(HelmholtzSpinFlipper_position_HSF1, 0, 0),
                            radius=R_HSF)

        # self.create_element(current=I_sf1,
        #                     element_class=RectangularCoil,
        #                     height=RECTANGULAR_COIL_HEIGHT,
        #                     length=RECTANGULAR_COIL_LENGTH,
        #                     position=(SpinFlipper_position1, 0, 0),
        #                     width=RECTANGULAR_COIL_WIDTH,
        #                     windings=WINDINGS,
        #                     wire_d=WIRE_D)

        self.create_element(element_class=CoilSet,
                            current=I_sf1,
                            position=SQUARE_COIL_POSITION_1ST)
