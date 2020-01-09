# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the MIEZE setup."""

from simulation.experiments.setup import Setup
from simulation.experiments.mieze.parameters import (
    I_hsf1, HelmholtzSpinFlipper_position_HSF1, POLARISATOR, R_HSF, L1, L2, default_beam_grid, CoilSet_position,
    distance_between_HSF1_coilset)

from utils.helper_functions import save_data_to_file


from simulation.elements.coils import Coil
from simulation.elements.coil_set import CoilSet
from simulation.elements.helmholtz_pair import HelmholtzPair
from simulation.elements.polariser import Polariser
# from elements.spin_flipper import SpinFlipper


class Mieze(Setup):
    def __init__(self, coil_type, sample_distance, detector_distance,
                 spin_flipper_distance=HelmholtzSpinFlipper_position_HSF1,
                 coil_set_distance=distance_between_HSF1_coilset):

        super(Mieze, self).__init__()

        self.sample_distance = sample_distance
        self.detector_distance = detector_distance
        self.spin_flipper_distance = spin_flipper_distance
        self.coil_set_distance = coil_set_distance

        self.coil_type = coil_type

    def create_setup(self, current):
        self.create_element(element_class=Polariser,
                            position=(POLARISATOR, 0, 0))

        self.create_element(coil_type=Coil,
                            current=I_hsf1,
                            element_class=HelmholtzPair,
                            position=(self.spin_flipper_distance, 0, 0),
                            radius=R_HSF)

        # self.create_element(current=I_sf1,
        #                     element_class=SpinFlipper,
        #                     height=RECTANGULAR_COIL_HEIGHT,
        #                     length=RECTANGULAR_COIL_LENGTH,
        #                     position=(SpinFlipper_position1, 0, 0),
        #                     width=RECTANGULAR_COIL_WIDTH,
        #                     windings=WINDINGS,
        #                     wire_d=WIRE_D)

        self.create_element(element_class=CoilSet,
                            current=100,  # [A]
                            position=CoilSet_position + self.coil_set_distance)


def main_mieze(grid_size=default_beam_grid,
               coil_set_distance=distance_between_HSF1_coilset,
               spin_flipper_distance=HelmholtzSpinFlipper_position_HSF1,
               filename='data/data_magnetic_field'):

    experiment = Mieze(coil_type=Coil,
                       spin_flipper_distance=spin_flipper_distance,
                       coil_set_distance=coil_set_distance,
                       detector_distance=L2-L1,
                       sample_distance=1.5)

    experiment.create_setup(current=5)

    experiment.initialize_computational_space(**grid_size)

    experiment.calculate_b_field()

    save_data_to_file(experiment.b, file_name=filename)
