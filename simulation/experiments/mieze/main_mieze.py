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
    I_hsf1, POLARISATOR, R_HSF, default_beam_grid, CoilSet_position,
    COIL_SET_PARAMETERS, ELEMENTS_POSITIONS, SPIN_FLIPPER_PARAMETERS)

from utils.helper_functions import save_data_to_file


from simulation.elements.coils import Coil
from simulation.elements.coil_set import CoilSet
from simulation.elements.helmholtz_pair import HelmholtzPair
from simulation.elements.polariser import Polariser
from simulation.elements.spin_flipper import SpinFlipper


class Mieze(Setup):

    def __init__(self, coil_type=Coil, **kwargs):

        super(Mieze, self).__init__()

        self.spin_flipper_distance = kwargs.get('spin_flipper_distance')
        self.coil_set_distance = kwargs.get('coil_set_distance')

        self.coil_type = coil_type

    def create_setup(self):
        """Create the elements for the setup."""

        self.create_element(element_class=Polariser,
                            position=(POLARISATOR, 0, 0))

        self.create_element(coil_type=Coil,
                            current=I_hsf1,
                            element_class=HelmholtzPair,
                            position=(self.spin_flipper_distance, 0, 0),
                            radius=R_HSF)

        self.create_element(current=SPIN_FLIPPER_PARAMETERS["I_sf1"],
                            element_class=SpinFlipper,
                            height=SPIN_FLIPPER_PARAMETERS["RECTANGULAR_COIL_HEIGHT"],
                            length=SPIN_FLIPPER_PARAMETERS["RECTANGULAR_COIL_LENGTH"],
                            position=(SPIN_FLIPPER_PARAMETERS["SpinFlipper_position1"], 0, 0),
                            width=SPIN_FLIPPER_PARAMETERS["RECTANGULAR_COIL_WIDTH"],
                            windings=SPIN_FLIPPER_PARAMETERS["WINDINGS"],
                            wire_d=SPIN_FLIPPER_PARAMETERS["WIRE_D"])

        self.create_element(element_class=CoilSet,
                            current=COIL_SET_PARAMETERS["CURRENT"],  # [A]
                            name='CoilSet',
                            position=CoilSet_position + self.coil_set_distance)


def main_mieze(grid_size=default_beam_grid,
               coil_set_distance=ELEMENTS_POSITIONS["coil_set_distance"],
               spin_flipper_distance=ELEMENTS_POSITIONS["spin_flipper_distance"],
               filename='data/data_magnetic_field'):

    # Initialize an object from the MIEZE class
    experiment = Mieze(coil_type=Coil,
                       spin_flipper_distance=spin_flipper_distance,
                       coil_set_distance=coil_set_distance)

    # Create the components of the beamline with their parameters
    experiment.create_setup()
    experiment.update_metadata()

    # Initialize the computational space (grid) and compute the magnetic field for it
    experiment.initialize_computational_space(**grid_size)
    experiment.calculate_b_field()

    # Compute the magnetic field for one point only
    # output = experiment.calculate_b_field(point=(0, 0, 0))
    # print(output)

    # Save the obtained data to a file
    experiment.save_data_to_file(filename=filename)

    # Plot results
    # experiment.set_plot_ticks(set_ticks=False)

    # experiment.plot_field_1d_scalar(component='x')
    # experiment.plot_field_1d_scalar(component='y')
    # experiment.plot_field_1d_scalar(component='z')

