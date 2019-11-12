# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the magnetic field for a given experiment."""

from experiments.mieze.main import Mieze
from elements.coils import Coil

from experiments.mieze.parameters import L1, L2

from utils.helper_functions import save_data_to_file


def main():
    experiment = Mieze(coil_type=Coil, sample_distance=1.5, coil_distance=L1, detector_distance=L2-L1,
                       increment=0.05)
    experiment.create_setup(current=5)

    # rho = 0.02
    # start = -0.25
    # end = experiment.sample_distance

    grid_size = {'x_start': -0.2, 'x_end': 1.2, 'y_start': -0.0, 'y_end': 0.0, 'z_start': -0.0, 'z_end': 0.0}
    experiment.initialize_computational_space(**grid_size)
    experiment.calculate_b_field()
    save_data_to_file(experiment.b, file_name='./data/data')


if __name__ == "__main__":
    main()
