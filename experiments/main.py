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

from experiments.mieze.parameters import L1, L2, startpoint, endpoint, step

from utils.helper_functions import save_data_to_file


def main():
    experiment = Mieze(coil_type=Coil, coil_distance=L1, detector_distance=L2-L1,
                       increment=step,  sample_distance=1.5)
    experiment.create_setup(current=5)

    grid_size = {'x_start': startpoint, 'x_end': endpoint, 'y_start': -0.0, 'y_end': 0.0, 'z_start': -0.0, 'z_end': 0.0}
    experiment.initialize_computational_space(**grid_size)

    experiment.calculate_b_field()

    save_data_to_file(experiment.b, file_name='../data/data')


if __name__ == "__main__":
    main()
