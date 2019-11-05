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


def main():
    experiment = Mieze(coil_type=Coil, sample_distance=1.5, coil_distance=0.53, detector_distance=2.503, increment=0.05)
    experiment.create_setup(current=5)

    # rho = 0.02
    # start = -0.25
    # end = experiment.sample_distance

    b_field_kwargs = {'coordinate_system': 'cylindrical', 'plane': 'xy', 'rho': 0.25,
                      'start': -0.25, 'end': 1.5}
    experiment.calculate_b_field(**b_field_kwargs)

    experiment.plot_field_1d_abs()
    # experiment.plot_field_2d_abs()
    # experiment.plot_field_1d_vec()
    # experiment.plot_2d_vectormap()


if __name__ == "__main__":
    main()
