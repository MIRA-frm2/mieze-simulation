# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the magnetic field for a given experiment."""

from simulation.experiments.mieze.main_mieze import main_mieze
from simulation.experiments.mieze.parameters import startpoint, beamend, step_x


if __name__ == "__main__":

    main_mieze()
