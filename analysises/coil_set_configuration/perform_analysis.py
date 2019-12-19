# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Over layer that performs the analysis."""

from analysises.coil_set_configuration.scripts.coil_set_positions import optimize_coils_positions

from simulation.elements.coils import Coil


if __name__ == "__main__":
    user_input = {"coil_type": Coil,
                  "n": 1,
                  "max_distance": 0.01}

    optimize_coils_positions(**user_input)

