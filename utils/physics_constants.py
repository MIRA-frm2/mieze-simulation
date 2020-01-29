# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Physically relevant constants."""

import numpy as np

pi = np.pi
MU_0 = 4e-7 * np.pi  # N/A, vacuum permeability
earth_field = np.array((0, 210, -436)) * 1e-3
factor_T_to_G = 1e4
