# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Physically relevant constants."""

import numpy as np
from pint import UnitRegistry, set_application_registry

unit = UnitRegistry()
set_application_registry(unit)

Q_ = unit.Quantity

pi = np.pi
MU_0 = 4e-7 * np.pi * unit.H / unit.m  # N/A, vacuum permeability
earth_field = np.array((0, 210, -436)) * 1e-3
