# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the MIEZE setup."""

from experiments.setup import Setup


class Mieze(Setup):
    def __init__(self, sample_distance, coils_dinstance, detector_distance, increment=0.001):

        super(Mieze, self).__init__(increment=increment)

        self.sample_distance = sample_distance
        self.detector_distance = detector_distance
        self.coil_distance = coils_dinstance

        self.coil_type = 'simple'

    def _create_coil_set(self, first_coil_pos=0, current=5):
        l_in = 0.086
        l_out = 0.05
        n_in = 168
        n_out = 48
        r_in = 0.177 / 2.
        r_out = 0.13

        self.create_coil(coil_type=self.coil_type, coil_mid_pos=first_coil_pos, length=l_out, windings=n_out, current=-current, r=r_out)
        self.create_coil(coil_type=self.coil_type, coil_mid_pos=first_coil_pos + 0.073, length=l_in, windings=n_in, current=current, r=r_in)
        self.create_coil(coil_type=self.coil_type, coil_mid_pos=first_coil_pos + 0.187, length=l_in, windings=n_in, current=current, r=r_in)
        self.create_coil(coil_type=self.coil_type, coil_mid_pos=first_coil_pos + 0.26, length=l_out, windings=n_out, current=-current, r=r_out)

    def create_mieze(self, current, l1=None, l2=None):
        if not l1:
            l1 = self.coil_distance
        
        if not l2:
            l2 = self.detector_distance - self.coil_distance
        
        i1 = current
        i2 = current*(l1+l2)/l2

        self._create_coil_set(current=i1)
        self._create_coil_set(current=i2, first_coil_pos=l1)
