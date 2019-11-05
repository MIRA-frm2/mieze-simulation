# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the MIEZE setup."""

from experiments.setup import Setup
from experiments.mieze.parameters import (DISTANCE_2ND_COIL, DISTANCE_3RD_COIL, DISTANCE_4TH_COIL,
                                          L_IN, L_OUT, N_IN, N_OUT, R_IN, R_OUT,
                                          SQUARE_COIL_POSITION_1ST, SQUARE_COIL_POSITION_2ND)


class Mieze(Setup):
    def __init__(self, sample_distance, coil_distance, detector_distance, increment=0.001):

        super(Mieze, self).__init__(increment=increment)

        self.sample_distance = sample_distance
        self.detector_distance = detector_distance
        self.coil_distance = coil_distance

        self.coil_type = 'simple'

    def _create_coil_set(self, first_coil_pos=0, current=5):
        """Crate the magnetic coil sets."""
        self.create_element(coil_type=self.coil_type, coil_mid_pos=first_coil_pos, length=L_OUT,
                            windings=N_OUT, current=-current, r=R_OUT)
        self.create_element(coil_type=self.coil_type, coil_mid_pos=first_coil_pos + DISTANCE_2ND_COIL,
                            length=L_IN, windings=N_IN, current=current, r=R_IN)
        self.create_element(coil_type=self.coil_type, coil_mid_pos=first_coil_pos + DISTANCE_3RD_COIL,
                            length=L_IN, windings=N_IN, current=current, r=R_IN)
        self.create_element(coil_type=self.coil_type, coil_mid_pos=first_coil_pos + DISTANCE_4TH_COIL,
                            length=L_OUT, windings=N_OUT, current=-current, r=R_OUT)

    def create_setup(self, current, l1=None, l2=None):
        if not l1:
            l1 = self.coil_distance
        
        if not l2:
            l2 = self.detector_distance - self.coil_distance
        
        current1 = current
        current2 = current * (l1+l2)/l2

        self.create_element(coil_type='square', coil_mid_pos=SQUARE_COIL_POSITION_1ST)

        self._create_coil_set(current=current1)
        self._create_coil_set(current=current2, first_coil_pos=l1)

        self.create_element(coil_type='square', coil_mid_pos=SQUARE_COIL_POSITION_2ND)
