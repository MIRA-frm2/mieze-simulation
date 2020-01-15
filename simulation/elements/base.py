# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Base element class."""

from abc import abstractmethod


class BasicElement(object):
    """Class implementing a basic experiment element."""

    def __init__(self, position, name):
        """Any physical element is supposed to have a position.

        Parameters
        ----------
        position: float, tuple
            If float, then it represents the x position.
            Otherwise, it is the 3d position as (x, y, z).
        """
        self.name = name

        if type(position) == tuple:
            self.position_x, self.position_y, self.position_z = position
        else:
            self.position_x = position
            self.position_y = 0
            self.position_y = 0

    @abstractmethod
    def b_field(self, r: '(x, y, z)'):
        raise NotImplementedError
