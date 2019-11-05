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

    def __init__(self):
        self.position = None

    @abstractmethod
    def b_field(self, x, y, z):
        raise NotImplementedError
