# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script to be executed."""


from analysises.neutron_polarisation_simulation.scripts.neutron_pol_sim import compute_neutron_beam
from analysises.neutron_polarisation_simulation.scripts.plotting_scripts import plot_polarisation_vector


def main():
    compute_neutron_beam()
    plot_polarisation_vector(normalize=True, length=0.025)


if __name__ == '__main__':
    main()
