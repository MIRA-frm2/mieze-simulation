Analysises
##########

The folder `analysis_scripts <https://github.com/dprelipcean/mieze-simulation/blob/master/analysises>`_ contains several analysis scripts that have been used to
find the ideal position for several elements of the experiment, or the check whether other properties are satisfied or
not.

* `Adiabatic Transition <https://github.com/dprelipcean/mieze-simulation/blob/master/analysises/adiabatic_check/scripts/adiabatic_check.py>`_: Plots the magnetic field and the adiabatic transition to check whether it is fulfilled or not.
* `CoilSet Position <https://github.com/dprelipcean/mieze-simulation/blob/master/analysises/coil_set_configuration/scripts/coil_set_positions.py>`_: Investigates what is the prefered position for the outer two coils such that the total magnetic field of the four coils is as close as possible to a unit function.
* `Neutron Polarisation Simulation <https://github.com/dprelipcean/mieze-simulation/blob/master/nalysises/neutron_polarisation_simulation/neutron_pol_sim.py>`_: Computes the neutrons trajectories along the beam and their polarisation. It saved the data in the `data_polarisation <./data/data_polarisation.csv>`_ file.
* `Polariser Adjustment <https://github.com/dprelipcean/mieze-simulation/blob/master/analysises/polariser_adjustment/polariser_adjustment.py>`_: Interpolates the measured polariser data points and optimises its magnetic dipole function.



.. include:: ../../analysises/adiabatic_check/README.md

.. include:: ../../analysises/adiabatic_polarisation/README.md

.. include:: ../../analysises/coil_set_configuration/README.md

.. include:: ../../analysises/coilset_influence_on_hsf/README.md

.. include:: ../../analysises/neutron_polarisation_simulation/README.md

.. include:: ../../analysises/numerical_inconsistencies/README.md

.. include:: ../../analysises/polariser_adjustment/README.md
