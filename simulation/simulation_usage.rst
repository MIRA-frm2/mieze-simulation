MIEZE Simulation
################

The parameter space for computing the magnetic field of the MIEZE setup is represented by 
`this file  <experimental_setup/mieze/parameters.py>`_.

Additional user input consists of defining the 3d computational grid in the `main file <simulate.py>`_.

In the following, each file is shortly described. For its usage, please consult also the documentation in each file.

The `main script <simulation/simulate.py>`_ computes the magnetic field for the MIEZE experiment, which is saved in the
`data_magnetifc_field <./data/data_magnetic_field.csv>`_ file.

The folder `analysis_scripts <analysises>`_ contains several analysis scripts that have been used to
find the ideal position for several elements of the experiment, or the check whether other properties are satisfied or
not.

* `Adiabatic Transition <analysises/adiabatic_check/scripts/adiabatic_check.py>`_: Plots the magnetic field and the adiabatic transition to check whether it is fulfilled or not.
* `CoilSet Position <analysises/coil_set_configuration/scripts/coil_set_positions.py>`_: Investigates what is the prefered position for the outer two coils such that the total magnetic field of the four coils is as close as possible to a unit function.
* `Neutron Polarisation Simulation <analysises/neutron_polarisation_simulation/neutron_pol_sim.py>`_: Computes the neutrons trajectories along the beam and their polarisation. It saved the data in the `data_polarisation <./data/data_polarisation.csv>`_ file.
* `Polariser Adjustment <analysises/polariser_adjustment/polariser_adjustment.py>`_: Interpolates the measured polariser data points and optimises its magnetic dipole function.


The folder `elements <simulation/elements>`_ contains the individual elements that are placed along the beam, such as
the beam trajectory, such as the:

* `Coils <simulation/elements/coils.py>`_: Circular (Simple and Real) and Rectangular
* `CoilSet <simulation/elements/coil_set.py>`_: The set of four coils for the mieze setup
* `HelmholtzPair <simulation/elements/helmholtz_pair.py>`_: The pair of two coils in Helmholtz condition.
* `Polariser <simulation/elements/coils.py>`_: The Polariser (similar to a dipole>)
* `SpinFlipper <simulation/elements/spin_flipper.py>`_: The Pi/2 Spin Flipper.