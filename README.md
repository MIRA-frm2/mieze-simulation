MIEZE Simulation
================

About
=====

This project simulates the MIEZE geometry together, focusing on the magnetic field computation and neutron polarization 
evolution along the trajectory.

Usage
=====

In the following, each file is shortly described. For its usage, please consult also the documentation in each file.

The [main script](simulation/simulate.py) computes the magnetic field for the MIEZE experiment, which is saved in the 
[data_magnetifc_field](./data/data_magnetic_field.csv) file.

The folder [analysis_scripts](analysises) contains several analysis scripts that have been used to 
find the ideal position for several elements of the experiment, or the check whether other properties are satisfied or 
not.

* [Adiabatic Transition](analysises/adiabatic_check/scripts/adiabatic_check.py): 
Plots the magnetic field and the adiabatic transition to check whether it is fulfilled or not.
* [CoilSet Position](analysises/coil_set_configuration/scripts/coil_set_positions.py): 
Investigates what is the prefered position for the outer two coils such that the total magnetic field of the four coils 
is as close as possible to a unit function.
* [Neutron Polarisation Simulation](analysises/neutron_polarisation_simulation/neutron_pol_sim.py): 
Computes the neutrons trajectories along the beam and their polarisation. It saved the data in the
[data_polarisation](./data/data_polarisation.csv) file. 
* [Polariser Adjustment](analysises/polariser_adjustment/polariser_adjustment.py):
Interpolates the measured polariser data points and optimises its magnetic dipole function.


The folder [elements](simulation/elements) contains the individual elements that are placed along the beam, such as 
the beam trajectory, such as the:

* [Coils](simulation/elements/coils.py): Circular (Simple and Real) and Rectangular 
* [CoilSet](simulation/elements/coil_set.py): The set of four coils for the mieze setup
* [HelmholtzPair](simulation/elements/helmholtz_pair.py): The pair of two coils in Helmholtz condition.
* [Polariser](simulation/elements/coils.py): The Polariser (similar to a dipole)
* [SpinFlipper](simulation/elements/spin_flipper.py): The Pi/2 Spin Flipper.
 

Documentation Style
===================

Most documentation is written withing the code base, such that it's easily accessible and it matches the actual scripts.
