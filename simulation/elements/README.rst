Elements
********

The folder `elements <simulation/elements>`_ contains the individual elements that are placed along the beam, such as
the beam trajectory, such as the:

* `Coils <https://github.com/MIRA-frm2/mieze-simulation/blob/master/simulation/elements/coils.py>`_: Circular (Simple and Real) and Rectangular
* `CoilSet <https://github.com/MIRA-frm2/mieze-simulation/blob/master/simulation/elements/coil_set.py>`_: The set of four coils for the mieze setup
* `HelmholtzPair <https://github.com/MIRA-frm2/mieze-simulation/blob/master/simulation/elements/helmholtz_pair.py>`_: The pair of two coils in Helmholtz condition.
* `Polariser <https://github.com/MIRA-frm2/mieze-simulation/blob/master/simulation/elements/coils.py>`_: The Polariser (similar to a dipole>)
* `SpinFlipper <https://github.com/MIRA-frm2/mieze-simulation/blob/master/simulation/elements/spin_flipper.py>`_: The Pi/2 Spin Flipper.

All elements are derived from a `Base class <https://github.com/MIRA-frm2/mieze-simulation/blob/master/simulation/elements/spin_flipper.py>`_,
containing the abstract method for computing the magnetic field.

To use these elements, their specific keyword arguments have to be provided. It is recommended that this done via one
exhaustive file per experiment, e.g. the `MIEZE parameters file <https://github.com/MIRA-frm2/mieze-simulation/blob/master/experiments/mieze/parameters.py>`_
