MIEZE Simulation
################

The parameter space for computing the magnetic field of the MIEZE setup is represented by 
`this file  <https://github.com/dprelipcean/mieze-simulation/blob/master/experiments/mieze/parameters.py>`_.

Additional user input consists of defining the 3d computational grid in the `main file <simulate.py>`_.

In the following, each file is shortly described. For its usage, please consult also the documentation in each file.

The `main script <https://github.com/dprelipcean/mieze-simulation/tree/master/simulation/simulate.py>`_ computes the magnetic field for the MIEZE experiment, which is saved in the
`data_magnetic_field <https://github.com/dprelipcean/mieze-simulation/blob/master/data/data_magnetic_field.csv>`_ file.


.. include:: ../../simulation/elements/README.rst

.. include:: ../../simulation/experimental_setup/README.rst

.. include:: ../../simulation/particles/README.rst
