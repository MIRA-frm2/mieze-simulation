MIEZE SetUp
***********

The implementation of the MIEZE setup is in the `main_mieze.py file <https://github.com/MIRA-frm2/mieze-simulation/blob/master/experiments/mieze/main_mieze.py>`_.
It allows to also compute the magnetic field value at specific points by passing a 3D array/list as the 'point_value'
keyword argument, eg. 'point_value = (0, 0, 0)'.

The `parameters file <https://github.com/MIRA-frm2/mieze-simulation/blob/master/experiments/mieze/parameters.py>`_
contains ALL the information specific to the MIEZE setup. This is where changes should be done with respect to distances,
current values, etc.

The parameters for each element are contained in a dictionary entitled 'PARAMETERS_<NAME>', e.g. 'PARAMETERS_COIL_SET'.
In particular, every dictionary has the 'USE' key that should either be 'True' or 'False' to indicate whether this
specific element should be used in the simulation or not.
