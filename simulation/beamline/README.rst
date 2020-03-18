Beamline
********

The flux of neutrons is simulated using the `beam  file <https://github.com/MIRA-frm2/mieze-simulation/blob/master/experiments/mieze/parameters.py>`_

Algorithm features:

* It creates a discrete computational space/grid where neutrons can fly through. Ideally, it should be identical with
  the one used for computing the magnetic field.

* Neutrons are then created at the beginning ("left side", indexes [0][0][0]) of the computational grid. Alternatively,
  the user can defined a 'starting position_x'.

* The neutrons have a starting velocity vector (also with radial component), if the 'distribution' parameter is passed.
  Otherwise, they only propagate in the beamline (x) direction.

* The time spent by each neutron in each cell is computed based on the cell size and the neutron velocity.

* The polarisation is adjusted based on the time spent in the cell.

* If a neutron flies out of the computational grid (it can either reach the end of the beam or diverge from the beam),
  then it is removed from the simulation.
