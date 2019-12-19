Coil Set Configuration
======================

Each of the inner pail of coils creates a magnetic field that is roughly similar to a Gaussian. Adding two such fields
with a small distance between them  almost produces a [resulting magnetic field](results/bfield_coils_inner.png) similar to a 
pedestal peak. In order to optimize this pedestal (that is, a flat top with straight edges), 
[two outer coils](results/bfield_coils_inner.png) with opposite current are used to remove and the tails.

This script investigates what should the position of the outer coils relative to the inner coils should be to achieve a
pedestal as good as possible, by iteratively moving the outer coils further away from the inner ones. The goodness of 
the fit is computed using a minimal squares approach, and the [result](results/coil_set_optimization.png) is plotted a 2d 
color map.

The reference function is taken to be a perfect pedestal, as in [two outer coils](results/bfield_coil_set_0.png).

Results
=======
The conclusion of this analysis is that the outer coils should be as close as possible to the inner coils.