Adiabatic Transition
====================

One crucial aspect of the MIEZE setup is that the polarization between the polariser and the helmholtz coils proceeds
adiabatically (theadiabatic transition condition can be found in the literature, but for convenience we provide 
[two slides](./adiabatic_transition_literature.pdf) about it). This investigation consists of two scripts:

1. A [check](./adiabatic_check.py) using the existing parameters.

    This condition is then plotted together with the actual change of the angle in the 
    [Adiabatic transition plot](./adiabatic_transition_condition.png).
    
    The magnetic field along x and y direction are plotted together with the angle change in [this plot](./by_bx.png).

2. An [investigation](./find_relative_distance.py) trying to find what the optimal parameters are.
    
    The result is plotted [here](./distance_investigation_for_the_adiabatic_transition_condition.png).
    
    Because the plot would be overcrowded if plotted in a similar fashion as for part 1, we define an adiabatic 
    transition discriminant from the adiabatic condition inequality as just:

    $2.65 \cdot B \cdot \lambda - \dfrac{d\theta}{dy}$

    This quantity has to always be greater than 0 to fulfill the adiabatic transition.
