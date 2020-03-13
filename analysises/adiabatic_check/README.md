Adiabatic Transition
********************

One crucial aspect of the MIEZE setup is that the polarization between the polariser and the helmholtz coils proceeds
adiabatically (the adiabatic transition condition can be found in the literature, but for convenience we provide 
[two slides](docs/adiabatic_transition_literature.pdf) about it). 

The way the scripts work is that they take existing magnetic field data values from .csv files, then compute all 
required variables and perform the plots.

This investigation is performed via an over-layer script called [perform_analysis.py](perform_analysis.py), where the 
user decides which analysis to perform. The two scripts to choose form are:

1. A [check](scripts/adiabatic_check.py) using the existing parameters.

    This condition is then plotted together with the actual change of the angle in the 
    [Adiabatic transition plot](results/adiabatic_transition_condition.png).
    
    The magnetic field along x and y direction are plotted together with the angle change in [this plot](results/by_bx.png).

    There is no user input for this analysis.

2. An [investigation](scripts/find_relative_distance.py) trying to find what the optimal parameters are.
    
    The result is plotted [here](results/distance_investigation_for_the_adiabatic_transition_condition.png).
    
    Because the plot would be overcrowded if plotted in a similar fashion as for part 1, we define an adiabatic 
    transition discriminant from the adiabatic condition inequality as just:

    $2.65 \cdot B \cdot \lambda - \dfrac{d\theta}{dy}$

    This quantity has to always be greater than 0 to fulfill the adiabatic transition.
    
    The user input for this analysis consits of an array of values to be iterated over for finding the ideal distance.
