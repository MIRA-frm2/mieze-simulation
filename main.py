from experiments.mieze.main_mieze import Mieze
from experiments.mieze.parameters import ELEMENTS_POSITIONS_RELATIVE

from simulation.simulate import simulate

if __name__ == '__main__':
    simulate(experiment_class=Mieze, experiment_parameters=ELEMENTS_POSITIONS_RELATIVE)
