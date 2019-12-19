from simulation.plotting_scripts.plotting_parameters import plot_parameters

from simulation.experiments.setup import Setup
from utils.helper_functions import save_data_to_file, save_obj


def compute_bfield(params):
    setup = Setup()

    setup.create_element(**params['element_kwargs'])

    setup.initialize_computational_space(**params['grid_size'])
    setup.calculate_b_field()

    save_data_to_file(setup.b, '../data/data')
    save_obj(setup.b, '../data/data')


if __name__ == "__main__":
    name = 'neutron_beam_axis'

    compute_bfield(params=plot_parameters[name])

