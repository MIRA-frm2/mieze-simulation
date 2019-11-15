from plotting_scripts.plotting_parameters import plot_parameters

from experiments.setup import Setup
from utils.helper_functions import save_data_to_file


def compute_bfield(params):
    setup = Setup(increment=0.01)
    setup.create_element(element_class=params['element'])

    setup.initialize_computational_space(**params['grid_size'])
    setup.calculate_b_field()

    save_data_to_file(setup.b, '../data/data')

    # if self.params['plot_dimension'] == '1d':
    #     if self.params['plot_args']['type'] == 'scalar':
    #         self.setup.plot_field_1d_scalar(self.params['plot_args']['component'])
    #     elif self.params['plot_args']['type'] == 'vec':
    #         self.setup.plot_field_1d_vec()
    # elif self.params['plot_dimension'] == '2d' and self.params['plot_type'] == 'abs':
    #     self.setup.plot_field_2d_abs(self.params['b_field_args']['plane'])


if __name__ == "__main__":
    name = 'coil_square_2d_xy'

    compute_bfield(params=plot_parameters[name])

