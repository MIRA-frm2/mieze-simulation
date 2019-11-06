from plotting_scripts.plotting_parameters import plot_parameters

from experiments.setup import Setup
from utils.helper_functions import save_data_to_file


class Plot:

    def __init__(self, params):
        self.params = params

        self.setup = Setup(increment=0.01)
        self.setup.create_element(element_class=params['element'])

    def plot(self):
        self.setup.calculate_b_field(**self.params['b_field_args'])
        if self.params['plot_dimension'] == '1d' and self.params['plot_type'] == 'abs':
            self.setup.plot_field_1d_abs()
        elif self.params['plot_dimension'] == '2d' and self.params['plot_type'] == 'abs':
            self.setup.plot_field_2d_abs(self.params['b_field_args']['plane'])

    def save_data(self, to_save=False):
        if to_save:
            save_data_to_file(self.setup.b, '../data/data')

    # def plot_1d_map(self, b_field_kwargs):
    #     self.setup.calculate_b_field(**b_field_kwargs)
    #
    #     save_data_to_file(self.setup.b, '../../data/data')
    #
    #     self.setup.plot_field_1d_abs()
    #
    # def plot_field_2d_abs(self):
    #     self.setup.plot_field_2d_abs(plane=self.plane)
    #
    # def plot_field_2d_vec(self):
    #     self.setup.plot_field_2d_vec()
    #
    # def get_coil_args(self):
    #     if self.coil_type == Coil:
    #         coil_kwargs = {'element_class': self.coil_type}
    #         b_field_kwargs = {'coordinate_system': 'cylindrical', 'coil_type': self.coil_type, 'plane': self.plane, 'rho': 0.25,
    #                           'start': -0.25, 'end': 1.5}
    #     elif self.coil_type == 'real':
    #         coil_kwargs = {'element_class': self.coil_type}
    #         b_field_kwargs = {'coordinate_system': 'cylindrical', 'coil_type': self.coil_type, 'plane': self.plane, 'rho': 0.25,
    #                           'start': -0.25, 'end': 1.5}
    #     elif self.coil_type == SquareCoil:
    #         coil_kwargs = {'element_class': self.coil_type, 'r': 0.2, 'length': 0.35}
    #         b_field_kwargs = {'coordinate_system': 'cartesian', 'coil_type': self.coil_type}
    #     else:
    #         coil_kwargs = None
    #         b_field_kwargs = None
    #     return coil_kwargs, b_field_kwargs
    #
    # def plot_simple_coil_field_1d_abs(self):
    #     b_field_kwargs = {'coordinate_system': 'cylindrical', 'coil_type': self.coil_type, 'plane': self.plane,
    #                       'rho': 0.25, 'start': -0.25, 'end': 1.5}
    #     self.setup.calculate_b_field(**b_field_kwargs)
    #
    #     self.setup.plot_field_1d_abs()
    #
    # def plot_simple_coil_field_2d_abs_xy(self):
    #     plane = 'xy'
    #     b_field_kwargs = {'coordinate_system': 'cylindrical', 'coil_type': self.coil_type, 'plane': plane,
    #                       'rho': 0.25, 'start': -0.25, 'end': 1.5}
    #     self.setup.calculate_b_field(**b_field_kwargs)
    #
    #     self.setup.plot_field_2d_abs(plane=plane)
    #
    # def plot_simple_coil_field_2d_abs_yz(self):
    #     plane = 'yz'
    #     b_field_kwargs = {'coordinate_system': 'cylindrical', 'coil_type': self.coil_type, 'plane': plane,
    #                       'rho': 0.25, 'start': -0.25, 'end': 1.5}
    #     self.setup.calculate_b_field(**b_field_kwargs)
    #
    #     self.setup.plot_field_2d_abs(plane=plane)


if __name__ == "__main__":
    name = 'coil_square_1d_abs'

    plot = Plot(params=plot_parameters[name])
    plot.plot()

