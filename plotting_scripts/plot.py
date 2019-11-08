from plotting_scripts.plotting_parameters import plot_parameters

from experiments.setup import Setup
from utils.helper_functions import save_data_to_file


class Plot:

    def __init__(self, params):
        self.params = params

        self.setup = Setup(increment=0.05)
        self.setup.create_element(element_class=params['element'])

    def plot(self):
        self.setup.initialize_computational_space(**self.params['b_field_args'])
        self.setup.calculate_b_field()

        if self.params['plot_dimension'] == '1d':
            if self.params['plot_type'] == 'abs':
                self.setup.plot_field_1d_abs()
            elif self.params['plot_type'] == 'vec':
                self.setup.plot_field_1d_vec()
        elif self.params['plot_dimension'] == '2d' and self.params['plot_type'] == 'abs':
            self.setup.plot_field_2d_abs(self.params['b_field_args']['plane'])

    def save_data(self, to_save=False):
        if to_save:
            save_data_to_file(self.setup.b, '../data/data')


if __name__ == "__main__":
    name = 'coil_square_1d_abs'

    plot = Plot(params=plot_parameters[name])
    plot.plot()

