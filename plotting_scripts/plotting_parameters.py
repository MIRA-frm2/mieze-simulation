from elements.coils import Coil, SquareCoil

plot_parameters = {
    'coil_simple_1d_abs': {
        'element': Coil,
        'plot_dimension': '1d',
        'plot_type': 'abs',
        'coil_args': None,
        'b_field_args': {'coordinate_system': 'cylindrical', 'plane': 'xy', 'rho': 0.25, 'start': -0.25, 'end': 1.5}
    },
    'coil_simple_1d_vec_beam': {
        'element': Coil,
        'plot_dimension': '1d',
        'plot_type': 'vec',
        'coil_args': None,
        'b_field_args': {'coordinate_system': 'beamline', 'plane': 'xy', 'rho': 0.25, 'start': -0.25, 'end': 1.5}
    },
    'coil_square_1d_abs': {
        'element': SquareCoil,
        'plot_dimension': '1d',
        'plot_type': 'abs',
        'coil_args': {'coil_mid_pos': 0.5, 'r': 0.2, 'length': 0.35},
        'b_field_args': {'coordinate_system': 'cartesian', 'plane': 'xy', 'start': 0, 'end': 1}
    },
    'coil_square_2d_abs_yz': {
        'element': SquareCoil,
        'plot_dimension': '2d',
        'plot_type': 'abs',
        'coil_args': {'r': 0.2, 'length': 0.35},
        'b_field_args': {'coordinate_system': 'cartesian', 'plane': 'yz'}
    },
    'coil_square_2d_abs_xy': {
        'element': SquareCoil,
        'plot_dimension': '2d',
        'plot_type': 'abs',
        'coil_args': {'r': 0.2, 'length': 0.35},
        'b_field_args': {'coordinate_system': 'cartesian', 'plane': 'xy'}
    },
    'coil_square_2d_abs_xz': {
        'element': SquareCoil,
        'plot_dimension': '2d',
        'plot_type': 'abs',
        'coil_args': {'r': 0.2, 'length': 0.35},
        'b_field_args': {'coordinate_system': 'cartesian', 'plane': 'xz'}
    }
}
