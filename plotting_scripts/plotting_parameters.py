from elements.coils import Coil, SquareCoil

plot_parameters = {
    'coil_simple_2d_xy': {
        'element': Coil,
        'plot_dimension': '2d',
        'plot_args': {'type': 'scalar', 'component': 'x'},
        'grid_size': {'x_start': -0.2, 'x_end': 0.2, 'y_start': -0.2, 'y_end': 0.2, 'z_start': -0.0, 'z_end': 0.0,
                      'rho': 0.2, 'zoom_factor': 1},
        'coil_args': None,
    },
    'coil_simple_2d_yz': {
        'element': Coil,
        'plot_dimension': '2d',
        'plot_args': {'type': 'scalar', 'component': 'x'},
        'grid_size': {'x_start': -0.0, 'x_end': 0.0, 'y_start': -0.2, 'y_end': 0.2, 'z_start': -0.2, 'z_end': 0.2,
                      'rho': 0.2, 'zoom_factor': 1},
        'coil_args': None,
    },
    'coil_simple_3d_xz_plane': {
        'element': Coil,
        'plot_dimension': '3d',
        'plot_args': {'type': 'scalar', 'component': 'x'},
        'grid_size': {'x_start': -0.2, 'x_end': 0.2, 'y_start': -0.0, 'y_end': 0.0, 'z_start': -0.2, 'z_end': 0.2,
                      'zoom_factor': 1},
        'coil_args': None,
    },
    'coil_simple_3d_beamline': {
        'element': Coil,
        'plot_dimension': '3d',
        'plot_args': {'type': 'scalar', 'component': 'x'},
        'grid_size': {'x_start': -0.2, 'x_end': 0.2, 'y_start': -0.0, 'y_end': 0.0, 'z_start': -0.0, 'z_end': 0.0,
                      'zoom_factor': 1},
        'coil_args': None,
    },
    'coil_square_3d': {
        'element': SquareCoil,
        'plot_dimension': '3d',
        'plot_args': {'type': 'scalar', 'component': 'x'},
        'grid_size': {'x_start': -0.2, 'x_end': 0.2, 'y_start': -0.0, 'y_end': 0.0, 'z_start': -0.2, 'z_end': 0.2,
                      'zoom_factor': 1},
        'coil_args': None,
    },
    'coil_square_3d_beamline': {
        'element': SquareCoil,
        'plot_dimension': '3d',
        'plot_args': {'type': 'scalar', 'component': 'x'},
        'grid_size': {'x_start': -0.0, 'x_end': 0.0, 'y_start': -0.0, 'y_end': 0.0, 'z_start': -0.2, 'z_end': 0.2,
                      'zoom_factor': 1},
        'coil_args': None,
    }
}

# plot_parameters = {
#     'coil_simple_1d_abs': {
#         'element': Coil,
#         'plot_dimension': '1d',
#         'plot_args': {'type': 'scalar', 'component': 'x'},
#         'coil_args': None,
#         'b_field_args': {'start': -0.25, 'end': 1.5}
#     },
#     'coil_simple_1d_z': {
#         'element': Coil,
#         'plot_dimension': '1d',
#         'plot_args': {'type': 'scalar', 'component': 'x'},
#         'grid_size': {'start': -0.2, 'end': 0.2, 'rho': 0.2, 'zoom_factor': 1},
#         'coil_args': None,
#         'b_field_args': {'start': -0.25, 'end': 1.5}
#     },
#     'coil_simple_1d_vec': {
#         'element': Coil,
#         'plot_dimension': '1d',
#         'plot_args': {'type': 'vec', 'component': 'x'},
#         'coil_args': None,
#         'b_field_args': {'start': -0.25, 'end': 1.5}
#     },
#     'coil_square_1d_z': {
#         'element': SquareCoil,
#         'plot_dimension': '1d',
#         'plot_args': {'type': 'scalar', 'component': 'x'},
#         'coil_args': {'coil_mid_pos': 0.5, 'r': 0.2, 'length': 0.35},
#         'b_field_args': {'coordinate_system': 'cartesian', 'plane': 'xy', 'start': 0, 'end': 1}
#     },
#     'coil_square_1d_vec': {
#         'element': SquareCoil,
#         'plot_dimension': '1d',
#         'plot_type': 'vec',
#         'coil_args': {'coil_mid_pos': 0.5, 'r': 0.2, 'length': 0.35},
#         'b_field_args': {'coordinate_system': 'beamline','start': 0, 'end': 1}
#     },
#     'coil_square_2d_abs_yz': {
#         'element': SquareCoil,
#         'plot_dimension': '2d',
#         'plot_type': 'abs',
#         'coil_args': {'r': 0.2, 'length': 0.35},
#         'grid_size': {'start': 0, 'end': 0, 'rho': 0.2, 'zoom_factor': 2},
#         'b_field_args': {'coordinate_system': 'cartesian', 'plane': 'yz'}
#     },
#     'coil_square_2d_abs_xy': {
#         'element': SquareCoil,
#         'plot_dimension': '2d',
#         'plot_type': 'abs',
#         'coil_args': {'r': 0.2, 'length': 0.35},
#         'b_field_args': {'coordinate_system': 'cartesian', 'plane': 'xy'}
#     },
#     'coil_square_2d_abs_xz': {
#         'element': SquareCoil,
#         'plot_dimension': '2d',
#         'plot_type': 'abs',
#         'coil_args': {'r': 0.2, 'length': 0.35},
#         'b_field_args': {'coordinate_system': 'cartesian', 'plane': 'xz'}
#     }
# }
