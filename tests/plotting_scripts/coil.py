from setup_elements.setup import Setup
from setup_elements.helper_functions import save_data_to_file


def main():
    coil_type = 'simple'
    plane = 'yz'

    setup = Setup(increment=0.05, coil_type=coil_type)

    if coil_type == 'simple':
        coil_kwargs = {}
        b_field_kwargs = {'coordinate_system': 'cylindrical', 'plane': plane, 'rho': 0.25, 'start': -0.25, 'end': 1.5}
    elif coil_type == 'square':
        coil_kwargs = {'r': 0.2, 'length': 0.35}
        b_field_kwargs = {'coordinate_system': 'cartesian', 'start': -0.25, 'plane': plane}

    else:
        coil_kwargs = None
        b_field_kwargs = None

    setup.create_coil(**coil_kwargs)

    setup.calculate_b_field(**b_field_kwargs)

    save_data_to_file(setup.b, '../../data/data')

    setup.plot_1d_abs()
    # setup.plot_2d_map(plane=plane)


if __name__ == "__main__":
    main()
