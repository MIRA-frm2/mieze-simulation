from experiments.setup import Setup
from utils.helper_functions import save_data_to_file


def main():
    coil_type = 'real'
    plane = 'yz'

    setup = Setup(increment=0.05)

    if coil_type == 'simple':
        coil_kwargs = {}
        b_field_kwargs = {'coordinate_system': 'cylindrical', 'coil_type': coil_type, 'plane': plane, 'rho': 0.25, 'start': -0.25, 'end': 1.5}
    elif coil_type == 'real':
        coil_kwargs = {}
        b_field_kwargs = {'coordinate_system': 'cylindrical', 'coil_type': coil_type, 'plane': plane, 'rho': 0.25, 'start': -0.25, 'end': 1.5}
    elif coil_type == 'square':
        coil_kwargs = {'r': 0.2, 'length': 0.35}
        b_field_kwargs = {'coordinate_system': 'cartesian', 'coil_type': coil_type, 'plane': plane, 'start': -0.25, }
    else:
        coil_kwargs = None
        b_field_kwargs = None

    setup.create_element(**coil_kwargs)

    setup.calculate_b_field(**b_field_kwargs)

    save_data_to_file(setup.b, '../../data/data')

    # setup.plot_1d_abs()
    setup.plot_field_2d_abs(plane=plane)


if __name__ == "__main__":
    main()
