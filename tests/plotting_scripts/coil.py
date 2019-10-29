from setup_elements.setup import Setup


def main():
    coil_type = 'simple'
    setup = Setup(increment=0.01, coil_type=coil_type)

    if coil_type == 'simple':
        coil_kwargs = {}
        b_field_kwargs = {'rho': 0.25, 'start': -0.25, 'end': 1.5, 'coordinate_system': 'cylindrical'}
    elif coil_type == 'square':
        coil_kwargs = {'r': 0.2, 'length': 0.35}
        b_field_kwargs = {'coordinate_system': 'cartesian', 'start': -0.25}

    else:
        coil_kwargs = None
        b_field_kwargs = None

    setup.create_coil(**coil_kwargs)

    setup.calculate_b_field(**b_field_kwargs)

    # setup.plot_1d_abs()
    setup.plot_2d_map(plane='yz')


if __name__ == "__main__":
    main()
