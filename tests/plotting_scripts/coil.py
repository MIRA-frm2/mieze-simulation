from setup_elements.setup import Setup


def main():
    setup = Setup(increment=0.01, coil_type='simple')
    setup.create_coil()

    rho = 0.02
    start = -0.25
    end = 1.5

    setup.calculate_b_field(start, end, rho=rho)

    # setup.plot_1d_abs()
    setup.plot_2d_map()
    # setup.plot_3d()


if __name__ == "__main__":
    main()
