from experiments.mieze.main import Mieze


def main():
    Experiment = Mieze(sample_distance=1.5, coils_dinstance=0.53, detector_distance=2.503, increment=0.01)
    Experiment.create_mieze(current=5)

    rho = 0.02
    start = -0.25
    end = Experiment.sample_distance

    b_field_kwargs = {'coordinate_system': 'cylindrical', 'plane': 'xy', 'rho': 0.25,
                      'start': -0.25, 'end': 1.5}
    Experiment.calculate_b_field(**b_field_kwargs)

    Experiment.plot_field_1d_abs()
    # Experiment.plot_field_2d_abs()
    # Experiment.plot_field_1d_vec()
    # Experiment.plot_2d_vectormap()


if __name__ == "__main__":
    main()
