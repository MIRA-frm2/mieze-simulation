from experiments.mieze import Mieze
import matplotlib.pyplot as plt


def main():
    x = 1000
    y = 1000

    Experiment = Mieze(sample_distance=1.5, coils_dinstance=0.53, detector_distance=2.503, increment=0.001)
    Experiment.create_mieze(5)
    #print(Experiment.MIEZE_setup.setup_changed)

    plt.plot(x/1000., y/100.)
    Experiment.plot_field_1d_abs()
    # Experiment.plot_2d_vectormap()


if __name__ == "__main__":
    main()
