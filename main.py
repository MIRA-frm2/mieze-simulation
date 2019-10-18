from experiments.mieze import Mieze
import matplotlib.pyplot as plt


def main():
    x = 1000
    y = 1000

    Experiment = Mieze(1.5, 0.53, 2.503, 0.001)
    Experiment.create_MIEZE(5)
    #print(Experiment.MIEZE_setup.setup_changed)
    plt.plot(x/1000.,y/100.)
    Experiment.plot_field_1D_abs()


if __name__ == "__main__":
    main()
