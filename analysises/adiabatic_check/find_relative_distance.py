import matplotlib.pyplot as plt
import numpy as np

from main import main
from analysises.adiabatic_check.adiabatic_check import MyPlotter



if __name__ == "__main__":
    iteration_values = np.arange(0.0, 0.2 + 0.025, step=0.025)

    diff_data_list = list()

    for val in iteration_values:
        data_file = f'./data/test_data_{int(val*100)}.csv'

        main(filename=data_file, spin_flipper_distance=val)

        plot = MyPlotter(data_file=data_file)

        plot.get_adiabatic_values(new=True)
        plot.get_adiabatic_difference_data()

        diff_data_list.append(plot.diff_data)

    for item in diff_data_list:
        plt.plot(plot.x_range, item)

    plt.xlabel('Neutron Trajectory (m)')
    plt.ylabel(r'Adiabatic condition discriminant: $2.65 \cdot B \cdot \lambda - \dfrac{d\theta}{dy}$')

    plt.yscale('log')
    plt.legend(np.asarray(iteration_values) * 100,
               title="Relative distance between polariser and helmholtz coils [cm]",
               loc=1)
    plt.show()

