import logging
import matplotlib.pyplot as plt
import numpy as np
from numpy import gradient

from elements.coils import RealCoil
# from elements.polariser import Polariser
from elements.helmholtz_spin_flipper import HelmholtzSpinFlipper
from elements.spin_flipper import SpinFlipper

from experiments.mieze.parameters import I_sf1, I_hsf1, lambda_n, I_real_coil, HSF1, absolute_x_position, step

from utils.physics_constants import unit


# Create a custom logger
logger = logging.getLogger(__name__)


def extra_coil_check(b_extra, b_polariser):
    if b_extra.units != b_polariser.units:
        raise RuntimeError(f'B field units do not match:\n extra: {b_extra.units}\n'
                           f'polariser: {b_polariser.units}')

    if b_extra > b_polariser:
        pass
        # raise RuntimeError("Field of the extra coils too strong. Try resetting the extra coils.")


def plot(bx_values, by_values, b_values, theta_values, dtheta_dy):
    unit.setup_matplotlib(True)
    index_first_hsf = np.argmin(abs(absolute_x_position-HSF1))

    fig1, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Neutron Trajectory (m)')
    ax1.set_ylabel('Magnetic field (G)', color=color)

    # logger.error(f'{absolute_x_position}\n]{Bx_values}\n{By_values}')
    ax1.yaxis.set_units(unit.tesla)
    ax1.xaxis.set_units(unit.meter)

    ax1.plot(absolute_x_position, bx_values)
    ax1.plot(absolute_x_position, by_values)

    ax1.legend((r'$B_x$', r'$B_y$'), loc=9)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:green'
    x_pos = absolute_x_position[index_first_hsf]
    # noinspection PyTypeChecker
    y_pos = theta_values[index_first_hsf]

    ax2.yaxis.set_units(unit.degree)
    ax2.xaxis.set_units(unit.meter)

    # we already handled the x-label with ax1
    ax2.set_ylabel(r'$\theta$ (degree)', color=color)
    ax2.plot(absolute_x_position, theta_values, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.plot(x_pos, y_pos, color='black', marker='o')
    ax2.text(x_pos, y_pos*0.9, '{:.1f}°'.format(y_pos))

    # fig1.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    # plt.savefig('By_Bx.pdf')
    # plt.close()

    fig1, ax = plt.subplots()

    ax.xaxis.set_units(unit.meter)
    ax.yaxis.set_units(unit.degree)

    # ax.set_yscale('log')
    ax.plot(absolute_x_position, np.asarray(dtheta_dy) * 1e-2)  # y from m to cm
    ax.plot(absolute_x_position, 2.65*lambda_n.magnitude * np.asarray(b_values) * 1e-1)  # B from Gauss to mT,
    ax.legend((r'$\frac{d\theta}{dy}$', r'$2.65\lambda B$'))

    ax.set_xlabel('Neutron Trajectory (m)')
    ax.set_ylabel("(degrees/cm)")
    ax.grid()
    # fig1.tight_layout()
    plt.show()
    # plt.savefig('Adiabatic_Check.pdf')
    # plt.close()


def main():
    # polariser = Polariser()
    real_coil = RealCoil(coil_mid_pos=0.05 * unit.m, length=0.1 * unit.m, windings=100,
                         current=I_real_coil, r=0.05 * unit.m)
    spin_flipper = SpinFlipper(current=I_sf1)
    helmholtz_spin_flipper = HelmholtzSpinFlipper(current=I_hsf1)

    bx_values = list()
    by_values = list()
    b_values = list()
    theta_values = list()

    for x in absolute_x_position:
        # b_field_polariser = polariser.b_field(x, 0, 0)
        b_field_realcoil = real_coil.b_field(x, 0, 0)[0]
        # Check
        # extra_coil_check(b_field_realcoil, b_field_polariser)

        b_field_spin_flipper = spin_flipper.b_field(x, 0, 0)
        b_field_helmholtz = helmholtz_spin_flipper.b_field(x, 0, 0)[0]

        bx = b_field_spin_flipper  # G
        by = b_field_helmholtz + b_field_realcoil  # G

        logger.error(f'Magnetic fields: \n'
                     # f'polariser: {b_field_polariser}\n'
                     f'real coil: {b_field_realcoil}\n'
                     f'spinflipper:{b_field_spin_flipper}\n'
                     f'helmholtz {b_field_helmholtz}')

        logger.error(f'bx: {bx}\n'
                     f'by: {by}\n')

        b = np.sqrt(bx ** 2 + by ** 2)

        theta = np.degrees(np.arctan(by.magnitude/bx.magnitude))

        bx_values.append(bx.magnitude)
        by_values.append(by.magnitude)
        b_values.append(b.magnitude)

        logger.error(f'theta: {theta}')
        theta_values.append(abs(theta))

    logger.error(f'theta_values: {theta_values}')
    dtheta_dy = gradient(theta_values, step.magnitude)
    # logger.error(f'dtheta_dy: {dtheta_dy}')

    plot(bx_values, by_values, b_values, theta_values, dtheta_dy)


if __name__ == '__main__':
    main()
