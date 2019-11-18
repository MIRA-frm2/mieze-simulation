import logging
import matplotlib.pyplot as plt
import numpy as np
from numpy import gradient

from elements.coils import RealCoil
from elements.polariser import Polariser
from elements.helmholtz_spin_flipper import HelmholtzSpinFlipper
from elements.spin_flipper import SpinFlipper

from experiments.mieze.parameters import I_sf1, I_hsf1, lambda_n, I_real_coil, endpoint, HSF1, absolute_x_position, step

from utils.physics_constants import unit


# Create a custom logger
logger = logging.getLogger(__name__)


def extra_coil_check(B_extra, B_polariser):
    if B_extra.units != B_polariser.units:
        raise RuntimeError(f'B field units do not match:\n extra: {B_extra[0].units}\n'
                           f'polariser: {B_polariser.units}')

    if B_extra > B_polariser:
        pass
        # raise RuntimeError("Field of the extra coils too strong. Try resetting the extra coils.")


def main():
    # Define Polariser
    polariser = Polariser()
    # Define Real Coil
    real_coil = RealCoil(coil_mid_pos=0.05 * unit.m, length=0.1 * unit.m, windings=100,
                         current=I_real_coil, r=0.05 * unit.m)
    # Define Spin Flippers and Helmholtzcoils
    spin_flipper = SpinFlipper(current=I_sf1)

    helmholtz_spin_flipper = HelmholtzSpinFlipper(current=I_hsf1)

    Bx_values = list()
    By_values = list()
    B_values = list()
    theta_values = list()
    dtheta_dy_values = list()

    for x in absolute_x_position:
        b_field_polariser = polariser.b_field(x, 0, 0)
        b_field_realcoil = real_coil.b_field(x, 0, 0)[1]
        # Check
        extra_coil_check(b_field_realcoil, b_field_polariser)

        b_field_spin_flipper = spin_flipper.b_field(x, 0, 0)
        b_field_helmholtz = helmholtz_spin_flipper.b_field(x, 0, 0)[1]

        Bx = b_field_polariser + b_field_spin_flipper  # G
        By = b_field_helmholtz + b_field_realcoil  # G

        logger.error(f'polariser{b_field_polariser}\nreal coil {b_field_realcoil}'
                     f'\nspinflipper{b_field_spin_flipper}\nheilholtz {b_field_helmholtz}')
        logger.error(f'bx{Bx}\nby{By}')

        B = np.sqrt(Bx**2+By**2)

        theta = np.degrees(np.arctan(By/Bx))
        dtheta_dy = gradient(theta, step)

        Bx_values.append(Bx)
        By_values.append(By)
        B_values.append(B)

        theta_values.append(theta)
        dtheta_dy_values.append(dtheta_dy)

    index_first_hsf = np.argmin(abs(absolute_x_position-HSF1))

    fig1, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Neutron Trajectory (m)')
    ax1.set_ylabel('Magnetic field (G)', color=color)
    ax1.plot(absolute_x_position, Bx_values)
    ax1.plot(absolute_x_position, By_values)
    ax1.legend((r'$B_x$', r'$B_y$'), loc=9)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:green'
    x_pos = absolute_x_position[index_first_hsf]
    y_pos = theta_values[index_first_hsf]
    ax2.set_ylabel(r'$\theta$ (degree)', color=color)  # we already handled the x-label with ax1
    ax2.plot(absolute_x_position, theta_values, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.plot(x_pos, y_pos, color='black', marker='o')
    ax2.text(x_pos, y_pos*0.9, '{:.1f}°'.format(y_pos))

    fig1.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    plt.savefig('By_Bx.pdf')
    # plt.close()

    plt.figure()
    plt.plot(absolute_x_position, dtheta_dy_values * 1e-2)  # y from m to cm
    plt.plot(absolute_x_position, 2.65*lambda_n*B_values*1e-1)  # B from Gauss to mT,
    plt.legend((r'$\frac{d\theta}{dy}$', r'$2.65\lambda B$'))
    plt.xlabel('Neutron Trajectory (m)')
    plt.ylabel("(degrees/cm)")
    plt.grid()
    plt.tight_layout()
    plt.show()
    plt.savefig('Adiabatic_Check.pdf')
    # plt.close()


if __name__ == '__main__':
    main()
