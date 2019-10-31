from positions import Positions
import numpy as np
from real_coil_skript import Real_Coil


class Helmholtz_Spin_Flipper:
    MU_0 = 4e-7 * np.pi #N/A, vacuum permeability
    R = 0.055 #Positions.R_HSF # Radius of each coil
    N = 33 # Windigs
    l = 0.01 # length of each coil
    d = 0.045 #coil distance
    def __init__(self):
        pass

    @classmethod
    def hsf(cls, x, y, z, mid_pos=Positions.HSF1, current=1.6):
        pos1 = mid_pos -cls.d/2.0

        pos2 = mid_pos +cls.d/2.0

        B1 = Real_Coil.B_field(x, y, z, pos1, l=cls.l, N=cls.N, R=cls.R )

        B2 = Real_Coil.B_field(x, y, z, pos2, l=cls.l, N=cls.N, R=cls.R )



        return (B1 + B2)

    