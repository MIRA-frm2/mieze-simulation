from setup_elements.setup import Setup


class Mieze:

    def __init__(self, sample_distance, coils_dinstance, detector_distance, increment=0.001):
        self.MIEZE_setup = Setup(increment)
        self.sample_distance = sample_distance
        self.detector_distance = detector_distance
        self.coil_distance = coils_dinstance

    def _create_coil_set(self, firstCoilPos=0, I=5):
        l_in = 0.086
        l_out = 0.05
        N_in = 168
        N_out = 48
        R_in = 0.177 / 2.
        R_out = 0.13

        self.MIEZE_setup.create_coil(firstCoilPos, length=l_out, windings=N_out, current=-I, r=R_out)
        self.MIEZE_setup.create_coil(firstCoilPos + 0.073, length=l_in, windings=N_in, current=I, r=R_in)
        self.MIEZE_setup.create_coil(firstCoilPos + 0.187, length=l_in, windings=N_in, current=I, r=R_in)
        self.MIEZE_setup.create_coil(firstCoilPos + 0.26, length=l_out, windings=N_out, current=-I, r=R_out)

    def create_MIEZE(self, I, L1 = None, L2 = None):
        if L1 == None:
            L1 = self.coil_distance
        
        if L2 == None:
            L2 = self.detector_distance - self.coil_distance
        
        I1 = I
        I2 = I*(L1+L2)/L2

        self._create_coil_set(I=I1)
        self._create_coil_set(I=I2, firstCoilPos=L1)

    def plot_field_1D_abs(self, rho=0):
        self.MIEZE_setup.plot_1D_abs(-0.25, self.sample_distance, rho)

    def plot_field_1D_vec(self, rho=0):
        self.MIEZE_setup.plot_1D_vector(-0.25, self.sample_distance, rho)
