import numpy as np
from vasprun import units

class IR:
    """
    A class to compute the IR intensity
    """
    def __init__(self, born_chgs, eigenvalues, eigenvectors, mass, vol):
        self.born_chgs = np.array(born_chgs)
        self.modes = eigenvectors
        self.freqs = eigenvalues
        self.Natom = len(born_chgs)
        IRs = []
        epsilons = []
        for mode, freq in zip(self.modes, self.freqs):
            IRx, IRy, IRz = 0, 0, 0
            mode = np.reshape(mode, [self.Natom, 3])
            for i in range(self.Natom):
                for j in range(3):
                    IRx += mode[i, j] * self.born_chgs[i,j,0]
                    IRy += mode[i, j] * self.born_chgs[i,j,1]
                    IRz += mode[i, j] * self.born_chgs[i,j,2]
            IR = IRx**2 + IRy**2 + IRz**2
            IRs.append(IR)
            if abs(IR) > 1e-2 and abs(freq) > 5e-3:
                # IR active, compute the epsilon
                epsilons.append(compute_epsilon_by_modes(mode, freq, self.born_chgs, vol, mass))
            else:
                epsilons.append(np.zeros([3,3]).flatten())
                #print('IR inactive, skip this mode')
        self.IRs = np.array(IRs)
        self.epsilons = epsilons
    def show(self):
        print("\n   Freq(cm-1)    IR Intensity     E_xx         E_yy         E_zz")
        for ir, freq, eps in zip(self.IRs, self.freqs, self.epsilons):
            print("{:12.3f} {:12.3f} {:12.3f} {:12.3f} {:12.3f}".format(freq*units.ev2cm, ir, eps[0], eps[4], eps[8]))
        eps_sum = np.sum(self.epsilons, axis=0)
        print("{:25s} {:12.3f} {:12.3f} {:12.3f}".format('Total', eps_sum[0], eps_sum[1], eps_sum[2]))
        print("{:25s} {:12.3f} {:12.3f} {:12.3f}".format('Total', eps_sum[3], eps_sum[4], eps_sum[5]))
        print("{:25s} {:12.3f} {:12.3f} {:12.3f}".format('Total', eps_sum[6], eps_sum[7], eps_sum[8]))

def compute_epsilon_by_modes(mode, freq, z, V, mass):
    """Compute the epsilon for the given mode"""
    #transform all units to hartree
    freq = freq/units.ev2hartree
    V = V/(units.a2bohr**3)  #bohr
    #mass = mass*units.proton_mass

    # compute the mode effiective charge tensors Z* (3 component)
    zt = np.zeros(3)
    for alpha in range(3):
        for beta in range(3):
            for i in range(len(z)):
                zt[alpha] += z[i, alpha, beta] * mode[i, beta] / np.sqrt(mass[i])
                #zt[alpha] += z[i, alpha, beta] * mode[i, beta]
                #zt[alpha] += z[i, alpha, beta] * mode[i, beta] / mass[i]
    epsilon = np.zeros([3,3])
    # compute the epsilon
    for alpha in range(3):
        for beta in range(3):
            epsilon[alpha, beta] = zt[alpha] * zt[beta] / ((freq)**2)
    factor = 4*units.pi/V/units.proton_mass
    epsilon *= factor
    return epsilon.flatten()


