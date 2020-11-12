import numpy as np
import scipy
import math
from numpy.random import rand
import matplotlib.pyplot as plt

class Poisson(object):

    def __init__(self, n, phi0, accuracy):
        self.n = n
        self.phi0 = phi0
        self.lattice = np.zeros((n,n,n))
        self.rho_lattice = np.full((n, n, n), 0)
        self.dx = 1
        self.dt = (self.dx * self.dx)/ 6
        self.accuracy = accuracy

    def set_boundary(self):
        self.lattice[0, :, :]           = 0
        self.lattice[:, 0, :]           = 0
        self.lattice[:, :, 0]           = 0
        self.lattice[self.n-1, :, :]    = 0
        self.lattice[:, self.n-1, :]    = 0
        self.lattice[:, :, self.n-1]    = 0

    def monopole(self):
        mid = self.n/2
        self.rho_lattice[mid, mid, mid] = 1

    def charged_wire(self):
        mid = self.n/2
        self.rho_lattice[mid, mid, :] = 1

    #function to compute sum of nearest nieghbours on 3-D lattice
    def neighbour_sum(self, i, j, k):
        nb = self.lattice[(i+1)%self.n, j, k] + self.lattice[(i-1)%self.n, j, k] + self.lattice[i, (j+1)%self.n, k] + self.lattice[i, (j-1)%self.n, k] + self.lattice[i,j,(k+1)%self.n] + self.lattice[i, j, (k-1)%self.n]
        return nb

    #function to compute the numerical solution of the poisson equation using jacobi algorithm
    def jacobi(self, nsweeps):

        for x in range(nsweeps):

            n_plus_1 = np.zeros([self.n, self.n, self.n])

            for i in range (1, self.n-1):
                for j in range (1, self.n-1):
                    for k in range (1, self.n-1):

                        nb = self.neighbour_sum(i,j,k)
                        rho = self.rho_lattice[i,j,k]

                        n_plus_1[i,j,k] = (nb + rho) / 6

            lattice_before = np.sum(self.lattice)
            self.lattice = np.copy(n_plus_1)
            lattice_sum = np.sum(self.lattice)
            convergence = abs(lattice_sum - lattice_before)
            if convergence < (lattice_sum*self.accuracy):
                break

    #function to compute the numerical solution of the poisson equation using gauss-seidel algorithm
    def gauss_seidel(self, nsweeps):

        for x in range(nsweeps):

            if (x % 10 == 0):
                print(x)

            lattice_before = np.sum(self.lattice)

            for i in range (1, self.n-1):
                for j in range (1, self.n-1):
                    for k in range (1, self.n-1):
                        
                        nb = self.neighbour_sum(i,j,k)
                        rho = self.rho_lattice[i,j,k]
                        phi = (nb + rho) / 6

                        self.lattice[i,j,k] = phi

            lattice_sum = np.sum(self.lattice)
            convergence = abs(lattice_sum - lattice_before)
            if convergence < (lattice_sum*self.accuracy):
                break

    
    def over_relaxation(self, nsweeps, omega):

        for x in range(nsweeps):

            if (x % 10 == 0):
                print(x)

            lattice_before = np.sum(self.lattice)

            for i in range (1, self.n-1):
                for j in range (1, self.n-1):
                    for k in range (1, self.n-1):

                        phi_before = self.lattice[i,j,k]
                        
                        nb = self.neighbour_sum(i,j,k)
                        rho = self.rho_lattice[i,j,k]
                        phi_after = (nb + rho) / 6

                        

                        self.lattice[i,j,k] = ((1-omega)* phi_before) + (omega * phi_after)

            lattice_sum = np.sum(self.lattice)
            convergence = abs(lattice_sum - lattice_before)
            if convergence < (lattice_sum*self.accuracy):
                break

        return x

    def over_relaxation_data(self,omega_min, omega_max, max_sweeps):
        omega = np.linspace(omega_min,omega_max,12)
        sweeps = []
        for i in range (len(omega)):
            self.lattice = np.zeros((self.n,self.n,self.n))
            num_sweeps = self.over_relaxation(max_sweeps, omega[i])
            sweeps.append(num_sweeps)

        np.savetxt('over_relaxation.dat', np.column_stack([omega, sweeps]))

    def get_potentials(self):
        k = self.n / 2
        
        potentials = self.lattice[:,:,k]

        return potentials

    def get_e_fields(self):
        k = self.n/2
        xgrad, ygrad = np.gradient(self.lattice[:, :, k])
        return xgrad, ygrad

    def get_b_fields(self):
        Bx = np.zeros((self.n, self.n))
        By = np.zeros((self.n, self.n))
        potential_slice = self.get_potentials()
        for i in range(self.n):
            for j in range(self.n):
                Bx[i,j] = (potential_slice[i, (j+1)%self.n] - potential_slice[i, (j-1)%self.n]) / (2 * self.dx)
                By[i,j] = -1 * (potential_slice[(i+1)%self.n, j] - potential_slice[(i-1)%self.n, j]) / (2 * self.dx)
        
        return Bx, By

    def potentials_and_fields(self, potential_filename, ex_filename, ey_filename, bx_filename, by_filename):
        
        k = self.n/2

        p = self.get_potentials()
        xgrad, ygrad = self.get_e_fields()
        Bx, By = self.get_b_fields()

        np.savetxt(potential_filename, p)
        np.savetxt(ex_filename, xgrad)
        np.savetxt(ey_filename, ygrad)
        np.savetxt(bx_filename, Bx)
        np.savetxt(by_filename, By)

    def normalise_field(self, x_filename, y_filename, x_savefile, y_savefile):

        x = np.loadtxt(x_filename)
        y = np.loadtxt(y_filename)

        for i in range(self.n):
            for j in range(self.n):
                normal = np.sqrt((x[i,j] * x[i,j])+(y[i,j] * y[i,j]))
                x[i,j] = x[i,j]/normal
                y[i,j] = y[i,j]/normal

        np.savetxt(x_savefile, x)
        np.savetxt(y_savefile, y)