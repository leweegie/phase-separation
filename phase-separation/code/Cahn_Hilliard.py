import numpy as np
import scipy
import math
from numpy.random import rand
import matplotlib.pyplot as plt

class Cahn_Hilliard(object):

    def __init__(self, n, phi0, dx, dt):
        self.n = n
        self.lattice = np.full((n, n), phi0)
        self.phi0 = phi0
        self.dx = dx
        self.dt = dt
        self.M = 0.1
        self.a = 0.1
        self.b = 0.1
        self.kappa = 0.1

    #add noise to lattice
    def noise(self, noise):
        self.lattice = np.random.normal(self.phi0, noise, (self.n, self.n))

    #function to compute the numerical solution of the differential equation using Euler  algorithm
    def numerical_solution_sweep(self):

        n_plus_1 = np.zeros([self.n, self.n])
        for i in range (self.n):
            for j in range (self.n):

                chemical_potentials = self.cp((i + 1)%self.n, j) + self.cp((i - 1)%self.n, j) + self.cp(i, (j + 1)%self.n) + self.cp(i, (j - 1)%self.n) - (4 * self.cp(i, j))
                diff_constants = (self.M * self.dt)/(self.dx * self.dx)

                n_plus_1[i,j] = self.lattice[i,j] + (diff_constants * chemical_potentials) 

        self.lattice = np.copy(n_plus_1)
        np.clip(self.lattice, -1, 1)


    #function to compute discretised chemical potential
    def cp(self, i, j):

        compositional_orders = self.lattice[(i + 1)%self.n, j] + self.lattice[(i - 1)%self.n, j] + self.lattice[i, (j + 1)%self.n] + self.lattice[i, (j - 1)%self.n] - (4 * self.lattice[i, j])
        c_o_cubed = self.lattice[i,j] * self.lattice[i,j] * self.lattice[i, j]
        diff_constants = self.kappa / (self.dx * self.dx)
        
        chemical_potential = (-self.a * self.lattice[i,j]) + (self.b * c_o_cubed) + (-diff_constants * compositional_orders)

        return chemical_potential

    #function to compute free energy of system
    def free_energy(self):

        f_e_d = []

        for i in range (self.n):
            for j in range (self.n):
                xgrad, ygrad = np.gradient(self.lattice)
                grad = np.sqrt((xgrad*xgrad) + (ygrad*ygrad))
                phi = self.lattice[i,j]
                f_e_d.append((-(self.a/2)*(phi * phi)) + ((self.a/4)*(phi * phi * phi * phi)) + ((self.kappa/2) * ((grad) * (grad))))
        
        f_e = np.sum(f_e_d) / (self.n * self.n)
        return f_e

    #function to compute free enerfy with respect to time
    def free_energy_vs_time(self, nsweeps, filename):

        times = []
        free_energies = []

        for i in range (nsweeps):
            if (i % 500 == 0):
                free_energy = self.free_energy()
                times.append(i)
                free_energies.append(free_energy)
                print(i)
            self.numerical_solution_sweep()
 
        np.savetxt(filename, np.column_stack([times, free_energies]))