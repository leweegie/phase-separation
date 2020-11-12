from Cahn_Hilliard import Cahn_Hilliard
from Poisson import Poisson
import matplotlib.pyplot as plt
import sys
import plot_functions

def animation(n, nsweeps, algortithm):

    d_e = algortithm
    fig = plt.figure()
    im=plt.imshow(d_e.lattice, animated=True)

    for i in range (nsweeps):
        d_e.numerical_solution_sweep() #perform flip
        if i > (nsweeps/5):
            if (i%10 == 0):
                plt.cla()
                im=plt.imshow(d_e.lattice, animated=True)
                plt.draw()
                plt.pause(0.0001)
                print(i)

def C_H_animate():

    n       = 50
    nsweeps = 100000

    dx      = 1
    dt      = 1
    noise   = 0.1

    phi0 = 0.5

    d_e = Cahn_Hilliard(n, phi0, dx, dt)
    d_e.noise(noise)
    animation(n, nsweeps, d_e)

def C_H_free_energy():
    
    n       = 50
    nsweeps = 100000

    dx      = 1
    dt      = 1
    noise   = 0.1

    phi0 = 0.5

    d_e = Cahn_Hilliard(n, phi0, dx, dt)
    d_e.noise(noise)
    d_e.free_energy_vs_time(nsweeps, 'energy_vs_time_0point5.dat')
    plot_functions.plot_function('energy_vs_time_0point5.dat', 'energy_vs_time_0point5.png', 'Energy vs Time (phi = 0.5)', 'Time (sweeps)', 'Free Energy')

def monopole():

    n           = 50
    nsweeps     = 100000
    phi0        = 0
    accuracy    = 0.00001
    
    p = Poisson(n, phi0, accuracy)
    p.set_boundary()
    p.monopole()

    p.gauss_seidel(nsweeps)
    p.potentials_and_fields('monopole_potentials_50.dat', 'monopole_e_fields_x_50.dat', 'monopole_e_fields_y_50.dat', 'monopole_b_fields_x_50.dat', 'monopole_b_fields_y_50.dat')
    p.normalise_field('monopole_e_fields_x_50.dat', 'monopole_e_fields_y_50.dat', 'normalized_monopole_e_fields_x_50.dat', 'normalized_monopole_e_fields_y_50.dat')
    
    plot_functions.plot_heatmap('monopole_potentials_50.dat', 'monopole_potentials_50.png', 'Heatmap of potentials for monopole (z = 25)', 'x', 'y')
    plot_functions.plot_quiver('normalized_monopole_e_fields_x_50.dat', 'normalized_monopole_e_fields_y_50.dat', 'electric_field_monopole_50.png', 'E Field (z = 25)', 'x', 'y', 50)

def func_charged_wire():

    n           = 50
    nsweeps     = 100000
    phi0        = 0
    accuracy    = 0.00001
    
    p = Poisson(n, phi0, accuracy)
    p.set_boundary()
    p.charged_wire()

    p.gauss_seidel(nsweeps)
    p.potentials_and_fields('charged_wire_potentials_50.dat', 'charged_wire_e_fields_x_50.dat', 'charged_wire_e_fields_y_50.dat', 'charged_wire_b_fields_x_50.dat', 'charged_wire_b_fields_y_50.dat')
    p.normalise_field('charged_wire_b_fields_x_50.dat', 'charged_wire_b_fields_y_50.dat', 'normalised_charged_wire_b_fields_x_50.dat', 'normalised_charged_wire_b_fields_y_50.dat')
    
    plot_functions.plot_heatmap('charged_wire_potentials_50.dat', 'charged_wire_potentials_50.png', 'Heatmap of potentials for charged wire (z = 25)', 'x', 'y')
    plot_functions.plot_quiver('normalised_charged_wire_b_fields_x_50.dat', 'normalised_charged_wire_b_fields_y_50.dat', 'magnetic_field_charged_wire_50.png', 'B Field', 'x', 'y', 50)

def over_relaxation():
    n           = 50
    nsweeps     = 100000
    phi0        = 0
    accuracy    = 0.0001

    p = Poisson(n, phi0, accuracy)
    p.set_boundary()
    p.monopole()
    p.over_relaxation_data(1, 2.1, 5000)
    plot_functions.plot_function('over_relaxation.dat', 'over_relaxation.png', 'No. of Sweeps vs Omega (N = 50, accuracy = 0.0001%)', 'Omega', 'Time(Sweeps)')

#C_H_animate()
#C_H_free_energy()
#monopole()
#func_charged_wire()
#over_relaxation()