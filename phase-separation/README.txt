###

The Code:

- All code can be run from main file.
- Three files:  main.py(Four functions set up to complete tasks specified in CP, by executing functions in other python files), 
                Cahn_Hilliard.py(contains all code used to complete 3.1 on CP)
                Poisson.py(contains all code used to complete 3.2 on CP)
- All variables passed into classes can be changed in respective functions in main.py file

###

main.py:

C_H_animate() dsiplays animation of oil/water emulsion.
C_H_free_energy() calculates free energy vs time graphs.

monopole() calculates all data for monopole questions and displays appropriate figures
func_charged_wire() calculates all data for charged wire question and displays appropriate figures
over_relaxation() calculates time to convergence for list of omegas

###

Checklist and respective filenames::::

3.1:

(Phi = 0)
-energy_vs_time_0.dat
-energy_vs_time_0.png

(Phi = 0.5)
-energy_vs_time_0point5.dat
-energy_vs_time_0point5.png

3.2:

(monopole)
-monopole_potentials_50.dat                 (2D 50x50 array of potential slice @ z = 25)
-monopole_potentials_50.png                 (heatmap of potential slice @ z = 25)

-monopole_e_fields_x_50.dat                 (2D 50 x 50 array of x components of electric field @ z = 25)
-monopole_e_fields_y_50.dat                 (2D 50 x 50 array of y components of electric field @ z = 25)
-normalized_monopole_e_fields_x_50.dat      (2D 50 x 50 array of normalised x components of electric field @ z = 25 for plotting)
-normalized_monopole_e_fields_y_50.dat      (2D 50 x 50 array of normalised y components of electric field @ z = 25 for plotting)
-electric_field_monopole_50.png             (quiver plot of electric field slice @ z = 25)

(charged wire)
-charged_wire_potentials_50.dat             (2D 50x50 array of potential slice @ z = 25)
-charged_wire_potentials_50.png             (heatmap of potential slice @ z = 25)

-charged_wire_b_fields_x_50.dat             (2D 50 x 50 array of x components of magnetic field @ z = 25)
-charged_wire_b_fields_y_50.dat             (2D 50 x 50 array of y components of magnetic field @ z = 25)
-normalised_charged_wire_b_fields_x_50.dat  (2D 50 x 50 array of normalised x components of magnetic field @ z = 25 for plotting)
-normalised_charged_wire_b_fields_y_50.dat  (2D 50 x 50 array of normalised y components of magnetic field @ z = 25 for plotting)
-magnetic_field_charged_wire_50.png         (quiver plot of magnetic field slice @ z = 25)

(SOR)
-over_relaxation.dat                        (numbers of sweeps to convergence for list of omegas)
-over_relaxation.png                        (Number of Sweeps vs omega)