"""
 " @file Get_Data.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Uses the main executable to get data
 " @version 3.0
 " @date 2023-04-01
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import os

# Declaration of the variables to control the data aquisition process
EQ = "spherical_wave"
Acc = 4
Bound = "even_0"
params = "1 1.0"
T = 3
IC = ["exp(-50x^2)", "(-100).exp(-50x^2)"]
NOut = 2
Out_Type = "solution"
W = [1, 1]
diss = 0.0

# Displays the information the aquisition process has started
os.system("echo Starting data aquisition process")

# Creates the folder to store the equation data
folder = ""
if(Acc == 2):
    folder = "Results/" + EQ + "-" + str(Acc) + "nd_order"
else:
    folder = "Results/" + EQ + "-" + str(Acc) + "th_order"

os.system("mkdir " + folder)

# Initializes the number of points
NPoints = 250

# Loops through all the number of points that the equation will be solved for
for i in range(0, 5):
    # Displays the information the equation is being solved for a specific number of points
    os.system("echo Solving the equation for initial conditions with " + str(NPoints) + " points")

    # Creates the subfolder to store the solution with this amount of points
    os.system("mkdir " + folder + "/" + str(NPoints) + "p")

    # Copies the initial conditions to be used to the file
    for i in range(0, len(IC)):
        os.system("cp 'Initial_Conditions/IC-" + IC[i] + "-" + str(NPoints) + "p.dat' " + folder + "/" + str(NPoints) + "p")

    # Writes the parameters file for the run
    # Opens the file
    FILE = open(folder + "/" + str(NPoints) + "p/Parameters-" + EQ + "-" + str(NPoints) + "p.txt", "w")

    # Writes the data to the file
    FILE.write("Eq: " + EQ + "\n")
    FILE.write("Acc: " + str(Acc) + "\n")
    FILE.write("Bound: " + Bound + "\n")
    FILE.write("params: " + params + "\n")
    FILE.write("tmax: " + str(T) + "\n")
    FILE.write("\n")
    FILE.write("IC: ")
    for i in range(0, len(IC)):
        FILE.write(folder + "/" + str(NPoints) + "p/IC-" + IC[i] + "-" + str(NPoints) + "p.dat ")
    FILE.write("\n")
    FILE.write("NPoints: " + str(NPoints) + "\n")
    FILE.write("step_x: " + str(5/NPoints) + "\n")
    FILE.write("\n")
    FILE.write("NOut: " + str(NOut) + "\n")
    FILE.write("Out_Type: " + Out_Type + "\n")
    FILE.write("Out_Filename: " + folder + "/" + str(NPoints) + "p/" + EQ + "_Phi-" + str(NPoints) + "p.dat " + folder + "/" + str(NPoints) + "p/" + EQ + "_Pi-" + str(NPoints) + "p.dat\n")
    FILE.write("W: " + str(W[0]) + " " + str(W[1]) + "\n")
    FILE.write("\n")
    FILE.write("dissipation: " + str(diss))

    # Closes the file
    FILE.close()

    # Solves the equation with the parameters given
    os.system("./main " + folder + "/" + str(NPoints) + "p/Parameters-" + EQ + "-" + str(NPoints) + "p.txt")

    os.system("echo \n")

    # Updates how many points the next initial conditions will have and how many timesteps will be written to disk
    NPoints *= 2
    W[0] *= 2
    W[1] *= 2

# Displays the information the aquisition process has ended
os.system("echo Data aquisition process finished")