"""
 " @file Get_Data.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Uses the main executable to get data
 " @version 4.0
 " @date 2023-04-01
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import os

# Chooses the equation to analyse
EQ = ""
Acc = 0

# Declaration of the variables to control the data aquisition process
if(EQ == "simple_wave"):
    Bound = "periodic"
    params = "1 1.0"
    T = 5
    IC = ["Sin(2pi.x)", "2pi.Cos(2pi.x)"]
    NOut = 2
    Out_Type = "solution"
    Out_Filename = ["Phi", "Pi"]
    W = [1, 1]
    diss = 0.0

    space_size = 1

if(EQ == "non_linear_simple_wave"):
    Bound = "periodic"
    params = "2 1.0 2.0"
    T = 5
    IC = ["Sin(2pi.x)", "2pi.Cos(2pi.x)"]
    NOut = 2
    Out_Type = "solution"
    Out_Filename = ["Phi", "Pi"]
    W = [1, 1]

    if(Acc == 2):
        diss = 0.0
    if(Acc == 4):
        diss = 0.00

    space_size = 1

if(EQ == "spherical_wave"):
    Bound = "even_0"
    params = "1 1.0"
    T = 3
    IC = ["exp(-50x^2)", "(-100).exp(-50x^2)"]
    NOut = 2
    Out_Type = "solution"
    Out_Filename = ["Phi", "Pi"]
    W = [1, 1]

    if(Acc == 2):
        diss = 0.02
    if(Acc == 4):
        diss = 0.007

    space_size = 5

if(EQ == "non_linear_spherical_wave"):
    Bound = "even_0"
    params = "2 1.0 4.0"
    T = 3
    IC = ["0.01.exp(-50x^2)", "(-1).exp(-50x^2)"]
    NOut = 2
    Out_Type = "solution"
    Out_Filename = ["Phi", "Pi"]
    W = [1, 1]

    if(Acc == 2):
        diss = 0.02
    if(Acc == 4):
        diss = 0.007

    space_size = 5

if(EQ == "adm_evolution"):
    Bound = "even_constant"
    params = "0"
    T = 5

    IC = ["0.002*exp(-2.0*x^2)+1", "(-0.008*x*exp(-2.0*x^2))(1+0.002*exp(-2.0*x^2))^(-1)", "0", "0.002*exp(-2.0*x^2)+1", "(-0.008*x*exp(-2.0*x^2))(1+0.002*exp(-2.0*x^2))^(-1)", "0", "0", "1", "0"]
    NOut = 13
    Out_Type = "solution reduction momentum hamiltonian"
    Out_Filename = ["A", "DA", "KA", "B", "DB", "KB", "lambda", "alpha", "Dalpha", "Reduction_A", "Reduction_B", "Momentum", "Hamiltonian"]
    W = [1, 1, 1, 1, 1, 1, 1, 1, 1, 10, 10, 10, 10]

    Acc = 2
    diss = 0.02

    space_size = 10

# Sets the length of the system and the point density
point_density = 50

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
NPoints = space_size*point_density

# Loops through all the number of points that the equation will be solved for
for j in range(0, 5):
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
    FILE.write("step_x: " + str(space_size/NPoints) + "\n")
    FILE.write("\n")
    FILE.write("NOut: " + str(NOut) + "\n")
    FILE.write("Out_Type: " + Out_Type + "\n")
    FILE.write("Out_Filename: ") 
    for i in range(0, len(Out_Filename)):
        FILE.write(folder + "/" + str(NPoints) + "p/" + EQ + "_" + Out_Filename[i] + "-" + str(NPoints) + "p.dat ")
    FILE.write("\n")
    FILE.write("W: ")
    for i in range(0, len(W)):
        FILE.write(str(W[0]) + " ")
    FILE.write("\n\n")
    FILE.write("dissipation: " + str(diss))

    # Closes the file
    FILE.close()

    # Solves the equation with the parameters given
    os.system("./main " + folder + "/" + str(NPoints) + "p/Parameters-" + EQ + "-" + str(NPoints) + "p.txt")

    # Updates how many points the next initial conditions will have and how many timesteps will be written to disk
    NPoints *= 2
    for i in range(0, len(W)):
        W[i] *= 2

# Displays the information the aquisition process has ended
os.system("echo Data aquisition process finished")