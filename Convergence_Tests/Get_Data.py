"""
 " @file Get_Data.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Uses the main executable to get data
 " @version 5.0
 " @date 2024-11-14
 " 
 " @copyright Copyright (c) 2024
 " 
"""

import os
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Process the arguments to get the data')

# Define the arguments
parser.add_argument('equation', type=str, help='The equation to analyze')
parser.add_argument('--N', type=int, default=100, help='Number of points (default: 100)')
parser.add_argument('--T', type=float, default=1, help='Time to run the simulation for (default: 1)')
parser.add_argument('--Size', type=int, default=1, help='Dimensions of the space (default: 1)')
parser.add_argument('--Runs', type=int, default=5, help='Number of different resolutions that will be run (default: 5)')
parser.add_argument('--W', type=int, default=1, help='Writing condition (default: 1)')
parser.add_argument('--Diss', type=float, default=0.02, help='Dissipation parameter (default: 0.02)')

# Parse the arguments
args = parser.parse_args()

EQ = args.equation # Chooses the equation to analyse
NPoints = args.N # Sets the initial number of points
space_size = args.Size # Sets the size of the space
Runs = args.Runs # Sets the number of different resolutions that will be run
T = args.T # Sets the time to run the simulation for
W = args.W # Sets the writing condition
diss = args.Diss # Sets the dissipation parameter

# Declaration of the variables to control the data aquisition process
if(EQ == "simple_wave"):
    Acc = 2
    Bound = "periodic"
    params = "1 1.0"

    IC = ["Sin(2pi.x)", "2pi.Cos(2pi.x)"]
    N_Ghosts = [2, 2]

    NOut = 2
    Out_Type = "solution"

    Field_Names = ["Phi", "Pi"]

if(EQ == "non_linear_simple_wave"):
    Acc = 2
    Bound = "periodic"
    params = "2 1.0 2.0"

    IC = ["Sin(2pi.x)", "2pi.Cos(2pi.x)"]
    N_Ghosts = [2, 2]

    NOut = 2
    Out_Type = "solution"

    Field_Names = ["Phi", "Pi"]


if(EQ == "spherical_wave"):
    Acc = 2
    Bound = "even_0"
    params = "1 1.0"

    IC = ["exp(-50x^2)", "(-100).exp(-50x^2)"]
    N_Ghosts = [2, 2]

    NOut = 2
    Out_Type = "solution"

    Field_Names = ["Phi", "Pi"]

if(EQ == "spherical_reduced_wave"):
    Acc = 2
    Bound = "poison"
    params = "0"

    IC = ["Gauss", "Gauss_Deriv", "0"]
    N_Ghosts = [1, 1]

    NOut = 3
    Out_Type = "solution"

    Field_Names = ["Psi", "Phi", "Pi"]

if(EQ == "non_linear_spherical_wave"):
    Acc = 2
    Bound = "even_0"
    params = "2 1.0 4.0"

    IC = ["0.01.exp(-50x^2)", "(-1).exp(-50x^2)"]
    N_Ghosts = [2, 2]

    NOut = 2
    Out_Type = "solution"

    Field_Names = ["Phi", "Pi"]

if(EQ == "adm_evolution"):
    Acc = 2
    Bound = "even_constant"
    params = "0"

    IC = ["0.002*exp(-2.0*x^2)+1", "(-0.008*x*exp(-2.0*x^2))(1+0.002*exp(-2.0*x^2))^(-1)", "0", "0.002*exp(-2.0*x^2)+1", "(-0.008*x*exp(-2.0*x^2))(1+0.002*exp(-2.0*x^2))^(-1)", "0", "0", "1", "0"]
    N_Ghosts = [2, 2]

    NOut = 9
    Out_Type = "solution"# reduction momentum hamiltonian"

    Field_Names = ["A", "DA", "KA", "B", "DB", "KB", "lambda", "alpha", "Dalpha"]

if(EQ == "compact_wave_equation"):
    Acc = 2
    Bound = "poison"
    params = "0"

    IC = ["Gauss", "Gauss_Deriv", "0", "H", "A"]
    N_Ghosts = [0, 0]

    NOut = 6
    Out_Type = "solution constraint"

    Field_Names = ["Psi", "Phi", "Pi", "H", "A"]


if(EQ == "spherical_compact_wave_equation"):
    Acc = 2
    Bound = "poison"
    params = "0"

    IC = ["Gauss", "Gauss_Deriv", "0", "H", "Omega", "L", "B"]
    N_Ghosts = [1, 1]

    NOut = 8
    Out_Type = "solution constraint"

    Field_Names = ["Psi", "Phi", "Pi", "H", "Omega", "L", "B"]

# Displays the information the aquisition process has started
os.system("echo Starting data aquisition process")

# Creates the folder to store the equation data
folder = ""
if(Acc == 2):
    folder = "Results/" + EQ + "-" + str(Acc) + "nd_order"
else:
    folder = "Results/" + EQ + "-" + str(Acc) + "th_order"

os.system("mkdir " + folder)

# Loops through all the number of points that the equation will be solved for
for j in range(0, Runs):
    # Displays the information the equation is being solved for a specific number of points
    os.system("echo Solving the equation for initial conditions with " + str(NPoints+1) + " points")

    # Creates the subfolder to store the solution with this amount of points
    current_folder = folder + "/" + str(NPoints+1) + "p"
    os.system("mkdir " + current_folder)

    # Copies the initial conditions to be used to the file
    for i in range(0, len(IC)):
        os.system("cp 'Initial_Conditions/IC-" + IC[i] + "-" + str(NPoints+1) + "p.dat' " + current_folder)

    # Writes the parameters file for the run
    # Opens the file
    FILE = open(current_folder + "/Parameters-" + EQ + "-" + str(NPoints+1) + "p.txt", "w")

    # Writes the data to the file
    FILE.write("Eq: " + EQ + "\n")
    FILE.write("Acc: " + str(Acc) + "\n")
    FILE.write("Bound: " + Bound + "\n")
    FILE.write("Params: " + params + "\n")
    FILE.write("TMax: " + str(T) + "\n")
    FILE.write("\n")
    FILE.write("IC: " + str(len(IC)))
    for i in range(0, len(IC)):
        FILE.write(" " + current_folder + "/IC-" + IC[i] + "-" + str(NPoints+1) + "p.dat")
    FILE.write("\n")
    FILE.write("NPoints: " + str(NPoints+1) + "\n")
    FILE.write("NGhosts: " + str(N_Ghosts[0] + N_Ghosts[1]) + "\n")
    FILE.write("Step_x: " + str(space_size/NPoints) + "\n")
    FILE.write("\n")
    FILE.write("Out: " + str(NOut) + " " + Out_Type + "\n")
    FILE.write("Folder: " + current_folder + "/\n")
    FILE.write("W: " + str(W) + "\n")
    FILE.write("\n")
    FILE.write("Dissipation: " + str(diss))

    # Closes the file
    FILE.close()

    # Solves the equation with the parameters given
    os.system("./main.exe " + folder + "/" + str(NPoints+1) + "p/Parameters-" + EQ + "-" + str(NPoints+1) + "p.txt")

    # Renames the output files with to the correct names of the fields
    for i in range(0, len(Field_Names)):
        os.system("mv " + current_folder + "/Field_" + str(i+1) + ".dat " + current_folder + "/" + Field_Names[i] + ".dat")

    # Updates how many points the next initial conditions will have and how many timesteps will be written to disk
    NPoints *= 2
    W *= 2

# Displays the information the aquisition process has ended
os.system("echo Data aquisition process finished")