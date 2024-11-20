"""
 " @file Self_Convergence_Test.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Compares the solutions given and checks their convergence
 " @version 4.0
 " @date 2024-11-14
 " 
 " @copyright Copyright (c) 2024
 " 
"""

import argparse
from numpy import sqrt, log2
from math import pi
import matplotlib.pyplot as plt

# Reads the file and stores its information in Data
def Read_File(filename, Data):
    IN = open(filename, "r")

    time = 0.0
    space = []

    for l in IN:
        if l[0] == "\"":
            aux = l.split()
            time = float(aux[2])
        elif l != "\n":
            aux = l.split()
            space.append((float(aux[0]), float(aux[1])))
        else:
            Data.append((time, space.copy()))
            space.clear()

    IN.close()

def Sol_Compare(low_res, high_res, norm, point_FILE, step_low, spherical, step):
    sum_error = 0.0

    # Calculates the pointwise error and the norm error for this time
    for l in range(0, step*NPoints_low+1, step):
        error = high_res[2*l][1] - low_res[l][1]
        if spherical:
            sum_error += low_res[l][0]*low_res[l][0]*error*error
        else:
            sum_error += error*error

        # Writes the pointwise error in the respective comparison file
        point_FILE.write(str(low_res[l][0]) + " " + str(error) + "\n")

    # Makes sure the pointwise comparison file is in the right format
    point_FILE.write("\n")

    # Saves the norm error for this time in the respective variable
    if spherical:
        norm.append(sqrt(4*pi*sum_error*step_low))
    else:
        norm.append(sqrt(sum_error*step_low))


# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Compare solutions and check their convergence.')

# Define the filenames argument
parser.add_argument('filenames', metavar='F', type=str, nargs='+', help='filenames to process')

# Define the optional flags
parser.add_argument('--N', type=int, default=100, help='Number of points (default: 100)')
parser.add_argument('--SX', type=float, help='Space step of the lower resolution solution')
parser.add_argument('--Runs', type=int, default=5, help='Number of different resolutions that will be run (default: 5)')
parser.add_argument('--Dir', type=str, default="", help='Directory the results will be saved to (default: "")')
parser.add_argument('--S', type=int, default=0, help='Flag to indicate if the solution is in spherical coordinates (default: 0)')
parser.add_argument('--CFL', type=float, default=0.25, help='CFL number (default: 0.25)')

# Parse the arguments
args = parser.parse_args()

# Access the parsed arguments
filenames = args.filenames # Sets the filenames to be processed
NPoints_low = args.N # Sets the number of points of the lower resolution
step_low = args.SX # Sets the space step of the lower resolution solution
Runs = args.Runs # Sets the number of different resolutions that will be run
directory = args.Dir # Sets the directory the results will be saved to
spherical = args.S # Sets the flag to indicate if the solution is in spherical coordinates
cfl = args.CFL # Sets the CFL number


# Declares a variable to store the data of the files
Data = []
aux = []

# Reads the files and saves their data
for file in filenames:
    aux.clear()
    Read_File(file, aux)
    Data.append(aux.copy())


# Declares the variables to store the time and the norm of the error
Time = []
Norm = []

# Calculates the norm of the error for each resolution
for resolution in range(Runs-1):
    # Opens the file that will be written
    FILE_point = open(directory + str(pow(2, resolution)*NPoints_low + 1) + "p," + str(pow(2, resolution+1)*NPoints_low + 1) + "p-Point_Comparison.dat", "w")
    
    # Clears the aux variable to store the norm of the error
    aux.clear()

    # Matches the time of the solutions
    for time in range(1, len(Data[0])):
        # Writes the time in the pointwise comparison file
        if resolution == 0:
            Time.append(Data[0][time][0])

        # Calculates the relative error between solutions in the current time step
        FILE_point.write("\"Time = " + str(Data[0][time][0]) + "\n")
        Sol_Compare(Data[resolution][int(time)][1], Data[resolution+1][int(time)][1], aux, FILE_point, step_low, spherical, pow(2, resolution))


    # Saves the norm to teh appropriate list
    Norm.append(aux.copy())

    # Closes the pointwise comparison file
    FILE_point.close()

legend = []

# Loops through the norm
for resolution in range(Runs-2):
    # Declares the variable to store the log2 of the norm error comparison
    log2_NormError = []

    # Calculates the log2 of the norm error comparison and writes it to the respective file
    for i in range(len(Time)):
        log2_NormError.append(log2(Norm[resolution][i]/Norm[resolution+1][i]))


    legend.append("|" + str(pow(2, resolution)*NPoints_low + 1) + "p - " + str(pow(2, resolution+1)*NPoints_low + 1) + "p| / |" + str(pow(2, resolution+1)*NPoints_low + 1) + "p - " + str(pow(2, resolution+2)*NPoints_low + 1) + "p|")

    # Plots the norm comparison 
    plt.plot(Time, log2_NormError)

# Writes the legend and saves it as a png
plt.ylim([0,4])
plt.xlabel("Time")
plt.ylabel("Log2(Norm_Comparison)")
plt.legend(legend)
plt.savefig(directory + "Norm_Convergence_Comparison.png")