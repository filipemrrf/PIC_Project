"""
 " @file Convergence_Test.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Compares the solutions given and checks their conversion
 " @version 3.0
 " @date 2023-05-17
 " 
 " @copyright Copyright (c) 2022
 " 
"""

from numpy import sqrt, log2
from math import pi
import sys
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
    for l in range(0, step*NPoints_low, step):
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

# Initializes the variables to hold directory and the name of the files
directory = ""
filenames = []

# Reads the command line arguments to define the files to compare to each other
for i in range(1, 6):
    filenames.append(sys.argv[i])

# Reads the command line arguments to define variables needed for the comparison
for i in range(6, len(sys.argv)):
    # Chooses the spatial step of the lower resolution solution
    if sys.argv[i] == "-SX":
        step_low = float(sys.argv[i+1])

    # Chooses the spatial step of each of the solutions
    if sys.argv[i] == "-CFL":
        cfl = float(sys.argv[i+1])

    # Chooses the number of points for the lower solution
    if sys.argv[i] == "-NP":
        NPoints_low = int(sys.argv[i+1])

    # Chooses the directory where the results of the comparison will be saved
    if sys.argv[i] == "-DIR":
        directory = sys.argv[i+1]

    # Sets the spherical flag
    if sys.argv[i] == "-S":
        if int(sys.argv[i+1]) == 1:
            spherical = True
        else:
            spherical = False


# Declares a variable to store the data of the files
Data = []

# Reads the files and saves their data
for i in range(len(filenames)):
    aux = []
    Read_File(filenames[i], aux)
    Data.append(aux)

Time = []
Norm = []

# Loops through the data
for resolution in range(4):
    # Opens the file that will be written
    FILE_point = open(directory + str(pow(2, resolution)*NPoints_low) + "p," + str(pow(2, resolution+1)*NPoints_low) + "p-Point_Comparison.dat", "w")
    
    # Declares an auxiliar variable to store the norm data temporarily
    aux = []

    # Matches the time of the solutions
    for time in range(1, len(Data[0])):
        # Writes the time in the pointwise comparison file
        if resolution == 0:
            Time.append(Data[0][time][0])

        # Calculates the relative error between solutions in the current time step
        FILE_point.write("\"Time = " + str(Data[0][time][0]) + "\n")
        Sol_Compare(Data[resolution][int(time)][1], Data[resolution+1][int(time)][1], aux, FILE_point, step_low, spherical, pow(2, resolution))


    # Saves the norm to teh appropriate list
    Norm.append(aux)

    # Closes the pointwise comparison file
    FILE_point.close()

# Loops through the norm
for resolution in range(3):
    # Declares the variable to store the log2 of the norm error comparison
    log2_NormError = []

    # Calculates the log2 of the norm error comparison and writes it to the respective file
    for i in range(len(Time)):
        log2_NormError.append(log2(Norm[resolution][i]/Norm[resolution+1][i]))

    # Plots the norm comparison 
    plt.plot(Time, log2_NormError)

# Writes the legend and saves it as a png
plt.ylim([0,5])
plt.xlabel("Time")
plt.ylabel("Log2(Norm_Comparison)")
plt.legend(["Low resolution", "Medium resolution", "High resolution"])
plt.savefig(directory + "Norm_Convergence_Comparison.png")