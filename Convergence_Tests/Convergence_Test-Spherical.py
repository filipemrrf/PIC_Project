"""
 " @file Convergence_Test.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Compares the solutions given and checks their conversion
 " @version 1.0
 " @date 2023-02-16
 " 
 " @copyright Copyright (c) 2023
 " 
"""

from numpy import sqrt, log2
from math import pi
import sys
import matplotlib.pyplot as plt

# Initializes the directory to a default
directory = ""

# Reads the command line arguments to define the files to compare to each other
input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]

# Scale factor declaration (for error comparison between the files)
S = 1

# Reads the command line arguments to define variables needed for the comparison
for i in range(4, len(sys.argv)):
    # Chooses the spatial step of each of the solutions
    if sys.argv[i] == "-SX":
        step_x1 = float(sys.argv[i+1])
        step_x2 = float(sys.argv[i+2])
        step_x3 = float(sys.argv[i+3])

    # Chooses the number of points for the files
    if sys.argv[i] == "-NP":
        NPoints1 = int(sys.argv[i+1])
        NPoints2 = int(sys.argv[i+2])
        NPoints3 = int(sys.argv[i+3])
        Ratio21 = int(NPoints2/NPoints1)
        Ratio31 = int(NPoints3/NPoints1)

    # Chooses the directory where the results of the comparison will be saved
    if sys.argv[i] == "-DIR":
        directory = sys.argv[i+1]

    # Sets the comparison constant (to scale for pointwise convergence check)
    if sys.argv[i] == "-S":
        S = int(sys.argv[i+1])


# Opens the file that will be written
FILE_point12 = open(directory + str(NPoints1) + "p," + str(NPoints2) + "p-Point_Comparison.dat", "w")
FILE_point23 = open(directory + str(NPoints2) + "p," + str(NPoints3) + "p-Point_Comparison.dat", "w")
FILE_norm = open(directory + str(NPoints1) + "p," + str(NPoints3) + "p-Norm_Comparison.dat", "w")

# Opens the file that will be compared
IN1 = open(input1, "r")
IN2 = open(input2, "r")
IN3 = open(input3, "r")

# Declares the variables that will store the solutions
time = 0.0
space = []

Data1 = []
Data2 = []
Data3 = []

# Reads the 1st input file and saves its data
for l in IN1:
    if l[0] == "\"":
        aux = l.split()
        time = float(aux[2])
    elif l != "\n":
        aux = l.split()
        space.append(float(aux[1]))
    else:
        auxl = space.copy()
        Data1.append((time, auxl))
        space.clear()

# Reads the 2nd input file and saves its data
for l in IN2:
    if l[0] == "\"":
        aux = l.split()
        time = float(aux[2])
    elif l != "\n":
        aux = l.split()
        space.append(float(aux[1]))
    else:
        auxl = space.copy()
        Data2.append((time, auxl))
        space.clear()

# Reads the 3rd input file and saves its data
for l in IN3:
    if l[0] == "\"":
        aux = l.split()
        time = float(aux[2])
    elif l != "\n":
        aux = l.split()
        space.append(float(aux[1]))
    else:
        Data3.append((time, list(space)))
        space.clear()


# Declares the variables needed to compare the first 2 solutions solutions
time = []
norm12 = []
norm23 = []
error = 0.0
sum_error = 0.0

# Matches the time of the solutions
for i in range(1, len(Data1)):
    # Declaration of the variable that tells if the matching indexes for the time were found
    found = False

    for j in range(1, len(Data2)):
        for k in range(1, len(Data3)):
            if Data1[i][0] == Data2[j][0] == Data3[k][0]:
                # Writes the time in the pointwise comparison file
                time.append(Data2[j][0])


                # Compares the 1st and the 2nd solution
                sum_error = 0.0
                FILE_point12.write("\"Time = " + str(Data2[j][0]) + "\n")

                # Calculates the pointwise error and the norm error for this time
                for l in range(NPoints1):
                    error = Data2[j][1][Ratio21*l] - Data1[i][1][l]
                    sum_error += (l*step_x1)*(l*step_x1)*error*error

                    # Writes the pointwise error in the respective comparison file
                    FILE_point12.write(str(l*step_x1) + " " + str(error) + "\n")

                # Makes sure the pointwise comparison file is in the right format
                FILE_point12.write("\n")

                # Saves the norm error for this time in the respective variable
                norm12.append(sqrt(4*pi*sum_error*step_x1))


                # Compares the 2nd and the 3rd solution
                sum_error = 0.0
                FILE_point23.write("\"Time = " + str(Data2[i][0]) + "\n")

                # Calculates the pointwise error and the norm error for this time
                for l in range(0, NPoints1):
                    error = Data3[k][1][Ratio31*l] - Data2[j][1][Ratio21*l]
                    sum_error += (l*step_x1)*(l*step_x1)*error*error
                    error *= S

                    # Writes the pointwise error in the respective comparison file
                    FILE_point23.write(str(l*step_x1) + " " + str(error) + "\n")

                # Makes sure the pointwise comparison file is in the right format
                FILE_point23.write("\n")

                # Saves the norm error for this time in the respective variable
                norm23.append(sqrt(4*pi*sum_error*step_x1))


                # Sets the variable saying the matching indexes were found
                found = True
                break

            # If the maching indexes were found, jump straight to the next time
            if found:
                break


# Declares the variable to store the log2 of the norm error comparison
log2_NormError = []

# Writes the header of the norm comparison file
FILE_norm.write("#Time log2_NormError\n")

# Calculates the log2 of the norm error comparison and writes it to the respective file
for i in range(len(time)):
    aux = log2(norm12[i]/norm23[i])
    log2_NormError.append(aux)
    FILE_norm.write(str(time[i]) + " " + str(aux) + "\n")

# Plots the norm comparison and saves it as a png
plt.plot(time, log2_NormError)
plt.title("Norm Convergence")
plt.ylim([0,5])
plt.xlabel("Time")
plt.ylabel("Log2(Norm_Comparison(" + str(NPoints1) + "p," + str(NPoints2) + "p," + str(NPoints3) + "p))")
plt.savefig(directory + str(NPoints1) + "p," + str(NPoints2) + "p," + str(NPoints3) + "p-Norm_Convergence_Comparison.png")


# Closes the files
FILE_point12.close()
FILE_point23.close()
FILE_norm.close()