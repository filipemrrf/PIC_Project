"""
 " @file Sol_Comparer-Sin.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Compares the data to the function sin(x+t)
 " @version 2.3
 " @date 2023-01-13
 " 
 " @copyright Copyright (c) 2023
 " 
"""

from numpy import sin, pi, sqrt, log2
import sys
import matplotlib.pyplot as plt

# Initializes the directory to a default
directory = ""

# Reads the command line arguments to define the 2 files to compare to the analytical solution (and with each other)
input1 = sys.argv[1]
input2 = sys.argv[2]

# Scale factor declaration (for error comparison between both files)
S = 1

# Reads the command line arguments to define variables needed for the comparison
for i in range(3, len(sys.argv)):
    # Chooses the spatial step of each of the solutions
    if sys.argv[i] == "-SX":
        step_x1 = float(sys.argv[i+1])
        step_x2 = float(sys.argv[i+2])

    # Chooses the number of points for both files
    if sys.argv[i] == "-NP":
        NPoints1 = int(sys.argv[i+1])
        NPoints2 = int(sys.argv[i+2])
        Ratio = int(NPoints2/NPoints1)

    # Chooses the directory of the files that save the results of the comparison
    if sys.argv[i] == "-DIR":
        directory = sys.argv[i+1]

    # Sets the comparison constant (that will be multiplied to the right)
    if sys.argv[i] == "-S":
        S = int(sys.argv[i+1])

# Opens the file that will be written
FILE_point1 = open(directory + str(NPoints1) + "p-Exact_Sol-Point_Error.dat", "w")
FILE_point2 = open(directory + str(NPoints2) + "p-Exact_Sol-Point_Error.dat", "w")
FILE_norm = open(directory + str(NPoints1) + "p," + str(NPoints2) + "p-Exact_Sol_Norm_Comparison.dat", "w")

# Opens the file that will be compared
IN1 = open(input1, "r")
IN2 = open(input2, "r")

# Declares the variables that will store the solutions
time = 0.0
space = []

Data1 = []
Data2 = []

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


# Declares the variables needed to compare the first 2 solutions solutions
time = []
norm1 = []
norm2 = []
error = 0.0
sum_error = 0.0

# Matches the time of the solutions
for i in range(1, len(Data1)):
    for j in range(1, len(Data2)):
        if Data1[i][0] == Data2[j][0]:
            # Writes the time in the pointwise comparison file
            time.append(Data1[i][0])


            # Compares the 1st solution to the analytical solution
            sum_error = 0.0
            FILE_point1.write("\"Time = " + str(Data2[j][0]) + "\n")

            # Calculates the pointwise error and the norm error for this time
            for l in range(NPoints1):
                error = sin(2*pi*(l*step_x1 + Data1[i][0])) - Data1[i][1][l]
                sum_error += error*error

                # Writes the pointwise error in the respective comparison file
                FILE_point1.write(str(l*step_x1) + " " + str(error) + "\n")

            # Makes sure the pointwise comparison file is in the right format
            FILE_point1.write("\n")

            # Saves the norm error for this time in the respective variable
            norm1.append(sqrt(sum_error*step_x1))


            # Compares the 2nd solution to the analytical solution
            sum_error = 0.0
            FILE_point2.write("\"Time = " + str(Data2[i][0]) + "\n")

            # Calculates the pointwise error and the norm error for this time
            for l in range(0, NPoints1):
                error = sin(2*pi*(l*step_x1 + Data2[j][0])) - Data2[j][1][l*Ratio]
                sum_error += error*error
                error *= S

                # Writes the pointwise error in the respective comparison file
                FILE_point2.write(str(l*step_x1) + " " + str(error) + "\n")

            # Makes sure the pointwise comparison file is in the right format
            FILE_point2.write("\n")

            # Saves the norm error for this time in the respective variable
            norm2.append(sqrt(sum_error*step_x1))

            # If the maching indexes were found, jump straight to the next time
            break

# Declares the variable to store the log2 of the norm error comparison
log2_NormError = []

# Writes the header of the norm comparison file
FILE_norm.write("#Time log2_NormError\n")

# Calculates the log2 of the norm error comparison and writes it to the respective file
for i in range(len(time)):
    aux = log2(norm1[i]/norm2[i])
    log2_NormError.append(aux)
    FILE_norm.write(str(time[i]) + " " + str(aux) + "\n")

# Plots the norm comparison and saves it as a png
plt.plot(time, log2_NormError)
plt.title("Norm Convergence")
plt.ylim([0,5])
plt.xlabel("Time")
plt.ylabel("Log2(Norm_Comparison(" + str(NPoints1) + "p,Exact_Sol," + str(NPoints2) + "p))")
plt.savefig(directory + str(NPoints1) + "p,Exact_Sol," + str(NPoints2) + "p-Norm_Convergence_Comparison.png")


# Closes the files
FILE_point1.close()
FILE_point2.close()
FILE_norm.close()