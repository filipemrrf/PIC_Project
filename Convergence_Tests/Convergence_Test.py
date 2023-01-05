import numpy
import sys
import matplotlib.pyplot as plt

# Initializes the filename and the number of points to a default
filename = ""

# Reads the command line arguments to define the 2 files to compare to the analytical solution (and with each other)
input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]

# Multiplicative factor declaration (for error comparison between both files)
E12 = 1
E23 = 1

# Reads the command line arguments to define variables needed for the comparison
for i in range(4, len(sys.argv)):
    # Chooses the spatial step of each of the solutions
    if sys.argv[i] == "-SX":
        step_x1 = float(sys.argv[i+1])
        step_x2 = float(sys.argv[i+2])
        step_x3 = float(sys.argv[i+3])

    # Chooses the number of points for both files
    if sys.argv[i] == "-NP":
        NPoints1 = int(sys.argv[i+1])
        NPoints2 = int(sys.argv[i+2])
        NPoints3 = int(sys.argv[i+3])

    # Chooses the filename of the files that save the results of the comparison
    if sys.argv[i] == "-FN":
        filename = sys.argv[i+1]

    # Sets the comparison constant (that will be multiplied to the right)
    if sys.argv[i] == "-E":
        E12 = int(sys.argv[i+1])
        E23 = int(sys.argv[i+2])


# Opens the file that will be written
FILE_point12 = open(filename + str(NPoints1) + "p_" + str(NPoints2) + "p-Point_Comparison.dat", "w")
FILE_point23 = open(filename + str(NPoints2) + "p_" + str(NPoints3) + "p-Point_Comparison.dat", "w")
FILE_norm = open(filename + str(NPoints1) + "p_" + str(NPoints2) + "p-Norm_Comparison.dat", "w")

# Opens the file that will be compared
IN1 = open(input1, "r")
IN2 = open(input2, "r")
IN3 = open(input3, "r")

# Reads the input data and saves it to lists
time = 0.0
space = []

Data1 = []
Data2 = []
Data3 = []

for l in IN1:
    if l[0] == "\"":
        aux = l.split()
        time = float(aux[2])
    elif l != "\n":
        aux = l.split()
        space.append(float(aux[1]))
    else:
        Data1.append((time, list(space)))
        space.clear()


for l in IN2:
    if l[0] == "\"":
        aux = l.split()
        time = float(aux[2])
    elif l != "\n":
        aux = l.split()
        space.append(float(aux[1]))
    else:
        Data2.append((time, list(space)))
        space.clear()

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


# Compares the solutions 1 and 2
time12 = []
norm12 = []
error = 0.0
sum_error = 0.0

for i in range(1, len(Data2)):
    for j in range(1, len(Data1)):
        if Data2[i][0] == Data1[j][0]:
            time12.append(Data2[i][0])
            sum_error = 0.0
            FILE_point12.write("\"Time = " + str(Data2[i][0]) + "\n")

            for k in range(NPoints1):
                error = Data2[i][1][2*k] - Data1[j][1][k]
                sum_error += error*error
                error *= E12

                FILE_point12.write(str(k*step_x1) + " " + str(error) + "\n")

            FILE_point12.write("\n")
            norm12.append(numpy.sqrt(sum_error*step_x1))
            break


# Compares the solutions 2 and 3
time23 = []
norm23 = []

for i in range(1, len(Data2)):
    for j in range(1, len(Data3)):
        if Data2[i][0] == Data3[j][0]:
            time23.append(Data2[i][0])
            sum_error = 0.0
            FILE_point23.write("\"Time = " + str(Data2[i][0]) + "\n")

            for k in range(NPoints1):
                error = Data3[j][1][4*k] - Data2[i][1][2*k]
                sum_error += error*error
                error *= E23

                FILE_point23.write(str(k*step_x1) + " " + str(error) + "\n")

            FILE_point23.write("\n")
            norm23.append(numpy.sqrt(sum_error*step_x1))
            break


# Divides the norm comparison of both sets of values
time = []
log2_NormError = []

FILE_norm.write("#Time log2_NormError\n")

for i in range(len(norm23)):
    for j in range(len(norm12)):
        if time23[i] == time12[j]:
            time.append(time12[j])
            log2_NormError.append(numpy.log2(norm12[j]/norm23[i]))

            FILE_norm.write(str(time12[j]) + " " + str(numpy.log2(norm12[j]/norm23[i])) + "\n")

# Plots the norm comparison and saves it as a png
plt.plot(time, log2_NormError)
plt.ylim([0,5])
plt.xlabel("Time")
plt.ylabel("Log2(Norm_Error-" + str(NPoints1) + "p/Norm_Error-" + str(NPoints3) + "p)")
plt.savefig("Convergence/Norm_Convergence_Comparison.png")


# Closes the files
FILE_point12.close()
FILE_point23.close()
FILE_norm.close()