import numpy
import sys
import matplotlib.pyplot as plt

# Initializes the filename and the number of points to a default
filename = "Error"

# Reads the command line arguments to define the 2 files to compare to the analytical solution (and with each other)
input1 = sys.argv[1]
input2 = sys.argv[2]

# Multiplicative factor declaration (for error comparison between both files)
E = 1

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

    # Chooses the filename of the files that save the results of the comparison
    if sys.argv[i] == "-FN":
        filename = sys.argv[i+1]

    # Sets the comparison constant (that will be multiplied to the right)
    if sys.argv[i] == "-E":
        E = int(sys.argv[i+1])

# Opens the file that will be written
FILE_point1 = open(filename + str(NPoints1) + "p-Point_Error.dat", "w")
FILE_point2 = open(filename + str(NPoints2) + "p-Point_Error.dat", "w")
FILE_norm = open(filename + "-Norm_Error.dat", "w")

# Opens the file that will be compared
IN1 = open(input1, "r")
IN2 = open(input2, "r")

# Writes in the output file the error of the parameter file compared to the analitical solution 
time = 0.0
sum_error = 0.0
Norm_error1 = []
Norm_error2 = []

for l in IN1:
    if l[0] == "\"":
        aux = l.split()
        time = float(aux[2])
        FILE_point1.write(l)
        sum_error = 0.0
    elif l != "\n":
        aux = l.split()
        error = numpy.sin(2*numpy.pi*(float(aux[0]) + time)) - float(aux[1])
        sum_error += error*error
        FILE_point1.write(aux[0] + " " + str(error) + "\n")
    else:
        FILE_point1.write("\n")
        sum_error = numpy.sqrt(sum_error*step_x1)
        Norm_error1.append((time, sum_error))


for l in IN2:
    if l[0] == "\"":
        aux = l.split()
        time = float(aux[2])
        FILE_point2.write(l)
        sum_error = 0.0
    elif l != "\n":
        aux = l.split()
        error = numpy.sin(2*numpy.pi*(float(aux[0]) + time)) - float(aux[1])
        sum_error += error*error
        FILE_point2.write(aux[0] + " " + str(E*error) + "\n")
    else:
        FILE_point2.write("\n")
        sum_error = numpy.sqrt(sum_error*step_x2)
        Norm_error2.append((time, sum_error))

# Determine which list is longer
if len(Norm_error1) > len(Norm_error2):
    longer_list = Norm_error1
    shorter_list = Norm_error2
else:
    longer_list = Norm_error2
    shorter_list = Norm_error1

time = []
log2_NormError = []

# Loops through the longer list
for i in range(1, len(longer_list)):
    for j in range(1, len(shorter_list)):
        if longer_list[i][0] == shorter_list[j][0]:
            time.append(longer_list[i][0])
            log2_NormError.append(numpy.log2((shorter_list[j][1])/(longer_list[i][1])))
            break


plt.plot(time, log2_NormError)
plt.ylim([0,4])
plt.xlabel("Time [s]")
if NPoints1 > NPoints2:
    plt.ylabel("Log2(Norm_Error-" + str(NPoints2) + "p/Norm_Error-" + str(NPoints1) + "p)")
else:
    plt.ylabel("Log2(Norm_Error-" + str(NPoints1) + "p/Norm_Error-" + str(NPoints2) + "p)")
plt.savefig("Convergence/Norm_Convergence_Comparison.png")

# Closes the files
FILE_point1.close()
FILE_point2.close()
FILE_norm.close()