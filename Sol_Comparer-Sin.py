import numpy
import sys

# Initializes the filename and the number of points to a default
filename = "Pointwise_Error.dat"

# Reads the command line arguments to define the file to compare to the analytical solution
input = sys.argv[1]

# Multiplicative factor declaration (for later error comparison)
E = 1

# Reads the command line arguments to define variables needed for the comparison
for i in range(2, len(sys.argv)):
    if sys.argv[i] == "-EQ":
        N_Eqs = int(sys.argv[i+1])

    if sys.argv[i] == "-NP":
        NPoints = int(sys.argv[i+1])

    if sys.argv[i] == "-FN":
        filename = sys.argv[i+1]

    if sys.argv[i] == "-E":
        E = int(sys.argv[i+1])

# Opens the file that will be written
FILE = open(filename, "w")

# Opens the file that will be compared
IN = open(input, "r")

# Writes in the output file the error of the parameter file compared to the analitical solution 
time = 0.0

for l in IN:
    if l[0] == "\"":
        aux = l.split()
        time = float(aux[2])
        FILE.write(l)

    elif l != "\n":
        aux = l.split()
        error = numpy.sin(2*numpy.pi*float(aux[0]) + time) - float(aux[1])
        FILE.write(aux[0] + " " + str(E*error) + "\n")
    else:
        FILE.write("\n")

# Closes the file
FILE.close()