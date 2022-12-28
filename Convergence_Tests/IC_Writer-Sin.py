import numpy
import sys

# Initializes the filename and the number of points to a default
filename = "IC.dat"
NPoints = 50

# Reads the command line arguments to define the previous variables
for i in range(1, len(sys.argv)):
    if sys.argv[i] == "-NP":
        NPoints = int(sys.argv[i+1])

    if sys.argv[i] == "-FN":
        filename = sys.argv[i+1]

# Opens the file that will be written
FILE = open(filename, "w")

# Writtes the number of equations in the ODE systemspace step, the number of points to the file and the parameters for the equation
FILE.write("#NEq: 2\n")
FILE.write("#step_x: ")
FILE.write(str(1/NPoints))
FILE.write("\n")
FILE.write("#NPoints: ")
FILE.write(str(NPoints+1))
FILE.write("\n")
FILE.write("#pars: 1 1.0\n")

# Writtes the initial conditions for the equation
for i in range( NPoints):
    FILE.write(str(numpy.sin(2*numpy.pi*i/NPoints)))
    FILE.write(" ")
    FILE.write(str(2*numpy.pi*numpy.cos(2*numpy.pi*i/NPoints)))
    FILE.write("\n")

FILE.write("0.0 ")
FILE.write(str(2*numpy.pi))

# Closes the file
FILE.close()