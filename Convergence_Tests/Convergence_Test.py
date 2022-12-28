import os

os.system("echo Beginning convergence analysis")

os.system("mkdir Convergence")

os.system("echo Writing the inital conditions for 50 points")
os.system("python3 Convergence_Tests/IC_Writer-Sin.py -NP 50 -FN Convergence/IC-Sin-50p.dat")

os.system("echo Solving the equation using this initial conditions")
os.system("./test.exe -IC Convergence/IC-Sin-50p.dat -EQ simple_wave -FN Convergence/Out-Sin-50p.dat")

os.system("echo \n")

os.system("echo Writing the inital conditions for 100 points")
os.system("python3 Convergence_Tests/IC_Writer-Sin.py -NP 100 -FN Convergence/IC-Sin-100p.dat")

os.system("echo Solving the equation using this initial conditions")
os.system("./test.exe -IC Convergence/IC-Sin-100p.dat -EQ simple_wave -FN Convergence/Out-Sin-100p.dat")

os.system("echo \n")

os.system("echo Calculating the convergence")
os.system("python3 Convergence_Tests/Sol_Comparer-Sin.py Convergence/Out-Sin-50p.dat Convergence/Out-Sin-100p.dat -SX 0.02 0.01 -NP 50 100 -FN Convergence/Sin- -E 4")

os.system("Creating tar file with the results")


os.system("echo Displaying pointwise convergence for both solutions")
os.system("echo Note: The error in the solution with 100 points is multiplied by 4")
os.system("muninn Convergence/Sin-50p-Point_Error.dat Convergence/Sin-100p-Point_Error.dat")