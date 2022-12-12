#include "Wave_Equation/Wave_Equation.h"

#include <fstream>
#include <math.h>
#include <string>

int main(){
    double* u0 = new double[100];
    double c = 1.0;

    for(int i = 0; i < 49; ++i){
        u0[i] = sin(2.0*M_PI*(double)i/49.0);
        u0[50+i] = 2.0*M_PI*cos(2.0*M_PI*(double)i/49.0);
    }

    u0[49] = u0[0];
    u0[99] = u0[50];

    Runge_Kutta_4(&Wave_Equation, u0, 100, 2, 1.0/49.0, (void*) &c, 1, 1e-3, "Tests/Test_10.dat");

    delete[] u0;

    return 0;
}