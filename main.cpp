#include "Wave_Equation/Wave_Equation.h"

#include <fstream>
#include <math.h>
#include <string>

int main(){
    double* u0 = new double[100];
    double c = 1.0;

    for(int i = 0; i < 49; ++i){
        u0[i] = sin(2.0*M_PI*(double)i/49.0);
        u0[50+i] = cos(2.0*M_PI*(double)i/49.0);
    }

    u0[49] = u0[0];
    u0[99] = u0[50];

    Runge_Kutta_4(&Wave_Equation, u0, 100, 1.0/49.0, (void*) &c, 0.3, 1e-3, "Test_7.dat");

    delete[] u0;

    return 0;
}