// rungekutta5thorder.cpp
// Thomas Watts, Oct. 20th 2020
// CPU: Intel Xeon E5-1607 @ 3 GHz
// Compiler: GCC via Mingw-w64 (VSCode)
//
// This program implements the 5th Runge-Kutta Method to solve
// dy/dt = f(y,t) for an arbitrary vector-valued function f(y,t)

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

using namespace std;

//Constants
const double me = 5.976e24, mm = me / 81, ms = 2.0e4; //Masses of the Earth, Moon, and Spaceship
const double G = 6.674e-11; //The Gravitation Constant

const int num_body = 3; //Number of bodies for the n-body equations

//This function computes f(y,t) = (y[1], (-w^2)y[0]) = k 
void frhs_harmonic(double y[], double t, double k[]) {
    double w = 1.0; //Set angular freqency to 1.0
    k[0] = y[1]; // k = f(y,t) 
    k[1] = -w * w * y[0];
    return;
}


//3 Body Problem Right-hand side function
void frhs_3body(double y[], double t, double k[]) {
    double temp_sum_x, temp_sum_y, dx, dy, dij;
    int x0 = 0, y0 = num_body, vx0 = 2 * num_body, vy0 = 3 * num_body;
    double mass[3];
    mass[0] = me, mass[1] = mm, mass[2] = ms;

    for (int i = 0; i < num_body; i++) {
        k[i + x0] = y[i + vx0];
        k[i + y0] = y[i + vy0];

        temp_sum_x = 0, temp_sum_y = 0;
        for (int j = 0; j < num_body; j++) {
            dx = y[j + x0] - y[i + x0];
            dy = y[j + y0] - y[i + y0];
            dij = sqrt(dx * dx + dy * dy);
            if (j != i) {
                temp_sum_x += G * mass[j] * dx / dij / dij / dij;
                temp_sum_y += G * mass[j] * dy / dij / dij / dij;
            }
        }

        k[i + vx0] = temp_sum_x;
        k[i + vy0] = temp_sum_y;
    }
}


//This function implements the 5th order Runge-Kutta Method
void ode5(double yold[], double ynew[], double error[], double dynew[], double scale[], double h, double t, int n, void frhs(double[], double, double[])) {
    int i; //Looping index
    double* k1, * k2, * k3, * k4, * k5, * k6, * temp; //Initialize pointers
    k1 = new double[7 * n];
    k2 = k1 + n;
    k3 = k2 + n;
    k4 = k3 + n;
    k5 = k4 + n;
    k6 = k5 + n;
    temp = k6 + n;

    //Computer k1,k2,k3,k4,k5,k6
    frhs(yold, t, k1);
    for (i = 0; i < n; i++) temp[i] = yold[i] + 0.2 * h * k1[i];
    frhs(temp, t + 0.2 * h, k2);
    for (i = 0; i < n; i++) temp[i] = yold[i] + (1.0 / 40.0) * h * (3.0 * k1[i] + 9.0 * k2[i]);
    frhs(temp, t + 0.3 * h, k3);
    for (i = 0; i < n; i++) temp[i] = yold[i] + 0.1 * h * (3.0 * k1[i] - 9.0 * k2[i] + 12.0 * k3[i]);
    frhs(temp, t + 0.6 * h, k4);
    for (i = 0; i < n; i++) temp[i] = yold[i] + (1.0 / 729.0) * h * (226.0 * k1[i] - 675.0 * k2[i] + 880.0 * k3[i] + 55.0 * k4[i]);
    frhs(temp, t + (2.0 / 3.0) * h, k5);
    for (i = 0; i < n; i++) temp[i] = yold[i] + (1.0 / 2970.0) * h * (-1991.0 * k1[i] + 7425.0 * k2[i] - 2660.0 * k3[i] - 10010.0 * k4[i] + 10206.0 * k5[i]);
    frhs(temp, t + h, k6);


    //Calculate ynew, dynew, error, scale using k1,k2,k3,k4,k5,k6
    for (i = 0; i < n; i++) {
        dynew[i] = (31.0 / 540.0 * k1[i] + 190.0 / 297.0 * k3[i] - 145.0 / 108.0 * k4[i] + 351.0 / 220.0 * k5[i] + 1.0 / 20.0 * k6[i]);
        ynew[i] = yold[i] + h * dynew[i];
        error[i] = (77.0 * k1[i] - 400.0 * k3[i] + 1925.0 * k4[i] - 1701.0 * k5[i] + 99.0 * k6[i]) / 2520.0;
        scale[i] = abs(ynew[i]) + abs(k1[i]) + abs(1e-30);
    }



    delete[] k1;
}

//Modified rk5 from Numerical Recipes
void Stepperode5(double y[], double dynew[], int n, double* t, double tol, double h_init, double* h_act, double scale[], void frhs(double[], double, double[])) {
    void ode5(double yold[], double ynew[], double error[], double dynew[], double scale[], double h, double t, int n, void frhs(double[], double, double[]));
    //Initialize variables
    int i; double* error, * ytemp;
    double h = h_init, e_max = 0.0;
    error = new double[n];
    ytemp = new double[n];

    //Infinite loop as long as sufficient time step h is found
    for (;;) {
        ode5(y, ytemp, error, dynew, scale, h, *t, n, frhs);

        //Max error
        e_max = 0.0;
        for (i = 0; i < n; i++) {
            if (abs(error[i] / scale[i]) > e_max) { e_max = abs(error[i] / scale[i]); }
        }

        //Set h to appropiate value given the following conditionals
        //exit loop if max error < tolerance
        if (e_max < .2 * tol) h *= pow((tol / e_max), .2);
        if (e_max > tol) h *= .5;
        if (e_max < tol) break;
    }
    //Update time step using new compute time step h
    *t += (*h_act = h);
    for (i = 0; i < n; i++) y[i] = ytemp[i];
    delete[] error, ytemp;

}




int main() {
    void Stepperode5(double y[], double dynew[], int n, double* t, double tol, double h_init, double* h_act, double scale[], void frhs(double[], double, double[]));
    void frhs_3body(double y[], double t, double k[]);
    void frhs_harmonic(double y[], double t, double k[]);

    
    ofstream body3;
    /*
    cos.open("cos.dat");  // open new file for output
    if (cos.fail()) { // or fp.bad()
        cout << "cannot open file" << endl;
        return(EXIT_SUCCESS);
    }
    */
    body3.open("body3.dat");  // open new file for output
    if (body3.fail()) { // or fp.bad()
        cout << "cannot open file" << endl;
        return(EXIT_SUCCESS);
    }
    

    //Initialize Variables
    const int n = 2 * 2 * num_body;
    double h_init = .01, h_act = 0, tol = 1e-6;
    int t_final = 5.184e6;
    

    //Initialize Variables for Harmonic Oscillator Test
    //double y_harmonic[2], dynew_harmonic[2], scale_harmonic[2];
    //Initial Conditions
    //y_harmonic[0] = 1;
    //y_harmonic[1] = 0;


    //Initialize Variables for 3 Body Problem
    double y[12], dynew[12], scale[12];
    //Initial Conditions
    y[0] = 0.0; //Earth x
    y[1] = 0.0; //Moon x
    y[2] = 1.1e7; //Spacecraft x
    y[3] = 0.0; //Earth y
    y[4] = 3.84e8; //Moon y
    y[5] = 1.0e7; //Spacecraft y
    y[6] = -12.593; //Earth vx
    y[7] = 1020.0; //Moon vx
    y[8] = 7.17e3; //Spacecraft vx
    y[9] = 0.0; //Earth vy
    y[10] = 0.0; //Moon vy
    y[11] = 6.08e2; //Spacecraft vy

    //Compute n times steps y(t_n) using 5-th order Runge-Kutta with automatic step size
    for (double t = 0.0; t < t_final; t += h_act) {
        //Stepperode5(y_harmonic, dynew_harmonic, 2, &t, tol, h_init, &h_act, scale_harmonic, frhs_harmonic);

        //cout << std::setprecision(9) << setw(20) << t << setw(20) << *y_harmonic << setw(20) << *dynew_harmonic << endl;
        //cos << std::setprecision(9) << setw(15) << t << setw(15) << *y_harmonic << setw(15) << *dynew_harmonic << endl;

        Stepperode5(y, dynew, n, &t, tol, h_init, &h_act, scale, frhs_3body);

        body3 << std::setprecision(9) << setw(15) << y[0] << setw(15) << y[1] << setw(15) << y[2] << setw(15) << y[3] << setw(15) << y[4] << setw(15) << y[5]  << endl;

        h_init = h_act;
    }

    //Close output files
    //cos.close();
    body3.close();


    return EXIT_SUCCESS;
} // end main()