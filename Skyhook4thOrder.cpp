// Skyhook4thOrder.cpp
// Thomas Watts, Dec. 18th 2020
// CPU: Intel Xeon E5-1607 @ 3 GHz
// Compiler: GCC via Mingw-w64 (VSCode)
//
// This program implements the 4th Runge-Kutta Method to solve
// the equations of motion of a simple model for a Skyhook (Coleman, Colin.  (2019).  Space Access for Future Planetary Science Missions.  10.5772/intechopen.88530)

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

using namespace std;

//Constants
const double mE = 5.976e24; //Masses of the Earth and Spaceships (kg)
const int n = 6;
const double G = 6.674e-11; //The Gravitation Constant
const double m=1000, L=1090000.0, mass=5.46, moment=5.06;

//Compute mT and I from mT/m and IT/mL*L
const double mT = m*mass, I = m*L*L*(moment+2);
const double c = G*mE*m/(2*m+mT);

//Returns polar coordinates of the tether masses as an double array {r1,theta1,r2,theta2}
void polarCoords(double r, double theta, double phi, double L, double result[]){
    double r1x, r1y, r2x, r2y;
    r1x=r*cos(theta)-L*cos(phi);
    r1y=r*sin(theta)-L*sin(phi);
    r2x=r*cos(theta)+L*cos(phi);
    r2y=r*sin(theta)+L*sin(phi);

    //Compute polar coordinates of the tether masses
    result[0]=sqrt(r1x*r1x+r1y*r1y);
    result[1]=atan2(r1y,r1x);
    result[2]=sqrt(r2x*r2x+r2y*r2y);
    result[3]=atan2(r2y,r2x);
}

//This function computes f(y,t) = k 
//where k = (r,theta,phi,r',theta',phi')
//mass = mT/m from materials table, moment = IT/m*L*L from materials table
void frhs(double y[], double t, double k[]) {
    //create list for polar coordinates of tether masses
    double polar[4];

    //Compute r1, theta1, r2, theta2
    polarCoords(y[0],y[1],y[2],L,polar);

    //set (k[0],k[1],k[2]) = (r',theta',phi')
    k[0]=y[3];
    k[1]=y[4];
    k[2]=y[5];

    //The expression for f(y,t) is quite long so I break it up into multiple terms here

    //term appearing in k[3]
    double term1 = (1/polar[0])+(1/polar[2]); //1/r1 + 1/r2

    //Compute r1^3
    double r1cubed = polar[0]*polar[0]*polar[0];
    //Compute r2^3
    double r2cubed = polar[2]*polar[2]*polar[2];
    //Compute r1^2 * r2^2
    double r1r2sqrd = polar[0]*polar[0]*polar[2]*polar[2];

    double term2 = cos(polar[1]-polar[3])*(r1cubed+r2cubed)/(r1r2sqrd); //cos(theta1-theta2)*(r1^3+r2^3/r1^2*r2^2)
    double term3 = (2.0*mT)/(y[0]*m); //2*mT/r*m

    k[3]=y[0]*y[4]*y[4]-c*(1.0/(2.0*y[0]))*(term1+term2+term3); //computer r'

    //term appearing in k[4] and k[5]
    double term4 = sin(polar[1]-polar[3])*(r1cubed-r2cubed)/(r1r2sqrd); //sin(theta1-theta2)*(r1^3+r2^3/r1^2*r2^2)
    k[4]=(-2.0*y[3]*y[4]/y[0])+c*(1.0/(2.0*y[0]*y[0]))*term4; //compute theta'
    k[5]=(-G*mE*m/(2.0*I))*term4; //compute phi'

}

//This function implements the 4th order Runge-Kutta Method
void ode4(double yold[], double ynew[], double h, double t, int n, void frhs(double[], double, double[])) {
    int i; //Looping index
    double* k1, * k2, * k3, * k4, * temp; //Initialize pointers
    k1 = new double[5 * n];
    k2 = k1 + n;
    k3 = k2 + n;
    k4 = k3 + n;
    temp = k4 + n;

    //Compute k1,k2,k3,k4
    frhs(yold, t, k1);
    for (i = 0; i < n; i++) temp[i] = yold[i] + 0.5 * h * k1[i]; 
    frhs(temp, t + 0.5 * h, k2);
    for (i = 0; i < n; i++) temp[i] = yold[i] + 0.5 * h * k2[i]; 
    frhs(temp, t + 0.5 * h, k3);
    for (i = 0; i < n; i++) temp[i] = yold[i] + h * k3[i]; 
    frhs(temp, t + h, k4);

    //Compute ynew = yold + h/6(k1+2k2+2k3+k4)
    for (i = 0; i < n; i++)  ynew[i] = yold[i] + h * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) / 6.0;

    //Delete k1 when ode4 is finished running
    delete[] k1;
}

int main() {
    void ode4(double yold[], double ynew[], double h, double t, int n, void frhs(double[], double, double[]));
    void frhs(double y[], double t, double k[]);

    
    ofstream fp1;
    fp1.open("traj.dat");  // open new file for output
    if (fp1.fail()) { // or fp.bad()
        cout << "cannot open file" << endl;
        return(EXIT_SUCCESS);
    }
    

    //Initialize Variables
    double y[6], ynew[6];

    //Initial Conditions
    y[0] = 8000000.0+6378100.0; //r
    y[1] = 0.0; //theta
    y[2] = 0.0; //phi
    y[3] = 0.0; //velocity r
    y[4] = sqrt(G*mE/y[0])/y[0]; //velocity theta (include velocity equation in report (equate Fg to cent. force)) (convert to rad)
    y[5] = 100; //velocity phi

    double h = 0.001; //Timp step
    int timeSteps = 17200000;
  
    //Save t,x(t),y(t),z(t) to text file to be loaded into python
    for (int i = 1; i < timeSteps; i++) {

        //Compute x(t), y(t),z(t)  using ode4 
        //and write them directly to a text file to avoid memory issues
        ode4(y, ynew, h, i * h, 6, frhs);

        //save data every so often
        if(i%1000==0){
            fp1.precision(5);
            for (int j = 0; j < 6; j++){
                fp1 << setw(15) << ynew[j];
            }
            fp1 << endl;
        }

        //Save the computed ynew as the yold of the next time step
        y[0] = ynew[0];
        y[1] = ynew[1];
        y[2] = ynew[2];
        y[3] = ynew[3];
        y[4] = ynew[4];
        y[5] = ynew[5];

    }
    
    //Close output files
    fp1.close();


    return EXIT_SUCCESS;
}
