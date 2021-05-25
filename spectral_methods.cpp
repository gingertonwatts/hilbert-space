// spectral_methods.cpp
// Thomas Watts, Nov. 29th 2020
// CPU: Intel Xeon E5-1607 @ 3 GHz
// Compiler: GCC via Mingw-w64 (VSCode)
//
// Run file using terminal command g++ -O -o hw10 hw10f.cpp -lfftw3
// Get output with terminal command ./hw10
//
// This program uses the FFT to solve the 2D wave equation

#include<iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include <complex.h>
using namespace std;

#include<fftw3.h>

/*
Define constants
*/
const double Lx=500.0, Ly=500.0, v=343.0;
const int Nx=512, Ny=512;

/*
Computes square Euclidean distance between the 2D vectors x and xn
*/
double sqrdist(double x[], double xn[]){
    return (x[0]-xn[0])*(x[0]-xn[0])+(x[1]-xn[1])*(x[1]-xn[1]);
}

/*
Computes cos(v*|kij|*t) for a given i,j and time t
*/
double factor(int i, int j, double t){
    double term=sqrt(((double)(i)/(double)(Lx))*((double)(i)/(double)(Lx))+((double)(j)/(double)(Ly))*((double)(j)/(double)(Ly)));
    return cos(v*term*t);
}

/*
Initial Condition psi(x,0) := three small Gaussians
*/
double psi0(double x[]){
    double xn1[2] = {0.4*Lx,0.4*Ly}, xn2[2] = {0.5*Lx,0.6*Ly}, xn3[2] = {0.6*Lx, 0.4*Ly}; 
    return exp(-sqrdist(x,xn1)/100.0)+2.0*exp(-sqrdist(x,xn2)/400.0)+exp(-sqrdist(x,xn3)/100.0);
}

int main(){
    // complex 2D waves
    fftw_complex *y2, *y2new;
    // FFTW plans
    fftw_plan planTf2, planTi2;

 
    y2 = (fftw_complex*) fftw_malloc(Nx*Ny*sizeof(fftw_complex)); // in place of C++ new
    if( NULL == y2 ) {
    //<check for NULL == y2 etc.>
        cout << "Cannot allocate array y2" << endl;
        exit( EXIT_FAILURE );
    }
    y2new = (fftw_complex*) fftw_malloc(Nx*Ny*sizeof(fftw_complex)); // in place of C++ new
    if( NULL == y2new ) {
    //<check for NULL == y2 etc.>
        cout << "Cannot allocate array y2" << endl;
        exit( EXIT_FAILURE );
    }

    // find a plan once for each size (for in-place FFT)
    // slower execution but fast plan
    planTf2 = fftw_plan_dft_2d(Nx, Ny, y2, y2, FFTW_FORWARD, FFTW_ESTIMATE );
    planTi2 = fftw_plan_dft_2d(Nx, Ny, y2new, y2new, FFTW_BACKWARD, FFTW_ESTIMATE );

    // initialize y2[j + i*Ny][0] = real part and y2[j + i*Ny][1] = imag. part etc.
    for(int i = 0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            double x[2]={Lx*(double)(i)/(double)(Nx),Ly*(double)(j)/(double)(Ny)};
            y2[j+i*Ny][0]=psi0(x);
            y2[j+i*Ny][1]=0.0;
        }
    }

    //Open output files
    ofstream fp1, fp2, fp3, fp4;
	fp1.open("t0.dat");  
	if (fp1.fail()) {	
		cout << "cannot open file" << endl;
		return(EXIT_SUCCESS);
	}
    fp2.open("t1.dat");  
	if (fp2.fail()) {	
		cout << "cannot open file" << endl;
		return(EXIT_SUCCESS);
	}
    fp3.open("t2.dat");  
	if (fp3.fail()) {	
		cout << "cannot open file" << endl;
		return(EXIT_SUCCESS);
	}
    fp4.open("t3.dat");  
	if (fp4.fail()) {	
		cout << "cannot open file" << endl;
		return(EXIT_SUCCESS);
	}

    for(int i = 0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            fp1 << std::setprecision(9) << setw(15) << y2[j+i*Ny][0] << endl;
        }
    }


    //Compute Forward FFT of psi0
    fftw_execute_dft( planTf2, y2, y2 );

    
    //t=0.1
    double t1=0.1;
    for(int i = 0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            double term = factor(i,j,t1);
            y2new[j+i*Ny][0]=y2[j+i*Ny][0]*term;
            y2new[j+i*Ny][1]=y2[j+i*Ny][1]*term;
        }
    }

    //Generate unnormalized solution at t=0.1 with Backwards FFT
    fftw_execute_dft(planTi2, y2new, y2new);

    //Save unnormalized solution at t=0.1
    for(int i = 0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            fp2 << std::setprecision(9) << setw(15) << y2new[j+i*Ny][0] << endl;
        }
    }

    //t=0.2
    double t2=0.2;
    for(int i = 0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            double term = factor(i,j,t2);
            y2new[j+i*Ny][0]=y2[j+i*Ny][0]*term;
            y2new[j+i*Ny][1]=y2[j+i*Ny][1]*term;
        }
    }

    //Generate unnormalized solution at t=0.2 with Backwards FFT
    fftw_execute_dft(planTi2, y2new, y2new);

    //Save unnormalized solution at t=0.2
    for(int i = 0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            fp3 << std::setprecision(9) << setw(15) << y2new[j+i*Ny][0] << endl;
        }
    }


    //t=0.4
    double t3=0.4;
    for(int i = 0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            double term = factor(i,j,t3);
            y2new[j+i*Ny][0]=y2[j+i*Ny][0]*term;
            y2new[j+i*Ny][1]=y2[j+i*Ny][1]*term;
        }
    }

    //Generate unnormalized solution at t=0.4 with Backwards FFT
    fftw_execute_dft(planTi2, y2new, y2new);


    //Save unnormalized solution at t=0.4
    for(int i = 0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            fp4 << std::setprecision(9) << setw(15) << y2new[j+i*Ny][0] << endl;
        }
    }

    //Close output files
    fp1.close();
    fp2.close();
    fp3.close();
    fp4.close();

    return (EXIT_SUCCESS);
}