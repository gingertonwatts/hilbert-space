// MonteCarloMethods.cpp
// Thomas Watts, Nov. 29th 2020
// Homework 9 Finding minimum energy configuration of atoms using Monte Carlo methods
// CPU: Intel Xeon E5-1607 @ 3 GHz
// Compiler: GCC via Mingw-w64 (VSCode)
//
// This program implements a Monte Carlo method for finding the minimum energy of a given potential function U(r)

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
using namespace std;

/*
Define constants
*/
const double dmax=0.25, Lx = 18.0, Ly = 18.0, kbT=0.15;
const int Na=250;

/*
Random Number Generator
*/
typedef double Doub; 
typedef int Int; 
typedef unsigned int Uint; 
typedef unsigned long long Ullong;

struct Ranq1 {Ullong v; Ranq1(Ullong j) : v(4101842887655102017LL) {
		v ^= j;
		v = int64();
	}
	inline Ullong int64() {
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}
	inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};

// Generate Initial Random Atomic Configuration and save to two 2D arrays a and anew
void initConfig(double a[2][Na], double anew[2][Na], unsigned int iseed1, unsigned int iseed2) {
	Ranq1 rx(iseed1), ry(iseed2);
	for (int i = 0; i < Na; i++) {
		a[0][i] = rx.doub() * Lx;
		a[1][i] = ry.doub() * Ly;
		anew[0][i] = a[0][i];
		anew[1][i] = a[1][i];
	}
}

//Calculates 6-12 Potential U(r)
double U(double r) {
	double inv = 1.0 / r;
	double result = inv * inv * inv * inv * inv * inv;
	return 4.0 * result * (result - 1);
}


//Calculates Total Energy of Atomic Configuration a
double Etot(double a[2][Na]) {
	double result = 0.0;
	for (int i = 0; i < Na; i++) {
		for (int j = i+1; j < Na; j++) {
			result += U(sqrt((a[0][i] - a[0][j]) * (a[0][i] - a[0][j]) + (a[1][i] - a[1][j]) * (a[1][i] - a[1][j])));
		}
	}
	return result;
}


//Generates new atomic configuration at random
void generateConfig(double a[2][Na],double r1, double r2, double r3) {
	double direction = r1;
	int atom = (int)(250.0 * r2);
	double h = dmax*(2.0*r3-1.0);
	if (direction >= 0.5) a[1][atom] += h;
	else a[0][atom] += h;
}

//Copy configuruation a to array anew
void changeConfig(double a[2][Na], double anew[2][Na]) {
	for (int i = 0; i < Na; i++) {
		a[0][i] = anew[0][i];
		a[1][i] = anew[1][i];
	}
}


int main() {
	//Create output files
	ofstream fp1, fp2, fp3, fp4, fp5;
	fp1.open("histo1.dat");  
	if (fp1.fail()) {	
		cout << "cannot open file" << endl;
		return(EXIT_SUCCESS);
	}
	fp2.open("histo2.dat");  
	if (fp2.fail()) {	
		cout << "cannot open file" << endl;
		return(EXIT_SUCCESS);
	}
	fp3.open("E.dat");  
	if (fp3.fail()) {	
		cout << "cannot open file" << endl;
		return(EXIT_SUCCESS);
	}
	fp4.open("x1.dat");  
	if (fp4.fail()) {	
		cout << "cannot open file" << endl;
		return(EXIT_SUCCESS);
	}
	fp5.open("y1.dat");  
	if (fp5.fail()) {	
		cout << "cannot open file" << endl;
		return(EXIT_SUCCESS);
	}


	//Collect RNG histogram data
	Ranq1 rand(17);
	const int num_bin = 100;
	double bin[num_bin+1];
	double num = 0, bin_size = 1.0/num_bin;
	for (int i = 0; i < 100; i++) bin[i] = 0;

	for (int i = 0; i < 100000; i++) {
		num = rand.doub();
		for (int j = 0; j < num_bin+1; j++) {
 			if (num < j * bin_size && num > (j - 1) * bin_size) {
				bin[j - 1] += 1;
				break;
			}
		}
	}

	//Save histogram data
	for (int i = 0; i < num_bin; i++) {
		fp1 << std::setprecision(9) << setw(15) << i * bin_size << setw(15) << (i + 1) * bin_size << setw(15) << bin[i] << endl;
	} fp1.close();

	//Save (x[i], x[i+1]) data
	double arr[10001];
	for (int i = 0; i < 10001; i++) arr[i] = rand.doub();	
	for (int i = 0; i < 10000; i++) fp2<< std::setprecision(9) << setw(15) << arr[i] << setw(15) << arr[i + 1] << endl;
	fp2.close();

	//Choose seeds for RNG
	unsigned int iseed1 = 13, iseed2 = 17, iseed3= 19, iseed4 = 15;

	//Create arrays to store atomic configurations
	double a[2][Na], anew[2][Na];

	//Choose random initial atomic configuration
	initConfig(a, anew, iseed1, iseed4);

	//Calculate initial energy E0
	double E0 = Etot(a), Enew;

	//Perform (num_itr) Monte Carlo iterations
	int num_itr = 1e+6;
	Ranq1 r1(iseed1), r2(iseed2), r3(iseed3), r4(iseed4);
	for (int i = 1; i <= num_itr; i++) {
		//Generate new configuration and calculate its energy Enew
		generateConfig(anew, r1.doub(),r2.doub(),r3.doub());
		Enew=Etot(anew);
		double dE = Enew-E0;
		//Decide whether to accept or reject new configuration
		if (dE <= 0.0) {
			E0 = Enew;
			changeConfig(a,anew);
		}
		else {
			if (exp(-dE/kbT) > r4.doub()) {
				E0 = Enew;
				changeConfig(a,anew);
			}
			else {
				//Reset atomic configuration anew to a if the new configuration is rejected
				changeConfig(anew,a);
			}
		}
		if(i % 2000 == 0){
			fp3 << std::setprecision(9) << setw(15) << E0 << endl;
		}
		if(i==1e+3 || i==1e+4 || i==1e+5 || i==1e+6){
			for(int j=0; j<Na; j++){
				fp4 << std::setprecision(9) << setw(15) << a[0][j] << endl;
				fp5 << std::setprecision(9) << setw(15) << a[1][j] << endl;
			}
		}

	} 

	cout << E0;

	fp3.close();
	fp4.close();
	fp5.close();

	return(EXIT_SUCCESS);
}    