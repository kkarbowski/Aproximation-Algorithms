#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include "Interpolation.h"

using namespace std;

int main(int argc, char *argv[]) {
	if (argc < 5) return 0;

	string inputFile = argv[1];
	string outputFile = argv[2];
	int delta = atoi(argv[3]);
	int delta2 = atoi(argv[4]);

	Interpolation inter;
	cout << "Liczba wybranych punktow: ";
	cout << inter.lagrange(inputFile.c_str(), (outputFile + ".lagr").c_str(), delta, delta2);
	inter.spline3Deg(inputFile.c_str(), (outputFile + ".spl3").c_str(), delta, delta2);
	
	getchar();
	return 0;
}