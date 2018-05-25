#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include "Interpolation.h"

using namespace std;

int main(int argc, char *argv[]) {
	//if (argc < 4) return 0;

	//int delta = atoi(argv[1]);
	//string inputFile = argv[2];
	//string outputFile = argv[3];
	string inputFile = "./przyk1.txt";
	string outputFile = "./przyk1";

	Interpolation inter;
	cout << "Liczba wybranych punktow: ";
	cout << inter.lagrange(inputFile.c_str(), (outputFile + ".lagr").c_str(), 0, 15);
	inter.spline3Deg(inputFile.c_str(), (outputFile + ".spl3").c_str(), 0, 15);
	
	getchar();
	return 0;
}