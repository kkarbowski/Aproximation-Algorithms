#include <cstdio>
#include <cstdlib>
#include <string>
#include "Interpolation.h"

using namespace std;

int main(int argc, char *argv[]) {
	//if (argc < 4) return 0;

	//int delta = atoi(argv[1]);
	//string inputFile = argv[2];
	//string outputFile = argv[3];
	string inputFile = "./data1.txt";
	string outputFile = "./data1";

	Interpolation inter;
	//inter.lagrange(inputFile.c_str(), (outputFile + ".lagr").c_str(), 50);
	inter.spline3Deg(inputFile.c_str(), (outputFile + ".spl3").c_str(), 5);


	return 0;
}