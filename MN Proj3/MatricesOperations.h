#pragma once
#include "Matrix.h"

class MatricesOperations
{
private:
	int iterations;
	double duration;

public:
	MatricesOperations();

	Matrix createBandMatrix(double a1, double a2, double a3, int n);
	Matrix createIdentityMatrix(int n);
	Matrix createSpecVectorB(int f, int n);

	Matrix solveJacobi(Matrix a, Matrix vectorB, double residuumNormLimit);
	Matrix solveGS(Matrix a, Matrix vectorB, double residuumNormLimit);
	Matrix solveLUfactorization(Matrix a, Matrix vectorB);

	double calculateVectorNorm(Matrix vector);

	double getDurationTime();
	int getIterations();

	~MatricesOperations();
};

