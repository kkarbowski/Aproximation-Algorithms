#include <cmath>
#include <chrono>
#include <iostream>
#include "MatricesOperations.h"
#include "Matrix.h"


MatricesOperations::MatricesOperations()
{
}

Matrix MatricesOperations::createBandMatrix(double a1, double a2, double a3, int n)
{
	Matrix mat(n, n, 0);
	mat(0, 0) = a1;
	mat(1, 1) = a1;
	mat(0, 1) = a2;
	mat(1, 0) = a2;

	for (int i = 2; i < n; ++i) {
		mat(i, i) = a1;
		mat(i - 1, i) = a2;
		mat(i, i - 1) = a2;
		mat(i - 2, i) = a3;
		mat(i, i - 2) = a3;
	}
	return mat;
}

Matrix MatricesOperations::createIdentityMatrix(int n)
{
	Matrix mat(n, n, 0);
	for (int i = 0; i < n; ++i)
		mat(i, i) = 1;
	return mat;
}

Matrix MatricesOperations::createSpecVectorB(int f, int n)
{
	Matrix vect(n, 1);
	for (int i = 0; i < n; ++i)
		vect(i) = sin((i + 1)*(f + 1));

	return vect;
}

Matrix MatricesOperations::solveJacobi(Matrix a, Matrix vectorB, double residuumNormLimit)
{
	int height = vectorB.getHeight();
	Matrix resultVect(height, 1);
	Matrix copyResultVect(height, 1, 1);
	Matrix residuum(height, 1);
	double tmp = 0;
	iterations = 0;

	auto begin = std::chrono::high_resolution_clock::now();
	while (true) {
		++iterations;
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < height; ++j)
				tmp += a(i, j)*copyResultVect(j);
			tmp -= a(i, i)*copyResultVect(i);	// (-) j==i element
			resultVect(i) = (vectorB(i) - tmp) / a(i, i);
			tmp = 0;
		}
		residuum = (a * resultVect) - vectorB;
		//{
		//	std::cout << calculateVectorNorm(residuum) << std::endl; //COUT
		//	if (iterations == 9) break;
		//}
		if (calculateVectorNorm(residuum) < residuumNormLimit)
			break;
		copyResultVect = resultVect;
	}
	auto end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()/(double)1000000000;
	return resultVect;
}

Matrix MatricesOperations::solveGS(Matrix a, Matrix vectorB, double residuumNormLimit)
{
	int height = vectorB.getHeight();
	Matrix resultVect(height, 1);
	Matrix copyResultVect(height, 1, 1);
	Matrix residuum(height, 1);
	double tmp = 0;
	iterations = 0;

	auto begin = std::chrono::high_resolution_clock::now();
	while (true) {
		++iterations;
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < i; ++j)
				tmp += (a(i, j)*resultVect(j));
			for (int j = i + 1; j < height; ++j)
				tmp += (a(i, j)*copyResultVect(j));
			resultVect(i) = (vectorB(i) - tmp) / a(i, i);
			tmp = 0;
		}
		residuum = (a * resultVect) - vectorB;
		//{
		//	std::cout << calculateVectorNorm(residuum) << std::endl; //COUT
		//	if (iterations == 9) break;
		//}
		if (calculateVectorNorm(residuum) < residuumNormLimit)
			break;
		copyResultVect = resultVect;
	}
	auto end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()/(double)1000000000;
	return resultVect;
}

Matrix MatricesOperations::solveLUfactorization(Matrix a, Matrix vectorB)
{
	int n = vectorB.getHeight();
	Matrix resultVect(n, 1);
	Matrix vectY(n, 1);
	Matrix matU = a;
	Matrix matL = createIdentityMatrix(n);
	double tmp = 0;
	iterations = -1;

	auto begin = std::chrono::high_resolution_clock::now();
	// Calculating LU
	for (int i = 0; i < n - 1; ++i) {
		for (int j = i + 1; j < n; ++j) {
			matL(j, i) = matU(j, i) / matU(i, i);
			for (int k = i; k < n; ++k)
				matU(j, k) = matU(j, k) - (matL(j, i)*matU(i, k));
		}
	}

	// LY=B
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < i; ++j)
			tmp += matL(i, j)*vectY(j);
		vectY(i) = vectorB(i) - tmp;
		tmp = 0;
	}

	// UX=Y
	for (int i = n - 1; i >= 0; --i) {
		for (int j = i + 1; j < n; ++j)
			tmp += matU(i, j)*resultVect(j);
		resultVect(i) = (vectY(i) - tmp) / matU(i, i);
		tmp = 0;	
	}

	auto end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()/(double)1000000000;
	return resultVect;
}

double MatricesOperations::calculateVectorNorm(Matrix vector)
{
	double tmp = 0;
	for (int i = 0; i < vector.getHeight(); ++i)
		tmp += vector(i) * vector(i);
	return sqrt(tmp);
}

double MatricesOperations::getDurationTime()
{
	return duration;
}

int MatricesOperations::getIterations()
{
	return iterations;
}

MatricesOperations::~MatricesOperations()
{
}
