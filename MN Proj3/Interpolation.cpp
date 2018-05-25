#include "Interpolation.h"
#include "MatricesOperations.h"
#include "Matrix.h"
#include <cstdio>
#include <vector>
#include <cmath>

Interpolation::Interpolation()
{
}

int Interpolation::lagrange(const char * inputFile, const char * outputFile, int delta)
{
	readData(inputFile);

	if (points.size() == 0)
		return 0;

	for (int i = 0; i < points.size(); ++i)
		result.push_back(Point{ points[i].x, interpolationLagrangeFunc(points[i].x, delta) });

	writeData(outputFile);
	clearData();
	return 1;
}

int Interpolation::spline3Deg(const char * inputFile, const char * outputFile, int delta)
{
	readData(inputFile);

	if (points.size() == 0)
		return 0;

	generateSplineFuncs(delta);

	for (int i = 0; i < points.size(); ++i)
		result.push_back(Point{ points[i].x, splineCalcValue(points[i].x, delta) });

	writeData(outputFile);
	pickedPoints.clear();
	return 1;
}

double Interpolation::fiFunc(double x, int index, int delta)
{
	double product = 1.0;
	for (int i = 0; i < points.size(); i += delta)
		if (i != index)
			product *= (x - points[i].x) / (points[index].x - points[i].x);
	return product;
}

double Interpolation::interpolationLagrangeFunc(double x, int delta)
{
	double sum = 0.0;
	for (int i = 0; i < points.size(); i += delta)
		sum += points[i].y*fiFunc(x, i, delta);
	return sum;
}

void Interpolation::generateSplineFuncs(int delta)
{
	//vector<Point> pickedPoints;
	int subintervals = points.size() / delta;
	Matrix aMat(4 * subintervals, 4 * subintervals, 0);
	Matrix bVect(4 * subintervals, 1, 0);

	for (int i = 0; i < subintervals + 1; ++i)
		pickedPoints.push_back(points[i*delta]);

	// 1 eq
	for (int i = 0; i < subintervals; ++i) {
		aMat(i * 4, i * 4) = 1;
		bVect(i * 4) = pickedPoints[i].y;
	}
	// 2 eq
	for (int i = 0; i < subintervals; ++i) {
		//int col = i - subintervals;
		double h = pickedPoints[i + 1].x - pickedPoints[i].x;

		aMat(i * 4 + 1, i * 4) = 1;
		aMat(i * 4 + 1, i * 4 + 1) = h;
		aMat(i * 4 + 1, i * 4 + 2) = h*h;
		aMat(i * 4 + 1, i * 4 + 3) = h*h*h;
		bVect(i * 4 + 1) = pickedPoints[i + 1].y;
	}
	// 3 eq
	for (int i = 0; i < subintervals - 1; ++i) {
		//int col = i - 2 * subintervals;
		double h = pickedPoints[i + 1].x - pickedPoints[i].x;

		aMat(i * 4 + 3, i * 4 + 1) = 1;
		aMat(i * 4 + 3, i * 4 + 2) = 2 * h;
		aMat(i * 4 + 3, i * 4 + 3) = 3 * h*h;
		aMat(i * 4 + 3, i * 4 + 5) = -1;
	}
	// 4 eq
	for (int i = 0; i < subintervals - 1; ++i) {
		//int col = i - 3 * subintervals + 1;
		double h = pickedPoints[i + 1].x - pickedPoints[i].x;

		aMat(i * 4 + 6, i * 4 + 2) = 2;
		aMat(i * 4 + 6, i * 4 + 3) = 6 * h;
		aMat(i * 4 + 6, i * 4 + 6) = -2;
	}
	// 5 eq
	{
		aMat(2, 2) = 1;
	}
	// 6 eq
	{
		int i = 4 * subintervals - 1;
		double h = pickedPoints[subintervals].x - pickedPoints[subintervals - 1].x;
		aMat(i, i - 1) = 2;
		aMat(i, i) = 6 * h;
	}
	resultMat = matOp.solveLUfactorization(aMat, bVect);
	//return matOp.solveLUfactorization(aMat, bVect);
}

double Interpolation::splineCalcValue(double x, int delta)
{
	int subintervals = points.size() / delta - 1;
	//int resSubInt = 0;

	int i = 0;
	for (; i < subintervals; ++i)
		if (x < points[(i + 1)*delta].x)
			break;

	double a = resultMat(4 * i),
		b = resultMat(4 * i + 1),
		c = resultMat(4 * i + 2),
		d = resultMat(4 * i + 3),
		xi = pickedPoints[i].x;

	return a + b*(x - xi) + c*pow((x - xi), 2) + d*pow((x - xi), 3);
}

Interpolation::~Interpolation()
{
}

void Interpolation::readData(const char * path)
{
	if (FILE *pFile = fopen(path, "r")) {
		double x, y;
		while (fscanf(pFile, "%lf\t%lf", &x, &y) > 0)
			points.push_back(Point{ x, y });
	}
}

void Interpolation::writeData(const char * path)
{
	if (FILE *pFile = fopen(path, "w"))
		for (int i = 0; i < result.size(); ++i)
			fprintf(pFile, "%f\t%f\n", result[i].x, result[i].y);
}

void Interpolation::clearData()
{
	points.clear();
	result.clear();
}
