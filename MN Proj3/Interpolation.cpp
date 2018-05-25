#include "Interpolation.h"
#include "MatricesOperations.h"
#include "Matrix.h"
#include <cstdio>
#include <vector>
#include <cmath>

Interpolation::Interpolation()
{
}

int Interpolation::lagrange(const char * inputFile, const char * outputFile, int delta, int delta2)
{
	readData(inputFile, delta, delta2);

	if (points.size() == 0)
		return 0;

	for (int i = 0; i < points.size(); ++i)
		result.push_back(Point{ points[i].x, interpolationLagrangeFunc(points[i].x) });

	writeData(outputFile);
	int pickedPointsNum = pickedPoints.size();
	clearData();
	return pickedPointsNum;
}

int Interpolation::spline3Deg(const char * inputFile, const char * outputFile, int delta, int delta2)
{
	readData(inputFile, delta, delta2);

	if (points.size() == 0)
		return 0;

	generateSplineFuncs();

	for (int i = 0; i < points.size(); ++i)
		result.push_back(Point{ points[i].x, splineCalcValue(points[i].x) });

	writeData(outputFile);
	int pickedPointsNum = pickedPoints.size();
	clearData();
	return pickedPointsNum;
}

double Interpolation::fiFunc(double x, int index)
{
	double product = 1.0;
	for (int i = 0; i < pickedPoints.size(); ++i)
		if (i != index)
			product *= (x - pickedPoints[i].x) / (pickedPoints[index].x - pickedPoints[i].x);
	return product;
}

double Interpolation::interpolationLagrangeFunc(double x)
{
	double sum = 0.0;
	for (int i = 0; i < pickedPoints.size(); ++i)
		sum += pickedPoints[i].y*fiFunc(x, i);
	return sum;
}

void Interpolation::generateSplineFuncs()
{
	Matrix aMat(4 * subintervals, 4 * subintervals, 0);
	Matrix bVect(4 * subintervals, 1, 0);

	// 1 eq
	for (int i = 0; i < subintervals; ++i) {
		aMat(i * 4, i * 4) = 1;
		bVect(i * 4) = pickedPoints[i].y;
	}
	// 2 eq
	for (int i = 0; i < subintervals; ++i) {
		double h = pickedPoints[i + 1].x - pickedPoints[i].x;

		aMat(i * 4 + 1, i * 4) = 1;
		aMat(i * 4 + 1, i * 4 + 1) = h;
		aMat(i * 4 + 1, i * 4 + 2) = h*h;
		aMat(i * 4 + 1, i * 4 + 3) = h*h*h;
		bVect(i * 4 + 1) = pickedPoints[i + 1].y;
	}
	// 3 eq
	for (int i = 0; i < subintervals - 1; ++i) {
		double h = pickedPoints[i + 1].x - pickedPoints[i].x;

		aMat(i * 4 + 3, i * 4 + 1) = 1;
		aMat(i * 4 + 3, i * 4 + 2) = 2 * h;
		aMat(i * 4 + 3, i * 4 + 3) = 3 * h*h;
		aMat(i * 4 + 3, i * 4 + 5) = -1;
	}
	// 4 eq
	for (int i = 0; i < subintervals - 1; ++i) {
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
}

double Interpolation::splineCalcValue(double x)
{
	int i = 0;
	for (; i < subintervals ? true : --i && false; ++i)
		if (x < pickedPoints[i + 1].x)
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

void Interpolation::readData(const char * path, int delta, int delta2)
{
	if (FILE *pFile = fopen(path, "r")) {
		double x, y;
		while (fscanf(pFile, "%lf\t%lf", &x, &y) > 0)
			points.push_back(Point{ x, y });

		int size = points.size();
		if (delta != 0)
			for (int i = 0; i < size; i += delta)
				pickedPoints.push_back(points[i]);
		else {
			int index = 0;
			while (index < size) {
				pickedPoints.push_back(points[index]);
				double multiplier = (-abs((1.6 / size) * (index - size / 2.0)) + 1.0);
				index += delta2 * multiplier;
			}
			if (index < size - 1)
				pickedPoints.push_back(points.back());
		}

		subintervals = pickedPoints.size() - 1;
	}
}

void Interpolation::writeData(const char * path)
{
	if (FILE *pFile = fopen(path, "w")) {
		for (int i = 0; i < result.size(); ++i)
			fprintf(pFile, "%f\t%f\n", result[i].x, result[i].y);
		fflush(pFile);
	}
}

void Interpolation::clearData()
{
	points.clear();
	result.clear();
	pickedPoints.clear();
}
