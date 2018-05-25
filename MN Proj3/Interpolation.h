#pragma once
#include <vector>
#include "MatricesOperations.h"
#include "Matrix.h"

using namespace std;

struct Point {
	double x, y;
};

class Interpolation
{
public:
	Interpolation();

	int lagrange(const char *inputFile, const char *outputFile, int delta, int delta2);
	int spline3Deg(const char *inputFile, const char *outputFile, int delta, int delta2);

	double fiFunc(double x, int index);
	double interpolationLagrangeFunc(double x);
	void generateSplineFuncs();
	double splineCalcValue(double x);

	~Interpolation();
private:
	vector<Point> points;
	vector<Point> result;
	vector<Point> pickedPoints;
	MatricesOperations matOp;
	Matrix resultMat;
	int subintervals;

	void readData(const char *path, int delta, int delta2);
	void writeData(const char *path);
	void clearData();
};
