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

	int lagrange(const char *inputFile, const char *outputFile, int delta);
	int spline3Deg(const char *inputFile, const char *outputFile, int delta);

	double fiFunc(double x, int index, int delta);
	double interpolationLagrangeFunc(double x, int delta);
	void generateSplineFuncs(int delta);
	double splineCalcValue(double x, int delta);

	~Interpolation();
private:
	vector<Point> points;
	vector<Point> result;
	vector<Point> pickedPoints;
	MatricesOperations matOp;
	Matrix resultMat;
	//int delta;

	void readData(const char *path);
	void writeData(const char *path);
	void clearData();
};
