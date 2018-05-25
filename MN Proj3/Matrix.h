#pragma once

class Matrix
{
private:
	double *values;
	int width, height, size;

public:
	Matrix();
	Matrix(int h, int w);
	Matrix(int h, int w, int fillValue);
	Matrix(const Matrix& m);

	void zero();

	void print(int n, int m);
	int getWidth();
	int getHeight();
	int getSize();

	double& operator() (int row, int col);
	double operator() (int row, int col) const;
	double& operator() (int elementIndex);
	double operator() (int elementIndex) const;
	Matrix operator+(const Matrix& m);
	Matrix operator-(const Matrix& m);
	Matrix operator*(const Matrix& m);
	Matrix operator*(const int scalar);
	Matrix& operator=(const Matrix& m);

	~Matrix();
};

