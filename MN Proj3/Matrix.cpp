#include "Matrix.h"
#include <iostream>
#include <cstdio>

using namespace std;

Matrix::Matrix()
	: values(nullptr), width(0), height(0), size(0) {}

Matrix::Matrix(int h, int w)
	: width(w), height(h), size(w*h)
{
	values = new double[size];
}

Matrix::Matrix(int h, int w, int fillValue)
	: Matrix(h, w)
{
	for (int i = 0; i < size; ++i)
		values[i] = fillValue;
}

Matrix::Matrix(const Matrix & m)
{
	width = m.width;
	height = m.height;
	size = m.size;
	values = new double[size];
	for (int i = 0; i < size; ++i)
		values[i] = m.values[i];
}

void Matrix::zero()
{
	for (int i = 0; i < size; ++i)
		values[i] = 0;
}

void Matrix::print(int n, int m)
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j)
			printf(" % .2f", (*this)(i, j));
		printf("\n");
	}
}

int Matrix::getWidth()
{
	return width;
}

int Matrix::getHeight()
{
	return height;
}

int Matrix::getSize()
{
	return size;
}

double& Matrix::operator()(int row, int col)
{
	return values[row*width + col];
}

double Matrix::operator()(int row, int col) const
{
	return values[row*width + col];
}

double& Matrix::operator()(int elementIndex)
{
	return values[elementIndex];
}

double Matrix::operator()(int elementIndex) const
{
	return values[elementIndex];
}

Matrix Matrix::operator+(const Matrix & m)
{
	Matrix matrix(height, width);
	for (int i = 0; i < size; ++i)
		matrix.values[i] = values[i] + m.values[i];
	return matrix;
}

Matrix Matrix::operator-(const Matrix & m)
{
	Matrix matrix(height, width);
	for (int i = 0; i < size; ++i)
		matrix.values[i] = values[i] - m.values[i];
	return matrix;
}

Matrix Matrix::operator*(const Matrix & m)
{
	Matrix matrix(height, m.width);
	double value = 0;
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < m.width; ++j) {
			for (int k = 0; k < width; ++k)
				value += (*this)(i, k) * m(k, j);
			matrix(i, j) = value;
			value = 0;
		}
	return matrix;
}

Matrix Matrix::operator*(const int scalar)
{
	Matrix matrix(height, width);
	for (int i = 0; i < size; ++i)
		matrix.values[i] = values[i] * scalar;
	return matrix;
}

Matrix& Matrix::operator=(const Matrix & m)
{
	if (this != &m) {
		if (width != m.width || height != m.height) {
			if (values != nullptr)
				delete[] values;
			width = m.width;
			height = m.height;
			size = m.size;
			values = new double[size];
		}
		for (int i = 0; i < size; ++i)
			values[i] = m.values[i];
	}
	return *this;
}

Matrix::~Matrix()
{
	delete[] values;
}
