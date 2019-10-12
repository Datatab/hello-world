/*
 * MatrixSolver.cpp
 *
 *  Created on: Oct 7, 2019
 *      Author: v1
 */

#include "MatrixSolver.h"

MatrixSolver::MatrixSolver(unsigned int N) {

	this->N = N;

	f = new double[N];
	x = new double[N];
	temp_x = new double[N];

	A = new double*[N];
	for(int i = 0; i < N; i++)
		A[i] = new double[N];
}

MatrixSolver::~MatrixSolver() {

	delete [] f;
	delete [] x;
	delete [] temp_x;

	for(int i = 0; i < N; i++)
	{
		delete [] A[i];
	}

	delete [] A;
}

int MatrixSolver::get_size()
{
	return N;
}

double* MatrixSolver::get_solution()
{
	return x;
}

void MatrixSolver::zeidel_method(double eps, unsigned int max_iter)
{
	for(int i = 0; i < N; i++)
	{
		temp_x[i] = f[i];
	}

	int k = 0;
	double norm;
	do
	{
		for(int i = 0; i < N; i++)
		{
			double S = 0;

			for(int j = 0; j < i; j++)
				S += A[i][j] * x[j];

			for(int j = i + 1; j < N; j++)
				S += A[i][j] * temp_x[j];

			x[i] = (f[i] - S) / A[i][i];
		}

		norm = norm_2(x, temp_x);

		for(int i = 0; i < N; i++)
		{
			temp_x[i] = x[i];
		}

		k++;
	}while( norm >= eps && k < max_iter);
}

double MatrixSolver::norm_2(double* u1, double* u2)
{
	double S = 0;
	for(int i = 0; i < N; i++)
	{
		S += (u1[i] - u2[i]) * (u1[i] - u2[i]);
	}

	return sqrt(S);
}

void MatrixSolver::clear()
{
	for(int i = 0; i < N; i++)
	{
		f[i] = 0;
	}

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			A[i][j] = 0;
		}
	}
}

void MatrixSolver::set_matrix(unsigned int i, unsigned int j, double val)
{
	A[i][j] = val;
}

void MatrixSolver::set_f(unsigned int i, double val)
{
	f[i] = val;
}

void MatrixSolver::add_matrix(unsigned int i, unsigned int j, double val)
{
	A[i][j] += val;
}

void MatrixSolver::subtract_matrix(unsigned int i, unsigned int j, double val)
{
	A[i][j] -= val;
}

void MatrixSolver::add_f(unsigned int i, double val)
{
	f[i] += val;
}

void MatrixSolver::subtract_f(unsigned int i, double val)
{
	f[i] -= val;
}


