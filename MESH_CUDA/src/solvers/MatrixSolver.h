/*
 * MatrixSolver.h
 *
 *  Created on: Oct 7, 2019
 *      Author: v1
 */

#ifndef MATRIXSOLVER_H_
#define MATRIXSOLVER_H_

#include <cmath>

class MatrixSolver {

public:
	MatrixSolver(unsigned int);
	virtual ~MatrixSolver();

	void zeidel_method(double, unsigned int);
	void clear();

	void set_matrix(unsigned int, unsigned int, double);
	void set_f(unsigned int, double);

	void add_matrix(unsigned int, unsigned int, double);
	void add_f(unsigned int, double);
	void subtract_matrix(unsigned int, unsigned int, double);
	void subtract_f(unsigned int, double);

	double* get_solution();
	int get_size();

private:
	double** A;
	double* f;
	double* x;
	double* temp_x;
	int N;

	double norm_2(double*, double*);
};

#endif /* MATRIXSOLVER_H_ */
