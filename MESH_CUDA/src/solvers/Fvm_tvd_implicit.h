/*
 * Fvm_tvd_implicit.h
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#ifndef FVM_TVD_IMPLICIT_H_
#define FVM_TVD_IMPLICIT_H_

#include "../mesh/Mesh.h"
#include "../iterators/MeshIterator.h"
#include "../iterators/FilterIterator.h"
#include "../iterators/BndIterator.h"
#include "Method.h"
#include "MatrixSolver.h"

#include <vector>

using namespace std;

class FVM_TVD_IMPLICIT: public Method
{
public:
		virtual void init(char* xmlFileName);
		virtual void run();
		virtual void done();

protected:
		void save(int);

private:
		static const int PLUS_JACOBIAN = 0;
		static const int MINUS_JACOBIAN = 1;

		double			TMAX;
		int				STEP_MAX;
		double			TAU;
		double			CFL;

private:
        Mesh* msh;
        vector<string> bndInletNames;
        vector<string> bndOutletNames;
        vector<string> bndRigidNames;

        double** allocate_mem();
        void free_mem(double**);
        void matrix_A(double**, double**, double*, double**, int);
        void eigen_values(double*, double, double, double, double, Point);
        void left_eigen_vecs(double**, double, double, double, double, double, Point);
        void right_eigen_vecs(double**, double, double, double, double, double, Point);
        void calc_F(double*, CellFluidDynamicsProps);
        void calc_H(double*, CellFluidDynamicsProps);
        void calc_G(double*, CellFluidDynamicsProps);
        void flux_Lax_Friedrichs(double*, double*, CellFluidDynamicsProps, double*, CellFluidDynamicsProps, Point);
};

#endif /* FVM_TVD_IMPLICIT_H_ */
