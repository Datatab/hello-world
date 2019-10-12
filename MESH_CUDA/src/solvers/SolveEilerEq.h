#ifndef SOLVEEILEREQ_H
#define SOLVEEILEREQ_H

#include "../mesh/Mesh.h"
#include "../iterators/MeshIterator.h"
#include "../iterators/FilterIterator.h"
#include "../iterators/BndIterator.h"
#include "../tinyxml/tinyxml2.h"
#include "MatrixSolver.h"

class SolveEilerEq
{
	protected:
		static const int PLUS_JACOBIAN = 0;
		static const int MINUS_JACOBIAN = 1;

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

    public:
       SolveEilerEq(Mesh* m) : msh(m) {}
       ~SolveEilerEq();

       void init(const char *filename);
       void calcEilerEquation(double t_max);
       void save(const char *filename, const char *header);

};

#endif // SOLVEEILEREQ_H
