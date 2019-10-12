#ifndef SOLVEHEATEQ_H
#define SOLVEHEATEQ_H

#include "../mesh/Mesh.h"
#include "../iterators/MeshIterator.h"
#include "../iterators/FilterIterator.h"
#include "../iterators/BndIterator.h"
#include "../tinyxml/tinyxml2.h"
#include "cuda.h"

class SolveHeatEq
{
    private:
        Mesh* msh;
        vector<string> bndNamesT;
        vector<string> bndNamesFlux;

    public:
        SolveHeatEq(Mesh* m) : msh(m) {}
        ~SolveHeatEq();

        void init(const char *filename);
        void calcHeatEquation(double t_max);
        void save(const char *filename, const char *header);
};

#endif // SOLVEHEATEQ_H
