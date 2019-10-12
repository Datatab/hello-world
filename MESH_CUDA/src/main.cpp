#include "mesh/Mesh.h"
#include "mesh_readers/MeshReaderUnv.h"
#include "solvers/SolveHeatEq.h"
#include "solvers/SolveEilerEq.h"
//#include "solvers/Solver.h"
#include <ctime>



int main()
{

    Mesh* msh = new Mesh();

    MeshReaderUnv mru;
    //mru.read(msh, "Kluch.unv");
    mru.read(msh, "Cube5_ForFDP.unv");

    SolveEilerEq* sheq = new SolveEilerEq(msh);

    sheq->init("bndFDPCube.xml");

    double te = clock();
    sheq->calcEilerEquation(0.008);
    cout << endl << "time of execution : " << (clock() - te) / CLOCKS_PER_SEC << endl;
    sheq->save("exampleXMLCuda.vtk", "ourMesh");
    // 463.062 - 0.1 сек


    //Method* m = Solver::initMethod("bndFDPCube.xml");


    return 0;
}
