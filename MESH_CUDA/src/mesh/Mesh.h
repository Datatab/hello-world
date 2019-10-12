#ifndef MESH_H
#define MESH_H

#include <vector>
#include <cmath>
#include <cstdio>
#include <string>
#include <map>
#include "../mesh_properties/BndFaceTemperature.h"
#include "../mesh_properties/CellTemperature.h"
#include "../mesh_properties/CellFluidDynamicsProps.h"


class PointIterator;

template <class T>
class MeshIterator;

template <class Predicate, class Iterator, class T>
class FilterIterator;

template <class T>
class BndIterator;

class Point;
class Edge;
class Face;
class Cell;
/* Callback - функция для iterateCells() */
typedef void(*iterateCellsFunc)(Cell*);
/* Callback - функция для iterateFaces() */
typedef void(*iterateFacesFunc)(Face*);
/* Callback - функция для iterateEdges() */
typedef void(*iterateEdgesFunc)(Edge*);
/* Callback - функция для iteratePoints() */
typedef void(*iteratePointsFunc)(Point*);

using namespace std;

class Point
{
public:
    double x;
    double y;
    double z;

    Point operator+(Point) const;
    Point operator/(double) const;
    static double scalar_product(Point*, Point*);
    static void normalize(Point*);
};

///////////////////////

class Cell;

class Edge // Класс некорректно работает
{
private:
    Point* p[2];

public:
    Edge(Point*, Point*);
    double getlength();
    static double getlength(Point*, Point*);
    static Point getVector(Point*, Point*);

    friend class MeshReaderUnv;
};

////////////////////////

class Face
{
 private:
    int type;
    Point** p;
    Edge** e;

    Cell* c[2]; // указатель на смежный Cell

    Point n; // единичный вектор нормали
    int in_cell; // 0 или 1, Выяснение направление вектора нормали по отношению к соседним клеткам
    int out_cell;
    int pCount;
    int eCount;
    int countOfUsing;

    double S;

    double h; // растояние между ячейками

    BndFaceTemperature bndT;

    CellFluidDynamicsProps faceFDP;

 public:
    Face(const int&, Point*, Point*, Point*);
    Face(const int&, Point*, Point*, Point*, Point*);

    ~Face();

    void initFace(Cell*, Edge*, Edge*, Edge*);
    void initFace(Cell*, Edge*, Edge*, Edge*, Edge*);

    void area();
    void calc_h();

    Cell** getCells()
    {
        return c;
    }

    friend class Cell;
    friend class Mesh;
    friend class MeshReaderUnv;
    friend class SolveHeatEq;
    friend class SolveEilerEq;
    friend class FVM_TVD_IMPLICIT;
};

///////////////////

class Mesh;

class Cell {

 private:
	int index; // индекс в массиве Mesh::cells
    int type;
    Point** p;
    Edge** e;
    Face** f;
    int pCount;
    int eCount;
    int fCount;
    double V;
    Point center;


    CellTemperature cellT; // Температура в клетке

    CellFluidDynamicsProps cellFDP; // Свойства газа в клетке

 public:
    Cell(const int&, Point*, Point*, Point*, Point*);
    Cell(const int&, Point*, Point*, Point*, Point*, Point*, Point*);
    Cell(const int&, Point*, Point*, Point*, Point*, Point*, Point*, Point*, Point*);
    ~Cell();

    void volume();

    friend class Face;
    friend class Mesh;
    friend class MeshReaderUnv;
    friend class SolveHeatEq;
    friend class SolveEilerEq;
    friend class FVM_TVD_IMPLICIT;
    //friend class Method;
};

//////////////

class IsBoundaryFace {
public:
    bool operator()(Face* f);
};

class IsInnerFace {
public:
    bool operator()(Face* f);
};

//////////////

class Mesh
{
private:
    Point* points;
    unsigned int pCount;
    vector<Edge*> edges;
    vector<Face*> faces;
    vector<Cell*> cells;
    vector<Face*> inner_faces;
    map<string, vector<Face*> > bnd_faces;

public:
    Mesh();
    ~Mesh();

    void createPoints(Point*, unsigned int);
    void createPoints(const vector<Point>&);

    friend class MeshReaderUnv;
    friend class SolveHeatEq;
    friend class SolveEilerEq;
    friend class FVM_TVD_IMPLICIT;

    typedef MeshIterator<Face> FaceIterator;
    typedef MeshIterator<Edge> EdgeIterator;
    typedef MeshIterator<Cell> CellIterator;
    typedef FilterIterator<IsBoundaryFace, FaceIterator, Face> BoundaryFaceIterator;
    typedef FilterIterator<IsInnerFace, FaceIterator, Face> InnerFaceIterator;
    typedef BndIterator<Face> BndFaceIterator;

    CellIterator beginCell();
    CellIterator endCell();

    FaceIterator beginFace();
    FaceIterator endFace();

    EdgeIterator beginEdge();
    EdgeIterator endEdge();

    PointIterator beginPoint();
    PointIterator endPoint();

    BoundaryFaceIterator beginBoundaryFace();
    BoundaryFaceIterator endBoundaryFace();
    ////
    FaceIterator beginInnerFace();
    FaceIterator endInnerFace();
    ////
    /*
    InnerFaceIterator beginInnerFace();
    InnerFaceIterator endInnerFace();
	*/
    FaceIterator beginBndFace(string);
    FaceIterator endBndFace(string);
    BndFaceIterator beginBndFace(map<std::string, vector<Face*>> *m, vector<string> *str);
    BndFaceIterator endBndFace(map<std::string, vector<Face*>> *m, vector<string> *str);

    void iterateCells(iterateCellsFunc);
    void iterateFaces(iterateFacesFunc);
    void iterateEdges(iterateEdgesFunc);
    void iteratePoints(iteratePointsFunc);
};

#endif // MESH_H
