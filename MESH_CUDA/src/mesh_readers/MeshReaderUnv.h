#ifndef MESHREADERUNV_H
#define MESHREADERUNV_H

#include "../mesh/Mesh.h"
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <cstring>
//#include <iostream>

using namespace std;

class MeshReaderUnv
{
    private:
        vector<set<Edge*> > pEdge;
        vector<map<int, set<Face*> > > pFace;
        vector<vector<int> > cells;
        vector<int> type_cells;
        map<int, Face*> bnd_cond;

        void read_block(vector<string>&, ifstream&);
        void parse_block(Mesh*, vector<string>&);
        void parse_block_164(Mesh*, vector<string>&);
        void parse_block_2420(Mesh*, vector<string>&);
        void parse_block_2411(Mesh*, vector<string>&);
        void parse_block_2412(Mesh*, vector<string>&);
        void parse_block_2477(Mesh*, vector<string>&);

        Face* find_face(const int&, const int&, const int&, const int&);
        Face* find_face(const int&, const int&, const int&, const int&, const int&);
        Edge* find_edge(const int&, const int&);

    public:
        MeshReaderUnv() { }
        ~MeshReaderUnv() { }

        virtual void read(Mesh* mesh, const char* fileName);
};

#endif // MESHREADERUNV_H
