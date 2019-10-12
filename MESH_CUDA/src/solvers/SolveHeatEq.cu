#include "SolveHeatEq.h"



SolveHeatEq::~SolveHeatEq()
{
    //dtor
}

void SolveHeatEq::init(const char *filename)
{
    tinyxml2::XMLDocument xmlDoc;
    tinyxml2::XMLError eResult = xmlDoc.LoadFile(filename);
    if (eResult != tinyxml2::XML_SUCCESS)
    {
        printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
        exit(eResult);
    }

    /* Температура для всех cell */
    tinyxml2::XMLElement* pNode1 = xmlDoc.FirstChildElement("regions");
    if (pNode1 == nullptr)
    {
        printf("XMLERROR: No such FirstChildElement.\n");
        exit(1);
    }
    tinyxml2::XMLElement* pNode2 = pNode1->FirstChildElement("region");
    if (pNode2 == nullptr)
    {
        printf("XMLERROR: No such FirstChildElement.\n");
        exit(1);
    }
    tinyxml2::XMLElement* pNode3 = pNode2->FirstChildElement("parameters");
    if (pNode3 == nullptr)
    {
        printf("XMLERROR: No such FirstChildElement.\n");
        exit(1);
    }
    tinyxml2::XMLElement* pNode4 = pNode3->FirstChildElement("T");
    if (pNode4 == nullptr)
    {
        printf("XMLERROR: No such FirstChildElement.\n");
        exit(1);
    }
    double t;
    double k;

    eResult = pNode4->QueryDoubleAttribute("value", &t);
    if (eResult != tinyxml2::XML_SUCCESS)
        {
            printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
            exit(eResult);
        }

    tinyxml2::XMLElement* pNode5 = pNode3->FirstChildElement("k");
    eResult = pNode5->QueryDoubleAttribute("value", &k);
    if (eResult != tinyxml2::XML_SUCCESS)
    {
            printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
            exit(eResult);
    }
    for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
    {
        it->cellT.T = t;
        it->cellT.k = k;
    }
    /* Температура для заданных face */
    tinyxml2::XMLElement* pBnd = xmlDoc.FirstChildElement("boundaries");
    if (pBnd == nullptr)
    {
        printf("XMLERROR: No such FirstChildElement.\n");
        exit(1);
    }
    tinyxml2::XMLElement* pBndElement = pBnd->FirstChildElement("boundCond");
    if (pBndElement == nullptr)
    {
        printf("XMLERROR: No such FirstChildElement.\n");
        exit(1);
    }
    while (pBndElement != nullptr)
    {
        string str;
        tinyxml2::XMLElement* pName = pBndElement->FirstChildElement("name");
        if (pName == nullptr)
        {
            printf("XMLERROR: No such FirstChildElement.\n");
            exit(1);
        }
        str = pName->GetText();

        tinyxml2::XMLElement* pPar = pBndElement->FirstChildElement("parameters");
        if (pPar == nullptr)
        {
            printf("XMLERROR: No such FirstChildElement.\n");
            exit(1);
        }
        tinyxml2::XMLElement* pT = pPar->FirstChildElement("T");
        tinyxml2::XMLElement* pFlux = pPar->FirstChildElement("Flux");
        if (pT == nullptr && pFlux == nullptr)
        {
            printf("XMLERROR: No such FirstChildElement.\n");
            exit(1);
        }
        double temp;

        if(pT != nullptr)
        {
        	eResult = pT->QueryDoubleAttribute("value", &temp);
        	if (eResult != tinyxml2::XML_SUCCESS)
        	{
        		printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
                exit(eResult);
            }

        	for (Mesh::FaceIterator it = msh->beginBndFace(str), ite = msh->endBndFace(str); it != ite; ++it)
        	{
        		it->bndT.T = temp;
        	}
        	bndNamesT.push_back(str);
        }
        else if(pFlux != nullptr)
        {
        	eResult = pFlux->QueryDoubleAttribute("value", &temp);
        	if (eResult != tinyxml2::XML_SUCCESS)
        	{
        		printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
                exit(eResult);
            }

        	for (Mesh::FaceIterator it = msh->beginBndFace(str), ite = msh->endBndFace(str); it != ite; ++it)
        	{
        		it->bndT.Flux = temp;
        	}
        	bndNamesFlux.push_back(str);
        }

        pBndElement = pBndElement->NextSiblingElement("boundCond");
    }
}



__constant__ int const_sizeBndFlux;
__constant__ int const_sizeBndT;
__constant__ int const_sizeInnerFaces;
__constant__ int const_sizeCell;
__constant__ double const_k;
__constant__ double const_tau;


 __global__ void BndFacesFlux( double* cell_Flux, double* flux, int* indCell) {

	 int tid = threadIdx.x + blockIdx.x * blockDim.x;
	 	
	 while(tid < const_sizeBndFlux)
	 {
		 cell_Flux[ indCell[tid] ] += flux[tid];
		
		 tid += blockDim.x * gridDim.x;
	 }
}


 

 __global__ void BndFacesT( double* T, double* cell_Flux, double* h, double* S, double* bnd_T, int* indCell) {

	 int tid = threadIdx.x + blockIdx.x * blockDim.x;
	 
	 while(tid < const_sizeBndT)
	 {
		 cell_Flux[ indCell[tid] ] += S[tid] * (bnd_T[tid] - T[ indCell[tid] ]) / h[tid];
		 
		 tid += blockDim.x * gridDim.x;
	 }
 }


 

__global__ void InnerFaces( double* T, double* cell_Flux, double* h, double* S, int* indCell) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	double temp_value;

	while(tid < const_sizeInnerFaces ) {
		 temp_value = S[tid] * (T[ indCell[2 * tid + 1] ] - T[ indCell[2 * tid] ]) / h[tid];
		 cell_Flux[ indCell[2 * tid] ] += temp_value;
		 cell_Flux[ indCell[2 * tid + 1] ] -= temp_value;
		 tid += blockDim.x * gridDim.x;
	 }
}


__global__ void Cells( double* T, double* cell_Flux, double* V) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	
	while(tid < const_sizeCell) {
		 T[tid] += const_k * cell_Flux[tid] * const_tau / V[tid];
		 cell_Flux[tid] = 0;
		 tid += blockDim.x * gridDim.x;
	 }
}




void SolveHeatEq::calcHeatEquation(double t_max)
{
    double min_volume = msh->cells[0]->V;

    for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
    {
        double vol = it->V;
        if (min_volume > vol)
            min_volume = vol;
    }

    double k = msh->cells[0]->cellT.k;

    double tau = min_volume / (2.1 * k);

    double temp_value;
    double t = 0;

	// Cell'ы
    int sizeCell = msh->cells.size();

    double* T = (double*)malloc( sizeCell * sizeof(double) ); // Температура в Cell'е
    double* cell_Flux = (double*)malloc( sizeCell * sizeof(double) ); // Суммарный поток в Cell'е
    double* V = (double*)malloc( sizeCell * sizeof(double) ); // Объем Cell'а

    for(int i = 0; i < sizeCell; i++)
    {
    	T[i] = msh->cells[i]->cellT.T;
    	V[i] = msh->cells[i]->V;
    }

    memset(cell_Flux, 0, sizeCell * sizeof(double) );

    double* dev_T;
    cudaMalloc( (void**)&dev_T, sizeCell * sizeof(double) );
    cudaMemcpy( dev_T, T, sizeCell * sizeof(double), cudaMemcpyHostToDevice );

    double* dev_V;
    cudaMalloc( (void**)&dev_V, sizeCell * sizeof(double) );
    cudaMemcpy( dev_V, V, sizeCell * sizeof(double), cudaMemcpyHostToDevice );

    double* dev_cell_Flux;
    cudaMalloc( (void**)&dev_cell_Flux, sizeCell * sizeof(double) );
    cudaMemcpy( dev_cell_Flux, cell_Flux, sizeCell * sizeof(double), cudaMemcpyHostToDevice );




	// внутренние face'ы
    int sizeInnerFaces = msh->inner_faces.size();

    double* inn_faces_h = (double*)malloc( sizeInnerFaces * sizeof(double) ); // расстояние между ячейками по данному Face'у
    double* inn_faces_S = (double*)malloc( sizeInnerFaces * sizeof(double) ); // Площадь
    int* inn_faces_indCell = (int*)malloc( 2 * sizeInnerFaces * sizeof(int) ); // inn_faces_indCell[ 2*i ], inn_faces_indCell[ 2*i + 1 ] - индексы смежных Cell'ов i-ого face'a в массивах T, cell_Flux, V
   
    for(int i = 0; i < sizeInnerFaces; i++)
    {
    	inn_faces_h[i] = msh->inner_faces[i]->h;
    	inn_faces_S[i] = msh->inner_faces[i]->S;
    	inn_faces_indCell[ 2*i ] = msh->inner_faces[i]->c[0]->index;
    	inn_faces_indCell[ 2*i + 1 ] = msh->inner_faces[i]->c[1]->index;
    }

    double* dev_inn_faces_h;
    cudaMalloc( (void**)&dev_inn_faces_h,  sizeInnerFaces * sizeof(double));
    cudaMemcpy( dev_inn_faces_h, inn_faces_h, sizeInnerFaces * sizeof(double), cudaMemcpyHostToDevice );

    double* dev_inn_faces_S;
    cudaMalloc( (void**)&dev_inn_faces_S,  sizeInnerFaces * sizeof(double));
    cudaMemcpy( dev_inn_faces_S, inn_faces_S, sizeInnerFaces * sizeof(double), cudaMemcpyHostToDevice );

    int* dev_inn_faces_indCell;
    cudaMalloc( (void**)&dev_inn_faces_indCell, 2 * sizeInnerFaces * sizeof(int));
    cudaMemcpy( dev_inn_faces_indCell, inn_faces_indCell, 2 * sizeInnerFaces * sizeof(int), cudaMemcpyHostToDevice );

   

	// поток на границе
    int sizeBndFlux = 0; 
    for(vector<string>::iterator it = bndNamesFlux.begin(); it != bndNamesFlux.end(); ++it)
	{
		sizeBndFlux += msh->bnd_faces[(*it)].size();
	}


	double* flux_bnd_h = (double*)malloc( sizeBndFlux * sizeof(double) );
	double* flux_bnd_S = (double*)malloc( sizeBndFlux * sizeof(double) );
	double* bnd_flux = (double*)malloc( sizeBndFlux * sizeof(double) ); // поток
	int* flux_bnd_indCell = (int*)malloc( sizeBndFlux * sizeof(int) );  // flux_bnd_indCell[ i ] - индекс смежного Cell'a для i-ого face'a 

	int q = 0;
	for(vector<string>::iterator it = bndNamesFlux.begin(); it != bndNamesFlux.end(); ++it)
	{
		vector<Face*> vec_temp = msh->bnd_faces[(*it)];
		for(int i = 0; i < vec_temp.size(); i++)
		{
			flux_bnd_h[q] = vec_temp[i]->h;
			flux_bnd_S[q] = vec_temp[i]->S;
			bnd_flux[q] = vec_temp[i]->bndT.Flux;
			flux_bnd_indCell[q] = vec_temp[i]->c[0]->index;
			q++;
		}
	}

	double* dev_flux_bnd_h;
	cudaMalloc( (void**)&dev_flux_bnd_h,  sizeBndFlux * sizeof(double));
	cudaMemcpy( dev_flux_bnd_h, flux_bnd_h, sizeBndFlux * sizeof(double), cudaMemcpyHostToDevice );

	double* dev_flux_bnd_S;
	cudaMalloc( (void**)&dev_flux_bnd_S,  sizeBndFlux * sizeof(double));
	cudaMemcpy( dev_flux_bnd_S, flux_bnd_S, sizeBndFlux * sizeof(double), cudaMemcpyHostToDevice );

	double* dev_bnd_flux;
	cudaMalloc( (void**)&dev_bnd_flux,  sizeBndFlux * sizeof(double));
	cudaMemcpy( dev_bnd_flux, bnd_flux, sizeBndFlux * sizeof(double), cudaMemcpyHostToDevice );

	int* dev_flux_bnd_indCell;
	cudaMalloc( (void**)&dev_flux_bnd_indCell, sizeBndFlux * sizeof(int) );
	cudaMemcpy( dev_flux_bnd_indCell, flux_bnd_indCell, sizeBndFlux * sizeof(int), cudaMemcpyHostToDevice );



	// Температура на границе
	int sizeBndT = 0;
	for(vector<string>::iterator it = bndNamesT.begin(); it != bndNamesT.end(); ++it)
	{
		sizeBndT += msh->bnd_faces[(*it)].size();
	}


	double* T_bnd_h = (double*)malloc( sizeBndT * sizeof(double) ); // расстояние до фиктивной ячейки
	double* T_bnd_S = (double*)malloc( sizeBndT * sizeof(double) );
	double* bnd_T = (double*)malloc( sizeBndT * sizeof(double) ); // температура face'a
	int* T_bnd_indCell = (int*)malloc( sizeBndT * sizeof(int) ); // индекс Cell'a

	q = 0;
	for(vector<string>::iterator it = bndNamesT.begin(); it != bndNamesT.end(); ++it)
	{
		vector<Face*> vec_temp = msh->bnd_faces[(*it)];
		for(int i = 0; i < vec_temp.size(); i++)
		{
			T_bnd_h[q] = vec_temp[i]->h;
			T_bnd_S[q] = vec_temp[i]->S;
			bnd_T[q] = vec_temp[i]->bndT.T;
			T_bnd_indCell[q] = vec_temp[i]->c[0]->index;
			q++;
		}
	}

	double* dev_T_bnd_h;
	cudaMalloc( (void**)&dev_T_bnd_h,  sizeBndT * sizeof(double) );
	cudaMemcpy( dev_T_bnd_h, T_bnd_h, sizeBndT * sizeof(double), cudaMemcpyHostToDevice );

	double* dev_T_bnd_S;
	cudaMalloc( (void**)&dev_T_bnd_S,  sizeBndT * sizeof(double) );
	cudaMemcpy( dev_T_bnd_S, T_bnd_S, sizeBndT * sizeof(double), cudaMemcpyHostToDevice );

	double* dev_bnd_T;
	cudaMalloc( (void**)&dev_bnd_T,  sizeBndT * sizeof(double) );
	cudaMemcpy( dev_bnd_T, bnd_T, sizeBndT * sizeof(double), cudaMemcpyHostToDevice );

	int* dev_T_bnd_indCell;
	cudaMalloc( (void**)&dev_T_bnd_indCell, sizeBndT * sizeof(int));
	cudaMemcpy( dev_T_bnd_indCell, T_bnd_indCell, sizeBndT * sizeof(int), cudaMemcpyHostToDevice );

	
	cudaMemcpyToSymbol(const_k, &k, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(const_tau, &tau, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(const_sizeCell, &sizeCell, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(const_sizeInnerFaces, &sizeInnerFaces, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(const_sizeBndFlux, &sizeBndFlux, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(const_sizeBndT, &sizeBndT, sizeof(int), 0, cudaMemcpyHostToDevice);
	

	int threads = 128;
	int blocks = 128;

    while (t < t_max)
    {	
        t += tau;
		
		// Поток на границе
		BndFacesFlux <<< blocks, threads >>>(dev_cell_Flux, dev_bnd_flux, dev_flux_bnd_indCell);
		// Поток на границе
		BndFacesT <<< blocks, threads >>>(dev_T, dev_cell_Flux, dev_T_bnd_h, dev_T_bnd_S, dev_bnd_T, dev_T_bnd_indCell);
		// Поток во внутренних face'ах
		InnerFaces <<< blocks, threads >>>(dev_T, dev_cell_Flux, dev_inn_faces_h, dev_inn_faces_S, dev_inn_faces_indCell);
		// Новое значение темрературы в ячейке
		Cells <<< blocks, threads >>>(dev_T, dev_cell_Flux, dev_V);
		
		/* Решение на CPU
        for (Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndNamesFlux), ite = msh->endBndFace(&(msh->bnd_faces), &bndNamesFlux); it != ite; ++it)
        {
           	it->c[0]->cellT.Flux += it->bndT.Flux;
        }
		 
        
        for (Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndNamesT), ite = msh->endBndFace(&(msh->bnd_faces), &bndNamesT); it != ite; ++it)
        {
        	it->c[0]->cellT.Flux += (it->S*(it->bndT.T - it->c[0]->cellT.T) / it->h);

        }
        
        for (Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite; ++it)
        {
        	temp_value = it->S*(it->c[1]->cellT.T - it->c[0]->cellT.T) / it->h;
        	it->c[0]->cellT.Flux += temp_value;
        	it->c[1]->cellT.Flux -= temp_value;

        }
		
        for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
        {
        	it->cellT.T += it->cellT.Flux * tau * k / it->V;
        	it->cellT.Flux = 0;
        }
   		 */
    }

	cudaMemcpy(T, dev_T, sizeCell * sizeof(double), cudaMemcpyDeviceToHost);
	for (int i = 0; i < sizeCell; i++)
	{
		msh->cells[i]->cellT.T = T[i];
	}

}


void SolveHeatEq::save(const char *filename, const char *header)
{
    FILE *out;
    out = fopen(filename, "w");
    fprintf(out, "# vtk DataFile Version 3.0\n");
    //The header can be used to describe the data
    fprintf(out, "%s\n", header);
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(out, "POINTS %d double\n", msh->pCount);
    for (int i = 0; i < msh->pCount; i++)
    {
        fprintf(out, "%f %f %f\n", msh->points[i].x, msh->points[i].y, msh->points[i].z);
    }

    int cellCount = msh->cells.size();

    int cellSize = 0;//the size of the cell list (count of points in all cells)
    for (int i = 0; i < cellCount; i++)
    {
        cellSize += msh->cells[i]->pCount;
    }

    /*
    cellSize + cellCount :
    cellSize + one number for each cell - count of points in this cell
    */
    fprintf(out, "CELLS %d %d\n", cellCount, cellSize + cellCount);
    for (int i = 0; i < cellCount; i++)
    {
        fprintf(out, "%d", msh->cells[i]->pCount);
        for (int k = 0; k < msh->cells[i]->pCount; k++)
        {
            int ind = -1;
            Point* addr = msh->cells[i]->p[k];
            for (int j = 0; j < msh->pCount; j++)
            {
                if (&(msh->points[j]) == addr)
                    ind = j;
            }
            fprintf(out, " %d", ind);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "CELL_TYPES %d\n", cellCount);
    for (int i = 0; i < cellCount; i++)
    {
        switch (msh->cells[i]->type)
        {
        case 111:
        {
            fprintf(out, "10\n"); //10 - VTK_TETRA
            break;
        }
        case 112:
        {
            fprintf(out, "13\n"); //13 - VTK_WEDGE
            break;
        }
        case 115:
        {
            fprintf(out, "12\n"); //12 - VTK_HEXAHEDRON
            break;
        }
        }

    }
    fprintf(out, "CELL_DATA %d\nSCALARS temperature double 1\nLOOKUP_TABLE default\n", cellCount);
    for (int i = 0; i < cellCount; i++)
    {
        fprintf(out, "%f\n", msh->cells[i]->cellT.T);
    }
    fclose(out);
}


