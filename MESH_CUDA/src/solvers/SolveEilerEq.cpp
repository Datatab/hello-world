#include "SolveEilerEq.h"

SolveEilerEq::~SolveEilerEq()
{
    //dtor
}


void SolveEilerEq::init(const char *filename)
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
    // ro
    tinyxml2::XMLElement* pNode4 = pNode3->FirstChildElement("ro");
    if (pNode4 == nullptr)
    {
        printf("XMLERROR: No such FirstChildElement.\n");
        exit(1);
    }

    double ro;
    eResult = pNode4->QueryDoubleAttribute("value", &ro);
    if (eResult != tinyxml2::XML_SUCCESS)
    {
        printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
        exit(eResult);
    }

    // u
    pNode4 = pNode3->FirstChildElement("u");
	if (pNode4 == nullptr)
	{
		printf("XMLERROR: No such FirstChildElement.\n");
		exit(1);
	}

	double u;
	eResult = pNode4->QueryDoubleAttribute("value", &u);
	if (eResult != tinyxml2::XML_SUCCESS)
	{
		printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
		exit(eResult);
	}

	// v
	pNode4 = pNode3->FirstChildElement("v");
	if (pNode4 == nullptr)
	{
		printf("XMLERROR: No such FirstChildElement.\n");
		exit(1);
	}

	double v;
	eResult = pNode4->QueryDoubleAttribute("value", &v);
	if (eResult != tinyxml2::XML_SUCCESS)
	{
		printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
		exit(eResult);
	}

	// w
	pNode4 = pNode3->FirstChildElement("w");
	if (pNode4 == nullptr)
	{
		printf("XMLERROR: No such FirstChildElement.\n");
		exit(1);
	}

	double w;
	eResult = pNode4->QueryDoubleAttribute("value", &w);
	if (eResult != tinyxml2::XML_SUCCESS)
	{
		printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
		exit(eResult);
	}

	// P
	pNode4 = pNode3->FirstChildElement("P");
	if (pNode4 == nullptr)
	{
		printf("XMLERROR: No such FirstChildElement.\n");
		exit(1);
	}

	double P;
	eResult = pNode4->QueryDoubleAttribute("value", &P);
	if (eResult != tinyxml2::XML_SUCCESS)
	{
		printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
		exit(eResult);
	}

    for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
    {
        it->cellFDP.ro = ro;
        it->cellFDP.ru = u * ro;
        it->cellFDP.rv = v * ro;
        it->cellFDP.rw = w * ro;
        it->cellFDP.gamma = 1.4;
        it->cellFDP.P = P;

        it->cellFDP.rE = CellFluidDynamicsProps::calc_rE(ro, P, u, v, w, it->cellFDP.gamma);
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
        string name;
        tinyxml2::XMLElement* pName = pBndElement->FirstChildElement("name");
        if (pName == nullptr)
        {
            printf("XMLERROR: No such FirstChildElement.\n");
            exit(1);
        }
        name = pName->GetText();

        string type;
		tinyxml2::XMLElement* pType = pBndElement->FirstChildElement("type");
		if (pType == nullptr)
		{
			printf("XMLERROR: No such FirstChildElement.\n");
			exit(1);
		}
		type = pType->GetText();

        tinyxml2::XMLElement* pPar = pBndElement->FirstChildElement("parameters");
        if (pPar == nullptr)
        {
            printf("XMLERROR: No such FirstChildElement.\n");
            exit(1);
        }

        if( type == "INLET" )
        {
        	// ro
			tinyxml2::XMLElement* pNodeBnd = pPar->FirstChildElement("ro");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &ro);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			// u
			pNodeBnd = pPar->FirstChildElement("u");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &u);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			// v
			pNodeBnd = pPar->FirstChildElement("v");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &v);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			// w
			pNodeBnd = pPar->FirstChildElement("w");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &w);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			// P
			pNodeBnd = pPar->FirstChildElement("P");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &P);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			for (Mesh::FaceIterator it = msh->beginBndFace(name), ite = msh->endBndFace(name); it != ite; ++it)
			{
				it->faceFDP.ro = ro;
				it->faceFDP.ru = u * ro;
				it->faceFDP.rv = v * ro;
				it->faceFDP.rw = w * ro;
				it->faceFDP.gamma = 1.4;
				it->faceFDP.P = P;

				it->faceFDP.rE = CellFluidDynamicsProps::calc_rE(ro, P, u, v, w, it->faceFDP.gamma);
			}

			bndInletNames.push_back(name);
        }
        else if( type == "OUTLET" )
        {
        	// ro
			tinyxml2::XMLElement* pNodeBnd = pPar->FirstChildElement("ro");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &ro);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			// u
			pNodeBnd = pPar->FirstChildElement("u");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &u);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			// v
			pNodeBnd = pPar->FirstChildElement("v");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &v);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			// w
			pNodeBnd = pPar->FirstChildElement("w");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &w);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			// P
			pNodeBnd = pPar->FirstChildElement("P");
			if (pNodeBnd == nullptr)
			{
				printf("XMLERROR: No such FirstChildElement.\n");
				exit(1);
			}

			eResult = pNodeBnd->QueryDoubleAttribute("value", &P);
			if (eResult != tinyxml2::XML_SUCCESS)
			{
				printf("XMLERROR is %d\nXML loading unsuccessfull.\n", eResult);
				exit(eResult);
			}

			for (Mesh::FaceIterator it = msh->beginBndFace(name), ite = msh->endBndFace(name); it != ite; ++it)
			{
				it->faceFDP.ro = ro;
				it->faceFDP.ru = u * ro;
				it->faceFDP.rv = v * ro;
				it->faceFDP.rw = w * ro;
				it->faceFDP.gamma = 1.4;
				it->faceFDP.P = P;

				it->faceFDP.rE = CellFluidDynamicsProps::calc_rE(ro, P, u, v, w, it->faceFDP.gamma);
			}

        	bndOutletNames.push_back(name);
        }
        else if( type == "RIGID" )
        {
        	bndRigidNames.push_back(name);
        }

        pBndElement = pBndElement->NextSiblingElement("boundCond");
    }
}


double** SolveEilerEq::allocate_mem()
{
	double** matrix = new double*[5];
	for(int i = 0; i < 5; i++)
	{
		matrix[i] = new double[5];
	}

	return matrix;
}

void SolveEilerEq::free_mem(double** matrix)
{
	for(int i = 0; i < 5; i++)
	{
		delete [] matrix[i];
	}

	delete [] matrix;
}

void SolveEilerEq::calc_F(double* F, CellFluidDynamicsProps cfdp)
{
	F[0] = cfdp.ru;
	F[1] = cfdp.P + cfdp.ru * cfdp.ru / cfdp.ro;
	F[2] = cfdp.ru * cfdp.rv / cfdp.ro;
	F[3] = cfdp.ru * cfdp.rw / cfdp.ro;
	F[4] = (cfdp.P + cfdp.rE) * cfdp.ru / cfdp.ro;
}

void SolveEilerEq::calc_G(double* G, CellFluidDynamicsProps cfdp)
{
	G[0] = cfdp.rv;
	G[1] = cfdp.ru * cfdp.rv / cfdp.ro;
	G[2] = cfdp.P + cfdp.rv * cfdp.rv / cfdp.ro;
	G[3] = cfdp.rv * cfdp.rw / cfdp.ro;
	G[4] = (cfdp.P + cfdp.rE) * cfdp.rv / cfdp.ro;
}

void SolveEilerEq::calc_H(double* H, CellFluidDynamicsProps cfdp)
{
	H[0] = cfdp.rw;
	H[1] = cfdp.ru * cfdp.rw / cfdp.ro;
	H[2] = cfdp.rv * cfdp.rw / cfdp.ro;
	H[3] = cfdp.P + cfdp.rw * cfdp.rw / cfdp.ro;
	H[4] = (cfdp.P + cfdp.rE) * cfdp.rw / cfdp.ro;
}

void SolveEilerEq::flux_Lax_Friedrichs(double* Flux, double* F1, CellFluidDynamicsProps cfdp1, double* F2, CellFluidDynamicsProps cfdp2, Point n)
{
	double eigen_val1 = sqrt(cfdp1.gamma * cfdp1.P / cfdp1.ro) + abs((cfdp1.ru * n.x + cfdp1.rv * n.y + cfdp1.rw * n.z)) / cfdp1.ro;
	double eigen_val2 = sqrt(cfdp2.gamma * cfdp2.P / cfdp2.ro) + abs((cfdp2.ru * n.x + cfdp2.rv * n.y + cfdp2.rw * n.z)) / cfdp2.ro;
	double alpha = max(eigen_val1, eigen_val2);

	Flux[0] = 0.5 * ( F1[0] + F2[0] - alpha * (cfdp2.ro - cfdp1.ro) );
	Flux[1] = 0.5 * ( F1[1] + F2[1] - alpha * (cfdp2.ru - cfdp1.ru) );
	Flux[2] = 0.5 * ( F1[2] + F2[2] - alpha * (cfdp2.rv - cfdp1.rv) );
	Flux[3] = 0.5 * ( F1[3] + F2[3] - alpha * (cfdp2.rw - cfdp1.rw) );
	Flux[4] = 0.5 * ( F1[4] + F2[4] - alpha * (cfdp2.rE - cfdp1.rE) );
}

void SolveEilerEq::eigen_values(double* eigen_val, double u, double v, double w, double c, Point n)
{
	double vel_n = u*n.x + v*n.y + w*n.z;

	eigen_val[0] = vel_n - c;
	eigen_val[1] = vel_n + c;
	eigen_val[2] = vel_n;
	eigen_val[3] = vel_n;
	eigen_val[4] = vel_n;
}


void SolveEilerEq::left_eigen_vecs(double** left_eigen_vecs, double u, double v, double w, double c, double gamma, Point n)
{
	double g_1 = gamma - 1;
	double g_q_2 = g_1 * 0.5 * (u*u + v*v + w*w);
	double vel_n = u*n.x + v*n.y + w*n.z;

	left_eigen_vecs[0][0] = 0.5 * (g_q_2 + c*vel_n);				left_eigen_vecs[0][1] = -0.5 * (g_1*u + c*n.x);		left_eigen_vecs[0][2] = -0.5 * (g_1*v + c*n.y);		left_eigen_vecs[0][3] = -0.5 * (g_1*w + c*n.z);		left_eigen_vecs[0][4] = 0.5 * g_1;
	left_eigen_vecs[1][0] = 0.5 * (g_q_2 - c*vel_n);				left_eigen_vecs[1][1] = -0.5 * (g_1*u - c*n.x);		left_eigen_vecs[1][2] = -0.5 * (g_1*v - c*n.y);		left_eigen_vecs[1][3] = -0.5 * (g_1*w - c*n.z);		left_eigen_vecs[1][4] = 0.5 * g_1;
	left_eigen_vecs[2][0] = (c*c - g_q_2)*n.x + c*(w*n.y - v*n.z);	left_eigen_vecs[2][1] = g_1*u*n.x;					left_eigen_vecs[2][2] = g_1*v*n.x + c*n.z;			left_eigen_vecs[2][3] = g_1*w*n.x - c*n.y;			left_eigen_vecs[2][4] = -g_1 * n.x;
	left_eigen_vecs[3][0] = (c*c - g_q_2)*n.y + c*(u*n.z - w*n.x);	left_eigen_vecs[3][1] = g_1*u*n.y - c*n.z;			left_eigen_vecs[3][2] = g_1*v*n.y;					left_eigen_vecs[3][3] = g_1*w*n.y + c*n.x;			left_eigen_vecs[3][4] = -g_1 * n.y;
	left_eigen_vecs[4][0] = (c*c - g_q_2)*n.z + c*(v*n.x - u*n.y);	left_eigen_vecs[4][1] = g_1*u*n.z + c*n.y;			left_eigen_vecs[4][2] = g_1*v*n.z - c*n.x;			left_eigen_vecs[4][3] = g_1*w*n.z;					left_eigen_vecs[4][4] = -g_1 * n.z;
}

void SolveEilerEq::right_eigen_vecs(double** right_eigen_vecs, double u, double v, double w, double c, double H, Point n)
{
	double q_2 = u*u + v*v + w*w;
	double vel_n = u*n.x + v*n.y + w*n.z;

	right_eigen_vecs[0][0] = 1;				right_eigen_vecs[0][1] = 1;				right_eigen_vecs[0][2] = n.x;								right_eigen_vecs[0][3] = n.y;								right_eigen_vecs[0][4] = n.z;
	right_eigen_vecs[1][0] = u - c*n.x;		right_eigen_vecs[1][1] = u + c*n.x;		right_eigen_vecs[1][2] = u*n.x;								right_eigen_vecs[1][3] = u*n.y - c*n.z;						right_eigen_vecs[1][4] = u*n.z + c*n.y;
	right_eigen_vecs[2][0] = v - c*n.y;		right_eigen_vecs[2][1] = v + c*n.y;		right_eigen_vecs[2][2] = v*n.x + c*n.z;						right_eigen_vecs[2][3] = v*n.y;								right_eigen_vecs[2][4] = v*n.z - c*n.x;
	right_eigen_vecs[3][0] = w - c*n.z;		right_eigen_vecs[3][1] = w + c*n.z;		right_eigen_vecs[3][2] = w*n.x - c*n.y;						right_eigen_vecs[3][3] = w*n.y + c*n.x;						right_eigen_vecs[3][4] = w*n.z;
	right_eigen_vecs[4][0] = H - c*vel_n;	right_eigen_vecs[4][1] = H + c*vel_n;	right_eigen_vecs[4][2] = 0.5*q_2*n.x + c*(v*n.z - w*n.y);	right_eigen_vecs[4][3] = 0.5*q_2*n.y + c*(w*n.x - u*n.z);	right_eigen_vecs[4][4] = 0.5*q_2*n.z + c*(u*n.y - v*n.x);
}

void SolveEilerEq::matrix_A(double** A, double** right_eigen_vecs, double* eigen_val_mat, double** left_eigen_vecs, int sign)
{
	double** temp_mat = allocate_mem();

	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			temp_mat[i][j] = left_eigen_vecs[i][j];
		}
	}

	for(int i = 0; i < 5; i++)
	{
		double temp_val = eigen_val_mat[i];

		if( (sign == 0 && eigen_val_mat[i] < 0) || (sign == 1 && eigen_val_mat[i] >= 0) )
		{
			temp_val = 0;
		}

		for(int j = 0; j < 5; j++)
			temp_mat[i][j] *= temp_val;
	}

	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			A[i][j] = 0;
			for(int k = 0; k < 5; k++)
			{
				A[i][j] += right_eigen_vecs[i][k] * temp_mat[k][j];
			}
		}
	}

	free_mem(temp_mat);
}


void SolveEilerEq::calcEilerEquation(double t_max)
{
	double tau = 0.00001;
	double t = 0;
	double eps = 1E-7;
	int max_iter = 25;

	double F_flux[5];
	double F1[5];
	double F2[5];

	double G_flux[5];
	double G1[5];
	double G2[5];

	double H_flux[5];
	double H1[5];
	double H2[5];

	MatrixSolver* matrix = new MatrixSolver( 5 * msh->cells.size() );

	double* solution;
	double* eigen_vals = new double[5];
	double** left_vecs = allocate_mem();
	double** right_vecs = allocate_mem();
	double** A_plus = allocate_mem();
	double** A_minus = allocate_mem();

	CellFluidDynamicsProps cfdp1;
	CellFluidDynamicsProps cfdp2;

	double temp_ro;
	double temp_u;
	double temp_v;
	double temp_w;
	double temp_gamma;
	double temp_H;
	double temp_P;
	double temp_rE;
	double temp_c;

	int ic, oc, index_cell1, index_cell2;


	while(t < t_max)
	{
		t += tau;

		matrix->clear();

		for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndInletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndInletNames); it != ite; ++it)
		{
			index_cell1 = it->c[0]->index;

			calc_F(F1, it->c[0]->cellFDP);
			calc_G(G1, it->c[0]->cellFDP);
			calc_H(H1, it->c[0]->cellFDP);

			calc_F(F2, it->faceFDP);
			calc_G(G2, it->faceFDP);
			calc_H(H2, it->faceFDP);

			flux_Lax_Friedrichs(F_flux, F1, it->c[0]->cellFDP, F2, it->faceFDP, it->n);
			flux_Lax_Friedrichs(G_flux, G1, it->c[0]->cellFDP, G2, it->faceFDP, it->n);
			flux_Lax_Friedrichs(H_flux, H1, it->c[0]->cellFDP, H2, it->faceFDP, it->n);

			for(int i = 0; i < 5; i++)
				matrix->subtract_f(5 * index_cell1 + i, (F_flux[i] * it->n.x + G_flux[i] * it->n.y + H_flux[i] * it->n.z) * it->S);

			temp_ro = CellFluidDynamicsProps::ro_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro);
			temp_u = CellFluidDynamicsProps::vel_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro, it->c[0]->cellFDP.ru / it->c[0]->cellFDP.ro, it->faceFDP.ru / it->faceFDP.ro);
			temp_v = CellFluidDynamicsProps::vel_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro, it->c[0]->cellFDP.rv / it->c[0]->cellFDP.ro, it->faceFDP.rv / it->faceFDP.ro);
			temp_w = CellFluidDynamicsProps::vel_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro, it->c[0]->cellFDP.rw / it->c[0]->cellFDP.ro, it->faceFDP.rw / it->faceFDP.ro);
			temp_gamma = 1.4;
			temp_H = CellFluidDynamicsProps::H_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro, (it->c[0]->cellFDP.rE + it->c[0]->cellFDP.P) / it->c[0]->cellFDP.ro, (it->faceFDP.rE + it->faceFDP.P) / it->faceFDP.ro);
			temp_rE = CellFluidDynamicsProps::rE_Rou_Avg(temp_ro, temp_u, temp_v, temp_w, temp_gamma, temp_H);
			temp_P = CellFluidDynamicsProps::P_Rou_Avg(temp_ro, temp_u, temp_v, temp_w, temp_gamma, temp_H);
			temp_c = sqrt(temp_gamma * temp_P / temp_ro);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_gamma, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, SolveEilerEq::PLUS_JACOBIAN);
			//matrix_A(A_minus, left_vecs, eigen_val_mat, right_vecs, 1);

			for(int i = 0; i < 5; i++)
				for(int j = 0; j < 5; j++)
				{
					matrix->add_matrix(5 * index_cell1 + i, 5 * index_cell1 + j, A_plus[i][j] * it->S);
				}
		}


		for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndOutletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndOutletNames); it != ite; ++it)
		{
			index_cell1 = it->c[0]->index;

			calc_F(F1, it->c[0]->cellFDP);
			calc_G(G1, it->c[0]->cellFDP);
			calc_H(H1, it->c[0]->cellFDP);

			calc_F(F2, it->faceFDP);
			calc_G(G2, it->faceFDP);
			calc_H(H2, it->faceFDP);

			flux_Lax_Friedrichs(F_flux, F1, it->c[0]->cellFDP, F2, it->faceFDP, it->n);
			flux_Lax_Friedrichs(G_flux, G1, it->c[0]->cellFDP, G2, it->faceFDP, it->n);
			flux_Lax_Friedrichs(H_flux, H1, it->c[0]->cellFDP, H2, it->faceFDP, it->n);

			for(int i = 0; i < 5; i++)
				matrix->subtract_f(5 * index_cell1 + i, (F_flux[i] * it->n.x + G_flux[i] * it->n.y + H_flux[i] * it->n.z) * it->S);

			temp_ro = CellFluidDynamicsProps::ro_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro);
			temp_u = CellFluidDynamicsProps::vel_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro, it->c[0]->cellFDP.ru / it->c[0]->cellFDP.ro, it->faceFDP.ru / it->faceFDP.ro);
			temp_v = CellFluidDynamicsProps::vel_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro, it->c[0]->cellFDP.rv / it->c[0]->cellFDP.ro, it->faceFDP.rv / it->faceFDP.ro);
			temp_w = CellFluidDynamicsProps::vel_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro, it->c[0]->cellFDP.rw / it->c[0]->cellFDP.ro, it->faceFDP.rw / it->faceFDP.ro);
			temp_gamma = 1.4;
			temp_H = CellFluidDynamicsProps::H_Rou_Avg(it->c[0]->cellFDP.ro, it->faceFDP.ro, (it->c[0]->cellFDP.rE + it->c[0]->cellFDP.P) / it->c[0]->cellFDP.ro, (it->faceFDP.rE + it->faceFDP.P) / it->faceFDP.ro);
			temp_rE = CellFluidDynamicsProps::rE_Rou_Avg(temp_ro, temp_u, temp_v, temp_w, temp_gamma, temp_H);
			temp_P = CellFluidDynamicsProps::P_Rou_Avg(temp_ro, temp_u, temp_v, temp_w, temp_gamma, temp_H);
			temp_c = sqrt(temp_gamma * temp_P / temp_ro);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_gamma, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, SolveEilerEq::PLUS_JACOBIAN);
			//matrix_A(A_minus, left_vecs, eigen_val_mat, right_vecs, 1);

			for(int i = 0; i < 5; i++)
				for(int j = 0; j < 5; j++)
				{
					matrix->add_matrix(5 * index_cell1 + i, 5 * index_cell1 + j, A_plus[i][j] * it->S);
				}
		}

		for(Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite; ++it)
		{
			oc = it->out_cell;
			ic = it->in_cell;
			index_cell1 = it->c[oc]->index;
			index_cell2 = it->c[ic]->index;

			calc_F(F1, it->c[oc]->cellFDP);
			calc_G(G1, it->c[oc]->cellFDP);
			calc_H(H1, it->c[oc]->cellFDP);

			calc_F(F1, it->c[ic]->cellFDP);
			calc_G(G1, it->c[ic]->cellFDP);
			calc_H(H1, it->c[ic]->cellFDP);

			flux_Lax_Friedrichs(F_flux, F1, it->c[oc]->cellFDP, F2, it->c[ic]->cellFDP, it->n);
			flux_Lax_Friedrichs(G_flux, G1, it->c[oc]->cellFDP, G2, it->c[ic]->cellFDP, it->n);
			flux_Lax_Friedrichs(H_flux, H1, it->c[oc]->cellFDP, H2, it->c[ic]->cellFDP, it->n);

			for(int i = 0; i < 5; i++)
				matrix->subtract_f(5 * index_cell1 + i, (F_flux[i] * it->n.x + G_flux[i] * it->n.y + H_flux[i] * it->n.z) * it->S);

			flux_Lax_Friedrichs(F_flux, F2, it->c[ic]->cellFDP, F1, it->c[oc]->cellFDP, it->n);
			flux_Lax_Friedrichs(G_flux, G2, it->c[ic]->cellFDP, G1, it->c[oc]->cellFDP, it->n);
			flux_Lax_Friedrichs(H_flux, H2, it->c[ic]->cellFDP, H1, it->c[oc]->cellFDP, it->n);

			for(int i = 0; i < 5; i++)
				matrix->add_f(5 * index_cell2 + i, (F_flux[i] * it->n.x + G_flux[i] * it->n.y + H_flux[i] * it->n.z) * it->S);

			temp_ro = CellFluidDynamicsProps::ro_Rou_Avg(it->c[ic]->cellFDP.ro, it->c[oc]->cellFDP.ro);
			temp_u = CellFluidDynamicsProps::vel_Rou_Avg(it->c[ic]->cellFDP.ro, it->c[oc]->cellFDP.ro, it->c[ic]->cellFDP.ru / it->c[ic]->cellFDP.ro, it->c[oc]->cellFDP.ru / it->c[oc]->cellFDP.ro);
			temp_v = CellFluidDynamicsProps::vel_Rou_Avg(it->c[ic]->cellFDP.ro, it->c[oc]->cellFDP.ro, it->c[ic]->cellFDP.rv / it->c[ic]->cellFDP.ro, it->c[oc]->cellFDP.rv / it->c[oc]->cellFDP.ro);
			temp_w = CellFluidDynamicsProps::vel_Rou_Avg(it->c[ic]->cellFDP.ro, it->c[oc]->cellFDP.ro, it->c[ic]->cellFDP.rw / it->c[ic]->cellFDP.ro, it->c[oc]->cellFDP.rw / it->c[oc]->cellFDP.ro);
			temp_gamma = 1.4;
			temp_H = CellFluidDynamicsProps::H_Rou_Avg(it->c[ic]->cellFDP.ro, it->c[oc]->cellFDP.ro, (it->c[ic]->cellFDP.rE + it->c[ic]->cellFDP.P) / it->c[ic]->cellFDP.ro, (it->c[oc]->cellFDP.rE + it->c[oc]->cellFDP.P) / it->c[oc]->cellFDP.ro);
			temp_rE = CellFluidDynamicsProps::rE_Rou_Avg(temp_ro, temp_u, temp_v, temp_w, temp_gamma, temp_H);
			temp_P = CellFluidDynamicsProps::P_Rou_Avg(temp_ro, temp_u, temp_v, temp_w, temp_gamma, temp_H);
			temp_c = sqrt(temp_gamma * temp_P / temp_ro);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_gamma, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, SolveEilerEq::PLUS_JACOBIAN);
			matrix_A(A_minus, right_vecs, eigen_vals, left_vecs, SolveEilerEq::MINUS_JACOBIAN);

			for(int i = 0; i < 5; i++)
				for(int j = 0; j < 5; j++)
				{
					matrix->add_matrix(5 * index_cell1 + i, 5 * index_cell1 + j, A_plus[i][j] * it->S);
					matrix->subtract_matrix(5 * index_cell2 + i, 5 * index_cell2 + j, A_plus[i][j] * it->S);

					matrix->add_matrix(5 * index_cell1 + i, 5 * index_cell2 + j, A_minus[i][j] * it->S);
					matrix->subtract_matrix(5 * index_cell2 + i, 5 * index_cell1 + j, A_minus[i][j] * it->S);
				}
		}


		for(Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
		{
			double V_tau = it->V / tau;
			for(int i = 0; i < 5; i++)
			{
				matrix->add_matrix(5 * it->index + i, 5 * it->index + i, V_tau);
			}
		}


		/*
		for(int i = 0; i < N; i++)
					 {
						//if(i % 5 == 124) //printf("%d\n\n",i / 5);
						 for(int j = 0; j < N; j++)
						 {
							 printf("%6.3f ", Matrix[i][j]);
						 }
						 printf(" | %0.3f\n", f[i]);
					 }
					 printf("\n");
		*/
		printf("%f\n",t);

		matrix->zeidel_method(eps, max_iter);
		solution = matrix->get_solution();

		for(Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
		 {
			 it->cellFDP.ro += solution[5 * it->index + 0];
			 it->cellFDP.ru += solution[5 * it->index + 1];
			 it->cellFDP.rv += solution[5 * it->index + 2];
			 it->cellFDP.rw += solution[5 * it->index + 3];
			 it->cellFDP.rE += solution[5 * it->index + 4];

			 it->cellFDP.P = CellFluidDynamicsProps::calc_P(it->cellFDP.ro, it->cellFDP.rE, it->cellFDP.ru, it->cellFDP.rv, it->cellFDP.rw, it->cellFDP.gamma);
		 }

		//if(t >= 1.0*tau) return;
	}


	delete matrix;
	delete [] eigen_vals;
	free_mem(left_vecs);
	free_mem(right_vecs);
	free_mem(A_plus);
	free_mem(A_minus);
}


void SolveEilerEq::save(const char *filename, const char *header)
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

    fprintf(out, "CELL_DATA %d\nSCALARS P double 1\nLOOKUP_TABLE default\n", cellCount);
           for (int i = 0; i < cellCount; i++)
           {
               fprintf(out, "%f\n", msh->cells[i]->cellFDP.P);
           }



    fprintf(out, "CELL_DATA %d\nVECTORS vectors double \n", cellCount);
		   for (int i = 0; i < cellCount; i++)
		   {
			   double ro = msh->cells[i]->cellFDP.ro;
			   fprintf(out, "%f %f %f\n", msh->cells[i]->cellFDP.ru/ro, msh->cells[i]->cellFDP.rv/ro, msh->cells[i]->cellFDP.rw/ro);
		   }

    fclose(out);
}
