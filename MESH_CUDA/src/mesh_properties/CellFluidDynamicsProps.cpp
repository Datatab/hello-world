#include "CellFluidDynamicsProps.h"

CellFluidDynamicsProps::CellFluidDynamicsProps(double ro, double u, double v, double w, double E,
												double P, double gamma) {
	this->ro = ro;
	this->ru = ro * u;
	this->rv = ro * v;
	this->rw = ro * w;
	this->rE = ro * E;
	this->P = P;
	this->gamma = gamma;
}

CellFluidDynamicsProps::CellFluidDynamicsProps() {
}


CellFluidDynamicsProps::~CellFluidDynamicsProps() {
}

double CellFluidDynamicsProps::calc_rE(double ro, double P, double u, double v, double w, double gamma)
{
	return P / (gamma - 1) + ro * 0.5 *(u*u + v*v + w*w);
}

double CellFluidDynamicsProps::calc_P(double ro, double rE, double ru, double rv, double rw, double gamma)
{
	return (rE - 0.5 *(ru*ru + rv*rv + rw*rw) / ro) * (gamma - 1);
}

double CellFluidDynamicsProps::ro_Rou_Avg(double ro1, double ro2)
{
	return sqrt(ro1*ro2);
}

double CellFluidDynamicsProps::vel_Rou_Avg(double ro1, double ro2, double u1, double u2)
{
	double ro1_sqrt = sqrt(ro1);
	double ro2_sqrt = sqrt(ro2);

	return ( ro1_sqrt*u1 + ro2_sqrt*u2 ) / ( ro1_sqrt + ro2_sqrt );
}

double CellFluidDynamicsProps::H_Rou_Avg(double ro1, double ro2, double H1, double H2)
{
	double ro1_sqrt = sqrt(ro1);
	double ro2_sqrt = sqrt(ro2);

	return (ro1_sqrt*H1 + ro2_sqrt*H2) / (ro1_sqrt + ro2_sqrt);
}

double CellFluidDynamicsProps::rE_Rou_Avg(double ro, double u, double v, double w, double gamma, double H)
{
	return ro * (H + (gamma - 1)*(u*u + v*v + w*w)/2) / gamma;
}

double CellFluidDynamicsProps::P_Rou_Avg(double ro, double u, double v, double w, double gamma, double H)
{
	return ro * (1 - 1.0/gamma) * (H - (u*u + v*v + w*w)/2);
}
