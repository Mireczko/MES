#include "GlobalData.h"

GlobalData::GlobalData(){}
GlobalData::~GlobalData(){}
GlobalData::GlobalData(int nH, int nL, double H, double L, double conductivity)
{
	this->nH = nH;
	this->nL = nL;
	this->H = H;
	this->L = L;
	this->conductivity = conductivity;
}