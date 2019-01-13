#include "VectorP.h"
#include "GlobalData.h"
#include <iostream>

extern GlobalData* data;

VectorP::VectorP(Jacobian *jacobian, MatrixC* matrixC, bool* isOnEdge, MatrixHBC* matrixHBC)
{
	P = new double*[4];
	P_pc = new double **[4];
	P_pc2 = new double **[4];
	sum = new double **[4];
	powN = new double **[4];


	double* Ksi = new double[2]{ -1, -1 };
	double* Eta = new double[2]{ 1, -1};


	for (int i = 0; i < 4; i++)
	{
		P[i] = new double[4];
	}

	for (int i = 0; i < 4; i++)
	{
		powN[i] = new double*[2];
		sum[i] = new double*[4];
		P_pc[i] = new double*[4];
		P_pc2[i] = new double*[4];
		for (int j = 0; j < 4; j++)
		{
			powN[i][j] = new double[4];
			sum[i][j] = new double[4];
			P_pc[i][j] = new double[4];
			P_pc2[i][j] = new double[4];
		}
	}


	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				P_pc[i][j][k] = 0;
				P_pc2[i][j][k] = 0;
			}
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			P[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				sum[i][j][k] = 0;
			}
		}
	}
	   
	//WARTOSCI FUNKCJI KSZTALTU W PUNKTACH CALKOWANIA
	double ksi[2] = { -1, 1 };
	double eta[2] = { -1, -1 };
	for (int i = 0; i < 2; i++) {
		powN[0][i][0] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		powN[0][i][1] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		powN[0][i][2] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		powN[0][i][3] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}


	ksi[0] = 1; ksi[1] = 1;
	eta[0] = -1; eta[1] = 1;
	for (int i = 0; i < 2; i++) {
		powN[1][i][0] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		powN[1][i][1] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		powN[1][i][2] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		powN[1][i][3] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}

	ksi[0] = 1; ksi[1] = -1;
	eta[0] = 1; eta[1] = 1;

	for (int i = 0; i < 2; i++) {
		powN[2][i][0] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		powN[2][i][1] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		powN[2][i][2] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		powN[2][i][3] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}

	ksi[0] = -1; ksi[1] = -1;
	eta[0] = 1; eta[1] = -1;

	for (int i = 0; i < 2; i++) {
		powN[3][i][0] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		powN[3][i][1] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		powN[3][i][2] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		powN[3][i][3] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				P_pc[i][j][k] = powN[i][0][k] * powN[i][0][j] * data->alfa*data->ambient_temperature;
				P_pc2[i][j][k] = powN[i][1][k] * powN[i][1][j] * data->alfa*data->ambient_temperature;
				sum[i][j][k] = (P_pc[i][j][k] + P_pc2[i][j][k]) * matrixHBC->detJ[j];
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				P[i][j] += sum[k][i][j] * isOnEdge[k];
			}
		}
	}

	//for (int i = 0; i <4; i++)
	//{
	//	for (int j = 0; j <4; j++)
	//	{
	//		delete[]P_pc[i][j];
	//		delete[]P_pc2[i][j];
	//		delete[]sum[i][j];
	//		delete[]powN[i][j];
	//	}
	//	delete[] P_pc[i];
	//	delete[] P_pc2[i];
	//	delete[] sum[i];
	//	//delete[] powN[i];
	//}
	//delete[] P_pc;
	//delete[] P_pc2;
	//delete[] sum;
	//delete[] powN;

};

VectorP::~VectorP() {};