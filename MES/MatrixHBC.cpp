#include "MatrixHBC.h"
#include "GlobalData.h"
#include <math.h>

extern GlobalData* data;

MatrixHBC::MatrixHBC(){}
MatrixHBC::~MatrixHBC(){}

MatrixHBC::MatrixHBC(Jacobian *jacobian, bool *isOnEdge)
{
	sideLength = new double[4];


	//DYNAMICZNA ALOKACJA TABLIC 3D
	pow1 = new double **[4];
	pow2 = new double **[4];
	pow3 = new double **[4];
	pow4 = new double **[4];
	powN = new double **[4];
	bc1 = new double **[4];
	bc2 = new double **[4];
	sum = new double **[4];

	for (int i = 0; i < 4; i++)
	{
		pow1[i] = new double*[4];
		pow2[i] = new double*[4];
		pow3[i] = new double*[4];
		pow4[i] = new double*[4];
		powN[i] = new double*[4];
		bc1[i] = new double*[4];
		bc2[i] = new double*[4];
		sum[i] = new double*[4];

		for (int j = 0; j < 4; j++)
		{
			pow1[i][j] = new double[4];
			pow2[i][j] = new double[4];
			pow3[i][j] = new double[4];
			pow4[i][j] = new double[4];
			powN[i][j] = new double[4];
			bc1[i][j] = new double[4];
			bc2[i][j] = new double[4];
			sum[i][j] = new double[4];
		}
	}

	//DYNAMICZNA ALOKACJA TABLIC 2D
	hbc = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		hbc[i] = new double[4];
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			hbc[i][j] = 0;
		}
	}

	detJ = new double[4];

	//DLUGOSCI BOKOW
	sideLength[0] = sqrt(pow(jacobian->x[0] - jacobian->x[1], 2) + pow(jacobian->y[0] - jacobian->y[1], 2));
	sideLength[1] = sqrt(pow(jacobian->x[2] - jacobian->x[1], 2) + pow(jacobian->y[2] - jacobian->y[1], 2));
	sideLength[2] = sqrt(pow(jacobian->x[3] - jacobian->x[2], 2) + pow(jacobian->y[3] - jacobian->y[2], 2));
	sideLength[3] = sqrt(pow(jacobian->x[3] - jacobian->x[0], 2) + pow(jacobian->y[3] - jacobian->y[0], 2));

	for (int i = 0; i < 4; i++)
	{
		detJ[i] = sideLength[i] / 2;
	}


	//WARTOSCI FUNKCJI KSZTALTU W PUNKTACH CALKOWANIA
	double ksi[2] = { -1 / sqrt(3), 1 / sqrt(3) };
	double eta[2] = { -1, -1 };
	for (int i = 0; i < 2; i++) {
		powN[0][i][0] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		powN[0][i][1] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		powN[0][i][2] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		powN[0][i][3] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}


	ksi[0] = 1; ksi[1] = 1;
	eta[0] = -1 / sqrt(3); eta[1] = 1 / sqrt(3);
	for (int i = 0; i < 2; i++) {
		powN[1][i][0] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		powN[1][i][1] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		powN[1][i][2] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		powN[1][i][3] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}

	ksi[0] = 1 / sqrt(3); ksi[1] = -1 / sqrt(3);
	eta[0] = 1; eta[1] = 1;

	for (int i = 0; i < 2; i++) {
		powN[2][i][0] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		powN[2][i][1] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		powN[2][i][2] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		powN[2][i][3] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}

	ksi[0] = -1; ksi[1]= -1;
	eta[0] = 1 / sqrt(3); eta[1] = -1 / sqrt(3);

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
				bc1[i][j][k] = powN[i][0][k] * powN[i][0][j] * data->alfa;
				bc2[i][j][k] = powN[i][1][k] * powN[i][1][j] * data->alfa;
				sum[i][j][k] = (bc1[i][j][k] + bc2[i][j][k]) * sideLength[i] / 2;
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				hbc[i][j] += sum[k][i][j] * isOnEdge[k];

			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			delete[]pow1[i][j];
			delete[]pow2[i][j];
			delete[]pow3[i][j];
			delete[]pow4[i][j];
			delete[]powN[i][j];
			delete[]bc1[i][j];
			delete[]bc2[i][j];
			delete[]sum[i][j];
		}
		delete[] pow1[i];
		delete[] pow2[i];
		delete[] pow3[i];
		delete[] pow4[i];
		delete[] powN[i];
		delete[] bc1[i];
		delete[] bc2[i];
		delete[] sum[i];
	}
	delete[] pow1;
	delete[] pow2;
	delete[] pow3;
	delete[] pow4;
	delete[] powN;
	delete[] bc1;
	delete[] bc2;
	delete[] sum;
}
