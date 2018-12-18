#include "MatrixC.h"
#include "GlobalData.h"

extern GlobalData *data;
MatrixC::MatrixC() {}
MatrixC::~MatrixC() {}

MatrixC::MatrixC(Jacobian *jacobian)
{
	this->jacobian = jacobian;


	//DYNAMICZNA ALOKACJA TABLICY 3D
	integral = new double **[4];

	for (int i = 0; i < 4; i++)
	{
		integral[i] = new double*[4];

		for (int j = 0; j < 4; j++)
		{
			integral[i][j] = new double[4];
		}
	}

	//DYNAMICZNA ALOKACJA TABLICY 2D
	C = new double *[4];
	for (int i = 0; i < 4; i++)
	{
		C[i] = new double[4];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			C[i][j] = 0;
		}
	}

	//c * ro * {N}{N}^T
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				integral[i][k][j] = (data->specific_heat) *(data->density) * jacobian->N[j][i] * jacobian->N[i][k] * jacobian->detJ[i];
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				C[i][j] += integral[k][i][j];
			}
		}
	}
}