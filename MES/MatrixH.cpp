#include "MatrixH.h"
#include "GlobalData.h"

extern GlobalData* data;
MatrixH::MatrixH(){}
MatrixH::~MatrixH(){}

MatrixH::MatrixH(Jacobian* jacobian)
{

	this->jacobian = jacobian;

	//DYNAMICZNA ALOKACJA TABLIC 3D
	MatrixX = new double **[4];
	MatrixY = new double **[4];
	sumXY = new double **[4];

	for (int i = 0; i < 4; i++)
	{
		MatrixX[i] = new double*[4];
		MatrixY[i] = new double*[4];
		sumXY[i] = new double*[4];

		for (int j = 0; j < 4; j++)
		{
			MatrixX[i][j] = new double[4];
			MatrixY[i][j] = new double[4];
			sumXY[i][j] = new double[4];
		}
	}


	//DYNAMICZNA ALOKACJA TABLIC 2D
	H = new double *[4];
	for (int i = 0; i < 4; i++)
	{
		H[i] = new double[4];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			H[i][j] = 0;
		}
	}

	// {dN/dX}{dN/dX}^T * detJ		i	  {dN/dY}{dN/dY}^T* detJ	
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				MatrixX[i][j][k] = jacobian->dNdX[i][j] * jacobian->dNdX[i][k] * jacobian->detJ[i];
				MatrixY[i][j][k] = jacobian->dNdY[i][j] * jacobian->dNdY[i][k] * jacobian->detJ[i];
			}
		}
	}


	//SUMA {dN/dx}{dN/dx}^T  +  {dN/dy}{dN/dy}^T PRZEMNO¯ONA PRZEZ PRZEWODNOŒÆ CIEPLN¥
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				sumXY[i][j][k] = (MatrixX[i][j][k] + MatrixY[i][j][k])*(data->conductivity);
			}
		}
	}

	//MACIERZ H - SUMA ODPOWIADAJ¥CYCH ELEMENTÓW DLA PUNTKÓW CA£KOWANIA
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k=0; k < 4; k++)
			{
				H[i][j] += sumXY[k][i][j];
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			delete[]MatrixX[i][j];
			delete[]MatrixY[i][j];
			delete[]sumXY[i][j];
		}
		delete[] MatrixX[i];
		delete[] MatrixY[i];
		delete[] sumXY[i];
	}
	delete[] MatrixX;
	delete[] MatrixY;
	delete[] sumXY;
}
