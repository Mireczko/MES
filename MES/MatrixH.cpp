#include "MatrixH.h"

MatrixH::MatrixH(){}
MatrixH::~MatrixH(){}

MatrixH::MatrixH(Jacobian jacobian, double conductivity)
{

	this->jacobian = jacobian;
	this->conductivity = conductivity;


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
				MatrixX[i][j][k] = jacobian.dNdX[i][j] * jacobian.dNdX[i][k] * jacobian.detJ[i];
				MatrixY[i][j][k] = jacobian.dNdY[i][j] * jacobian.dNdY[i][k] * jacobian.detJ[i];
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				sumXY[i][j][k] = (MatrixX[i][j][k] + MatrixY[i][j][k])*conductivity;
			}
		}
	}

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
}