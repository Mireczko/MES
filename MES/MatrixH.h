#pragma once
#include "Jacobian.h"
class MatrixH
{
	public:
		Jacobian jacobian;
		double ***MatrixX;
		double ***MatrixY;
		double ***sumXY;
		double conductivity=30;
		double **H;

	public:
		//Jacobian getJacobian();
		//double*** getMatrixX();
		//double*** getMatrixY();
		MatrixH(Jacobian, double);
		MatrixH();
		~MatrixH();
};