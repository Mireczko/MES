#pragma once
#include "Jacobian.h"
class MatrixH
{
	public:
		Jacobian* jacobian;
		double ***MatrixX;
		double ***MatrixY;
		double ***sumXY;
		double **H;

	public:
		MatrixH(Jacobian*);
		MatrixH();
		~MatrixH();
};