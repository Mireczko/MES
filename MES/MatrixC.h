#pragma once
#include "Jacobian.h"
class MatrixC
{
	public:
		Jacobian *jacobian;
		double **C;
		double ***integral;
		double c, ro;
	
	public:
		MatrixC();
		~MatrixC();
		MatrixC(Jacobian*, double, double);
};