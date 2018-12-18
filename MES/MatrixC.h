#pragma once
#include "Jacobian.h"
class MatrixC
{
	public:
		Jacobian *jacobian;
		double **C;
		double ***integral;

	public:
		MatrixC();
		~MatrixC();
		MatrixC(Jacobian*);
};