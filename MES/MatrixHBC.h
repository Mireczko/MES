#pragma once
#include "Jacobian.h"

class MatrixHBC
{	
	public:
		double ***pow1;
		double ***pow2;
		double ***pow3;
		double ***pow4;
		double ***powN;
		double ***bc1;
		double ***bc2;
		double ***sum;
		double **hbc;

		double *sideLength;
		//Jacobian *jacobian;

	public:
		MatrixHBC(Jacobian *, bool*);
		MatrixHBC();
		~MatrixHBC();
};