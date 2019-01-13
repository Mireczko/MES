#pragma once
#include "Jacobian.h"
#include "MatrixC.h"
#include "MatrixHBC.h"

class VectorP
{
	public:
		double **P;
		double ***P_pc;
		double ***P_pc2;
		double ***sum;
		double ***powN;

	public:
		VectorP(Jacobian *, MatrixC*, bool*, MatrixHBC*);
		~VectorP();
};