#pragma once
#include "MatrixH.h"

class Element
{
	public:
		Node* nodes;
		double conductivity;
		Jacobian *jacobian;
		MatrixH* matrixh;
		//MatrixH* matrixh = new MatrixH(jacobian, conductivity);

	public:
		Element();
		Element(Node*, double);
		~Element();
		Node* getNodes();
		void printMatrixH();
};