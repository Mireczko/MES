#pragma once
#include "MatrixH.h"
#include "MatrixC.h"

class Element
{
	public:
		Node* nodes;
		double conductivity;
		Jacobian *jacobian;
		MatrixH* matrixh;
		MatrixC* matrixc;
		//MatrixH* matrixh = new MatrixH(jacobian, conductivity);

	public:
		Element();
		Element(Node*, double);
		~Element();
		Node* getNodes();
		void printMatrixH();
		void printMatrixC();
};