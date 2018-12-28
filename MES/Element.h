#pragma once
#include "MatrixH.h"
#include "MatrixC.h"
#include "MatrixHBC.h"

class Element
{
	public:
		Node* nodes;
		double conductivity;
		Jacobian *jacobian;
		MatrixH* matrixh;
		MatrixC* matrixc;
		MatrixHBC* matrixhbc;
		bool *isOnEdge;
		//MatrixH* matrixh = new MatrixH(jacobian, conductivity);

	public:
		Element();
		Element(Node*, double);
		~Element();
		Node* getNodes();
		void printMatrixH();
		void printMatrixC();
		void printMatrixHBC();
};