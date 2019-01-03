#pragma once
#include "MatrixH.h"
#include "MatrixC.h"
#include "MatrixHBC.h"
#include "VectorP.h"

class Element
{
	public:
		Node* nodes;
		double conductivity;
		Jacobian *jacobian;
		MatrixH* matrixh;
		MatrixC* matrixc;
		MatrixHBC* matrixhbc;
		VectorP * vectorP;
		bool *isOnEdge;
		int number;

	public:
		Element();
		Element(Node*, double);
		~Element();
		Node* getNodes();
		void printMatrixH();
		void printMatrixC();
		void printMatrixHBC();
};