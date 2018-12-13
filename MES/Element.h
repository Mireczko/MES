#pragma once
#include "MatrixH.h"

class Element
{
	public:
		Node* nodes;
		double conductivity;
		Jacobian jacobian;
		MatrixH matrixh;
		//Jacobian* jacobian = new Jacobian(nodes);
		//MatrixH* matrixh = new MatrixH(jacobian, conductivity);

	public:
		Element();
		Element(Node*, double);
		~Element();
		Node* getNodes();
};