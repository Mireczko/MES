#include "Element.h"
#include <iostream>
Element::Element() 
{
	this->nodes = new Node[4];
};

Element::Element(Node* arg, double conductivity) 
{ 
	this->nodes = arg; 
	this->conductivity = conductivity;
}

Element::~Element() {};

Node* Element::getNodes() 
{
	return this->nodes;
}

void Element::printMatrixH()
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			std::cout << this->matrixh->H[i][j] << "  ";
		}
		std::cout << std::endl;
	}
}