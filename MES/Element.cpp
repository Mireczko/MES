#include "Element.h"

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