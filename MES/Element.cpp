#include "Element.h"

Element::Element() 
{
	this->nodes = new Node[4];
};

Element::Element(Node* arg) 
{ 
	this->nodes = arg; 
}

Element::~Element() {};

Node* Element::getNodes() 
{
	return this->nodes;
}