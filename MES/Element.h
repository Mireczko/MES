#pragma once
#include "Node.h"

class Element
{
	public:
		Node* nodes;
	public:
		Element();
		Element(Node*);
		~Element();
		Node* getNodes();
};