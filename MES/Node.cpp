#pragma once
#include "Node.h"
#include <iostream>

Node::Node() {}
Node::~Node() {}
Node::Node(double x, double y)
{
	this->x = x;
	this->y = y;
}

double Node::getx()
{
	return this->x; 
}
double Node::gety()
{
	return this->y;
}

void Node::setx(double x)
{
	this->x = x;
}
void Node::sety(double y)
{
	this->y = y;
}

int Node::getNumer()
{
	return numer;
}