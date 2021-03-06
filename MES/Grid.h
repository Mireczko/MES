#pragma once
#include "Element.h"
#include <iostream>

class Grid
{ 
	public:
		int nH, nL, nN, nE;
		long double L, H;
		Node* nodes;
		Element* elements;
		double** globalMatrixH;
		double** globalMatrixC;
		double** globalMatrixHBC;
		double* globalVectorP;

	public:
		Grid();
		Grid(int, int, double, double);
		~Grid();
		Element* getElements();
		void setnH(int);
		void setnL(int);
		void setnN(int);
		void setnE(int);
		int getnH();
		int getnL();
		int getnN();
		int getnE();
		void printGrid();	
};