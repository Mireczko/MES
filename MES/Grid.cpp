#include "Grid.h"

Grid::Grid() {}
Grid::~Grid(){}
Grid::Grid(int nH, int nL, double H, double L)
{
	this->nE = (nH - 1)*(nL - 1);
	this->nH = nH;
	this->nL = nL;
	this->H = H;
	this->L = L;
	this->nN = nH * nL;
	this->nodes = new Node[nN];



	//W Y P E £ N I A N I E    T A B L I C Y   N O D Ó W
	int k = 0;
	double iIncr = L / (nL - 1);
	double jIncr = H / (nH - 1);
	for (double i = 0; i < nL*iIncr; i += iIncr)
	{
		for (double j = 0; j < nH*jIncr; j += jIncr)
		{
			nodes[k].numer = k;
			nodes[k].setx(i);
			nodes[k].sety(j);
			k++;
		}
	}

	this->elements = new Element[nE];

	//W Y P E £ N I A N I E    T A B L I C Y   E L E M E N T Ó W
	k = 0;

	for (int i = 0; i < (nL-1); i++)
	{
		for (k; k < (i+1)*(nH - 1); k++)
		{
			elements[k].getNodes()[0] = nodes[k+i];
			elements[k].getNodes()[1] = nodes[k + nH +i];
			elements[k].getNodes()[2] = nodes[k + nH + 1+i];
			elements[k].getNodes()[3] = nodes[k + 1+i];
		}
	}

	//OBLICZANIE MACIERZY JACOBIEGO ORAZ MACIERZY H DLA KA¯DEGO ELEMENTU

	for (int i = 0; i < nE; i++)
	{
		elements[i].jacobian = new Jacobian(elements[i].getNodes());
		elements[i].matrixh = new MatrixH(elements[i].jacobian, 30);
		elements[i].matrixc = new MatrixC(elements[i].jacobian, 700, 7800);
	}
}

void Grid::printGrid()
{
	for (int i = 0; i < nE; i++)
	{
		std::cout << "Element nr. " << i << ": (" << this->getElements()[i].getNodes()[0].getx() << ",   " << this->getElements()[i].getNodes()[0].gety() << "), "
			<< "(" << this->getElements()[i].getNodes()[1].getx() << ",   " << this->getElements()[i].getNodes()[1].gety() << "), "
			<< "(" << this->getElements()[i].getNodes()[2].getx() << ",   " << this->getElements()[i].getNodes()[2].gety() << "), "
			<< "(" << this->getElements()[i].getNodes()[3].getx() << ",   " << this->getElements()[i].getNodes()[3].gety() << ")"
			<< std::endl << "Numery nodow: " << this->getElements()[i].getNodes()[0].getNumer() << ", "
			<< this->getElements()[i].getNodes()[1].getNumer() << ", "
			<< this->getElements()[i].getNodes()[2].getNumer() << ", "
			<< this->getElements()[i].getNodes()[3].getNumer();

		//MACIERZ H KONTROLNIE
		std::cout <<std::endl <<std::endl;
		elements[i].printMatrixC();
		std::cout << std::endl << std::endl;
	}
}


// G E T T E R Y     /    S E T T E R Y
Element* Grid::getElements()
{
	return this->elements;
}
void Grid::setnH(int a)
{
	this->nH = a;
}
void Grid::setnL(int a)
{
	this->nL = a;
}
void Grid::setnN(int a)
{
	this->nN = a;
}
void Grid::setnE(int a)
{
	this->nE = a;
}
int Grid::getnH()
{
	return this->nH;
}
int Grid::getnL()
{
	return this->nL;
}
int Grid::getnN()
{
	return this->nN;
}
int Grid::getnE()
{
	return this->nE;
}
