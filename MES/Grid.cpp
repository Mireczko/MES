#include "Grid.h"
#include "GlobalData.h"

extern GlobalData *data;
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

	//DYNAMICZNA ALOKACJA TABLIC 2D
	globalMatrixH = new double *[nH*nL];
	globalMatrixC = new double *[nH*nL];
	globalMatrixHBC = new double *[nH*nL];
	for (int i = 0; i < nH*nL; i++)
	{
		globalMatrixH[i] = new double[nH*nL];
		globalMatrixC[i] = new double[nH*nL];
		globalMatrixHBC[i] = new double[nH*nL];
	}	
	for (int i = 0; i < nH*nL; i++)
	{
		for (int j = 0; j < nH*nL; j++)
		{
			globalMatrixH[i][j] = 0;
			globalMatrixC[i][j] = 0;
			globalMatrixHBC[i][j] = 0;
		}
	}

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

	//WARUNKI BRZEGOWE
	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			//WARUNEK BRZEGOWY Z DO£U
			if (elements[i].nodes[j].y == 0)
			{
				elements[i].isOnEdge[0] = true;
			}

			//WARUNEK BRZEGOWY Z PRAWEJ STRONY
			if (elements[i].getNodes()[j].x == L)
			{
				elements[i].isOnEdge[1] = true;
			}

			//WARUNEK BRZEGOWY Z GÓRY
			if (elements[i].getNodes()[j].y == H)
			{
				elements[i].isOnEdge[2] = true;
			}

			//WARUNEK BRZEGOWY Z LEWEJ STRONY
			if (elements[i].getNodes()[j].x == 0)
			{
				elements[i].isOnEdge[3] = true;
			}
		}

	}

	//OBLICZANIE MACIERZY JACOBIEGO, H, C, BC DLA KA¯DEGO ELEMENTU

	for (int i = 0; i < nE; i++)
	{
		elements[i].jacobian = new Jacobian(elements[i].getNodes());
		elements[i].matrixh = new MatrixH(elements[i].jacobian);
		elements[i].matrixc = new MatrixC(elements[i].jacobian);
		elements[i].matrixhbc = new MatrixHBC(elements[i].jacobian, elements[i].isOnEdge);
	}


	//GLOBALNA MACIERZ H
	//GLOBALNA MACIERZ C
	//GLOBALNA MACIERZ BC
	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				globalMatrixH[elements[i].nodes[j].numer][elements[i].nodes[k].numer] += elements[i].matrixh->H[j][k];
				globalMatrixC[elements[i].nodes[j].numer][elements[i].nodes[k].numer] += elements[i].matrixc->C[j][k];
				globalMatrixHBC[elements[i].nodes[j].numer][elements[i].nodes[k].numer] += elements[i].matrixhbc->hbc[j][k];
			}
		}
	}

	//C/simulation_step_time
	//H + BC + C
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			globalMatrixC[i][j] /= data->simulation_step_time;
			globalMatrixH[i][j] += globalMatrixC[i][j] + globalMatrixHBC[i][j];
		}
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

		std::cout << std::endl;
		std::cout << "Warunki: ";
		for (int j = 0; j < 4; j++)
		{
			std::cout << this->getElements()[i].isOnEdge[j] << "  ";
		}
		std::cout << std::endl <<std::endl;
	}
	//MACIERZ H/C GLOBALNA KONTROLNIE
	std::cout<<std::endl << "Matrix H + C" << std::endl<<std::endl;
	for (int i = 0; i < nL*nH; i++)
	{
		for (int j = 0; j < nL*nH; j++)
		{
			std::cout << globalMatrixH[i][j] << "  ";
		}
		std::cout << std::endl;
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
