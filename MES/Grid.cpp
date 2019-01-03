#include "Grid.h"
#include "GlobalData.h"
#include <cmath>

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

	double pulpa = this->L;
	double **globalVectorP2d;
	//DYNAMICZNA ALOKACJA TABLIC 2D i 1D
	globalMatrixH = new double *[nH*nL];
	globalMatrixC = new double *[nH*nL];
	globalMatrixHBC = new double *[nH*nL];
	globalVectorP2d = new double *[nH*nL];
	globalVectorP = new double[nH*nL];
	temperatures = new double[nH*nL];
	for (int i = 0; i < nN; i++)
	{
		globalVectorP[i] = 0;
	}

	for (int i = 0; i < nH*nL; i++)
	{
		globalMatrixH[i] = new double[nH*nL];
		globalMatrixC[i] = new double[nH*nL];
		globalMatrixHBC[i] = new double[nH*nL];
		globalVectorP2d[i] = new double[nH*nL];
	}	
	for (int i = 0; i < nH*nL; i++)
	{
		globalVectorP[i] = 0;
		temperatures[i] = 0;
		for (int j = 0; j < nH*nL; j++)
		{
			globalMatrixH[i][j] = 0;
			globalMatrixC[i][j] = 0;
			globalMatrixHBC[i][j] = 0;
			globalVectorP2d[i][j] = 0;
		}
	}

	//W Y P E £ N I A N I E    T A B L I C Y   N O D Ó W
	int k = 0;
	double iIncr = (L / (nL - 1.0)*1.0);
	double jIncr = (H / (nH - 1.0)*1.0);
	for (double i = 0; i < nL*iIncr - (0.00000000000000015); i += iIncr )
	{
		for (double j = 0; j < nH*jIncr - (0.00000000000000015); j += jIncr)
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
			elements[i].number = k;
		}
	}

	//WARUNKI BRZEGOWE
	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			//WARUNEK BRZEGOWY Z DO£U
			if ((elements[i].nodes[j].y) == 0)
			{
				elements[i].isOnEdge[0] = true;
			}

			//WARUNEK BRZEGOWY Z PRAWEJ STRONY
			if ((elements[i].nodes[j].x) == L || (elements[i].nodes[j].x) > (L - 0.000000000001))
			{
				elements[i].isOnEdge[1] = true;
			}

			//WARUNEK BRZEGOWY Z GÓRY
			if ((elements[i].nodes[j].y) == H || (elements[i].nodes[j].y) > H - 0.000000000001)
			{
				elements[i].isOnEdge[2] = true;
			}

			//WARUNEK BRZEGOWY Z LEWEJ STRONY
			if ((elements[i].nodes[j].x) == 0)
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
		elements[i].vectorP = new VectorP(elements[i].jacobian, elements[i].matrixc, elements[i].isOnEdge, elements[i].matrixhbc);
	}

	//GLOBALNA MACIERZ H
	//GLOBALNA MACIERZ C
	//GLOBALNA MACIERZ BC
	//GLOBALNY WEKTOR P 2D
	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int z = 0; z < 4; z++)
			{
				globalMatrixC[elements[i].nodes[j].numer][elements[i].nodes[z].numer] += elements[i].matrixc->C[j][z];
				globalMatrixH[elements[i].nodes[j].numer][elements[i].nodes[z].numer] += elements[i].matrixh->H[j][z];
				globalMatrixHBC[elements[i].nodes[j].numer][elements[i].nodes[z].numer] += elements[i].matrixhbc->hbc[j][z];
			}
		}
	}
	//temperature table
	for (int i = 0; i < nN; i++)
	{
		temperatures[i] = data->initial_temperature;
	}

	//Aggregation of P vector
	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				globalVectorP2d[elements[i].nodes[j].numer][elements[i].nodes[k].numer] += elements[i].vectorP->P[j][k];
			}
		}
	}

	double** HC = new double*[nN];
	double**globalVectorP2dTmp = new double*[nN];
	double* vectorP = new double[nN];
	double step_time = data->simulation_step_time;

	for (int i = 0; i < nN; i++)
	{
		HC[i] = new double[nN];
		globalVectorP2dTmp[i] = new double[nN];
	}

	for (int a = 0; a < data->simulation_time / data->simulation_step_time; a++)
	{
		//C/simulation_step_time
		//H + BC + C
		for (int i = 0; i < nN; i++)
		{
			vectorP[i] = 0;
			for (int j = 0; j < nN; j++)
			{
				HC[i][j] = globalMatrixH[i][j]+(globalMatrixC[i][j]/ data->simulation_step_time) +globalMatrixHBC[i][j];
				globalVectorP2dTmp[i][j] = globalVectorP2d[i][j] + (globalMatrixC[i][j] / data->simulation_step_time) * temperatures[j];
				vectorP[i] += globalVectorP2dTmp[i][j];
			}
		}

		double m, s;
		double **tmpGlobalMatrixH = new double*[nN];
		for (int i = 0; i < nN; i++)
		{
			tmpGlobalMatrixH[i] = new double[nN + 1];
		}

		//Gauss elimination
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN+1; j++) {
				if (j == nN) {
					tmpGlobalMatrixH[i][j] = vectorP[i];
				}
				else tmpGlobalMatrixH[i][j] = HC[i][j];
			}
		}

		for (int i = 0; i < nN - 1; i++)
		{
			for (int j = i + 1; j < nN; j++)
			{
				m = -tmpGlobalMatrixH[j][i] / tmpGlobalMatrixH[i][i];
				for (int k = i + 1; k <= nN; k++)
				{
					tmpGlobalMatrixH[j][k] += m * tmpGlobalMatrixH[i][k];
				}
			}
		}

		for (int i = nN - 1; i >= 0; i--)
		{
			s = tmpGlobalMatrixH[i][nN];
			for (int j = nN - 1; j >= i + 1; j--)
			{
				s -= tmpGlobalMatrixH[i][j] * vectorP[j];
			}
			vectorP[i] = s / tmpGlobalMatrixH[i][i];
			temperatures[i]	= vectorP[i];
		}

		double minTemp = 0, maxTemp = 0;
		for (int i = 0; i < nN; i++)
		{
			if (i == 0)
			{
				minTemp = temperatures[i];
				maxTemp = temperatures[i];
			}
			if (temperatures[i] > maxTemp)
			{
				maxTemp = temperatures[i];
			}
			if (temperatures[i] < minTemp)
			{
				minTemp = temperatures[i];
			}
		}

		std::cout << "Time: " << step_time << " seconds  " << "Min: " << minTemp << ", Max: " << maxTemp << std::endl << std::endl;
		step_time += data->simulation_step_time;
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
	//std::cout<<std::endl << "Matrix H + C" << std::endl<<std::endl;
	//for (int i = 0; i < nL*nH; i++)
	//{
	//	for (int j = 0; j < nL*nH; j++)
	//	{
	//		std::cout << globalMatrixH[i][j] << "  ";
	//	}
	//	std::cout << std::endl;
	//}
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
