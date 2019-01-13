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

	double **globalVectorP2d;
	//DYNAMICZNA ALOKACJA TABLIC 2D i 1D
	globalMatrixH = new double *[nH*nL];
	globalMatrixC = new double *[nH*nL];
	globalMatrixHBC = new double *[nH*nL];
	globalVectorP2d = new double *[nH*nL];
	globalVectorP = new double[nH*nL];
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
	double iIncr = (L / (nL - 1.0)*1.0);	//nL-1 to ilosc odcinkow po dlugosci
	double jIncr = (H / (nH - 1.0)*1.0);	//nL-H to ilocs odcinkow po wysokosci
	for (double i = 0; i < nL*iIncr - (0.00000000000000015); i += iIncr )		//petla po dlugosci
	{
		for (double j = 0; j < nH*jIncr - (0.00000000000000015); j += jIncr)	//petla po wysokosci
		{
			nodes[k].numer = k;
			nodes[k].setx(i);			//ustawianie wartosci x punktu
			nodes[k].sety(j);			//ustawianie wartosci y punktu
			k++;
		}
	}

	this->elements = new Element[nE];

	//W Y P E £ N I A N I E    T A B L I C Y   E L E M E N T Ó W
	k = 0;

	for (int i = 0; i < (nL-1); i++)			//petla po ilosci odcinkow po dlugosci
	{
		for (k; k < (i+1)*(nH - 1); k++)		//petla po ilosci odcinkow po wysokosci
		{
			elements[k].getNodes()[0] = nodes[k+i];				//lewy dolny
			elements[k].getNodes()[1] = nodes[k + nH +i];		//prawy dolny
			elements[k].getNodes()[2] = nodes[k + nH + 1+i];	//prawy górny
			elements[k].getNodes()[3] = nodes[k + 1+i];			//lewy górny
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
			for (int k = 0; k < 4; k++)
			{
				globalMatrixC[elements[i].nodes[j].numer][elements[i].nodes[k].numer] += elements[i].matrixc->C[j][k];
				globalMatrixH[elements[i].nodes[j].numer][elements[i].nodes[k].numer] += elements[i].matrixh->H[j][k];
				globalMatrixHBC[elements[i].nodes[j].numer][elements[i].nodes[k].numer] += elements[i].matrixhbc->hbc[j][k];
				globalVectorP2d[elements[i].nodes[j].numer][elements[i].nodes[k].numer] += elements[i].vectorP->P[j][k];
			}
		}
	}
	//Ustawienie temperatur wêz³owych na temperaturê pocz¹tkow¹
	for (int i = 0; i < nN; i++)
	{
		nodes[i].temperature = data->initial_temperature;
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
				globalVectorP2dTmp[i][j] = globalVectorP2d[i][j] + (globalMatrixC[i][j] / data->simulation_step_time) * nodes[j].temperature;
				vectorP[i] += globalVectorP2dTmp[i][j];
			}
		}

		double m, s;
		double **tmpGlobalMatrixH = new double*[nN];
		for (int i = 0; i < nN; i++)
		{
			tmpGlobalMatrixH[i] = new double[nN + 1];	//macierz wspó³czynnikow z dopisanym wektorem P
		}

		//Eliminacja gaussa
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN+1; j++) {
				if (j == nN) {
					tmpGlobalMatrixH[i][j] = vectorP[i];//dopisujemy wektor P
				}
				else tmpGlobalMatrixH[i][j] = HC[i][j];//macierz wspolczynnikow
			}
		}
		
		//Eliminacja wspó³czynników
		for (int i = 0; i < nN - 1; i++)
		{
			for (int j = i + 1; j < nN; j++)
			{
				m = -tmpGlobalMatrixH[j][i] / tmpGlobalMatrixH[i][i];	//m mnoznik przez który mnozone sa elementy macierzy
				for (int k = i + 1; k <= nN; k++)
				{
					tmpGlobalMatrixH[j][k] += m * tmpGlobalMatrixH[i][k];
				}
			}
		}

		//Wyliczanie niewiadomych
		for (int i = nN - 1; i >= 0; i--)
		{
			s = tmpGlobalMatrixH[i][nN];								//s zlicza sume iloczynow
			for (int j = nN - 1; j >= i + 1; j--)
			{
				s -= tmpGlobalMatrixH[i][j] * vectorP[j];
			}
			vectorP[i] = s / tmpGlobalMatrixH[i][i];
			nodes[i].temperature = vectorP[i];
		}

		double minTemp = nodes[0].temperature;
		double maxTemp = nodes[0].temperature;
		for (int i = 1; i < nN; i++)
		{
			if (nodes[i].temperature > maxTemp)
			{
				maxTemp = nodes[i].temperature;
			}
			if (nodes[i].temperature < minTemp)
			{
				minTemp = nodes[i].temperature;
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
