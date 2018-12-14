#include <Eigen/Dense>
#include "Grid.h"

using namespace Eigen;


int main()
{
	int nH = 5;
	int nL = 3;
	double H = 0.1;
	double L = 0.050;

	Grid* siatka = new Grid(nH, nL, H, L);
	siatka->printGrid();

	system("pause");
}