#include <Eigen/Dense>
#include "Grid.h"

using namespace Eigen;


int main()
{
	int nH = 4;
	int nL = 4;
	double H = 0.1;
	double L = 0.1;

	Grid* siatka = new Grid(nH, nL, H, L);
	siatka->printGrid();

	system("pause");
}