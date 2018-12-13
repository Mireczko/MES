#include <Eigen/Dense>
#include "Grid.h"

using namespace Eigen;


int main()
{
	int nH = 2;
	int nL = 2;
	double H = 0.025;
	double L = 0.025;

	Grid* siatka = new Grid(nH, nL, H, L);
	siatka->printGrid();

	system("pause");
}