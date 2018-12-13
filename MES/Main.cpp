#include <Eigen/Dense>
#include "Grid.h"
#include "Jacobian.h"

using namespace Eigen;


int main()
{
	int nH = 3;
	int nL = 3;
	double H = 0.025;
	double L = 0.025;

	Grid* siatka = new Grid(nH, nL, H, L);
	siatka->printGrid();

	system("pause");
}