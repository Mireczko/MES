#include "Grid.h"
#include "GlobalData.h"
extern GlobalData* data;

int main()
{
	Grid* siatka = new Grid(data->nH, data->nL, data->H, data->L);
	siatka->printGrid();

	system("pause");
}