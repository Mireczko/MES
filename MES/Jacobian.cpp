#include "Jacobian.h"
#include "Element.h"
#include <math.h>

Jacobian::Jacobian(Node* nodes)
{
	this->eta = new double[4];
	eta[0] = -(1 / (sqrt(3)));
	eta[1] = -(1 / (sqrt(3)));
	eta[2] = (1 / (sqrt(3)));
	eta[3] = (1 / (sqrt(3)));

	this->ksi = new double[4];
	ksi[0] = -(1 / (sqrt(3)));
	ksi[1] = (1 / (sqrt(3)));
	ksi[2] = (1 / (sqrt(3)));
	ksi[3] = -(1 / (sqrt(3)));

	// DYNAMICZNA ALOKACJA PAMIÊCI DLA PÓL BED¥CYCH TABLICAMI 1D
	x = new double[4];
	y = new double[4];

	dXdKsi = new double[4];
	dXdEta = new double[4];
	dYdKsi = new double[4];
	dYdEta = new double[4];

	detJ = new double[4];

	for (int i = 0; i < 4; i++)
	{
		x[i] = nodes[i].x;
		y[i] = nodes[i].y;
	}

	// DYNAMICZNA ALOKACJA PAMIÊCI DLA PÓL BÊD¥CYCH TABLICAMI 2D
	N = new double* [4];

	dNdKsi = new double*[4];
	dNdEta = new double*[4];

	dNdX = new double*[4];
	dNdY = new double*[4];

	for (int i = 0; i < 4; i++)
	{
		N[i] = new double[4];
		dNdKsi[i] = new double[4];
		dNdEta[i] = new double[4];
		dNdX[i] = new double[4];
		dNdY[i] = new double[4];
	}


	// FUNKCJE KSZTA£TU W KA¯DYM PUNKCIE
	for (int i = 0; i < 4; i++)
	{
		N[i][0] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		N[i][1] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		N[i][2] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		N[i][3] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}

	// POCHODNE FUNKCJI KSZTA£TU WZGLÊDEM OSI KSI ORAZ ETA     dN/dKsi   ORAZ   dN/dEta
	for (int i = 0; i < 4; i++)
	{
		dNdKsi[i][0] = -0.25*(1 - eta[i]);
		dNdKsi[i][1] = 0.25*(1 - eta[i]);
		dNdKsi[i][2] = 0.25*(1 + eta[i]);
		dNdKsi[i][3] = -0.25*(1 + eta[i]);

		dNdEta[i][0] = -0.25*(1 - ksi[i]);
		dNdEta[i][1] = -0.25*(1 + ksi[i]);
		dNdEta[i][2] = 0.25*(1 + ksi[i]);
		dNdEta[i][3] = 0.25*(1 - ksi[i]);
	}

	// JAKOBIAN PRZEKSZTA£CENIA      dX/dKsi, dX/dEta, dYdKsi, dYdEta
	for (int i = 0; i < 4; i++)
	{
		dXdKsi[i] = dNdKsi[i][0] * x[0] + dNdKsi[i][1] * x[1] + dNdKsi[i][2] * x[2] + dNdKsi[i][3] * x[3];
		dXdEta[i] = dNdEta[i][0] * x[0] + dNdEta[i][1] * x[1] + dNdEta[i][2] * x[2] + dNdEta[i][3] * x[3];
		dYdKsi[i] = dNdKsi[i][0] * y[0] + dNdKsi[i][1] * y[1] + dNdKsi[i][2] * y[2] + dNdKsi[i][3] * y[3];
		dYdEta[i] = dNdEta[i][0] * y[0] + dNdEta[i][1] * y[1] + dNdEta[i][2] * y[2] + dNdEta[i][3] * y[3];
	}

	// WYZNACZNIK MACIERZY J
	for (int i = 0; i < 4; i++)
	{
		detJ[i] = (dXdKsi[i] * dYdEta[i]) - (dYdKsi[i] * dXdEta[i]);
	}

	// POCHODNE FUNKCJI KSZTA£TU PO X I Y
	for (int i = 0; i < 4; i++)
	{
		dNdX[i][0] = (1.0 / detJ[i]) * ((dYdEta[i] * dNdKsi[i][0]) - (dYdKsi[i] * dNdEta[i][0]));
		dNdX[i][1] = (1.0 / detJ[i]) * ((dYdEta[i] * dNdKsi[i][1]) - (dYdKsi[i] * dNdEta[i][1]));
		dNdX[i][2] = (1.0 / detJ[i]) * ((dYdEta[i] * dNdKsi[i][2]) - (dYdKsi[i] * dNdEta[i][2]));
		dNdX[i][3] = (1.0 / detJ[i]) * ((dYdEta[i] * dNdKsi[i][3]) - (dYdKsi[i] * dNdEta[i][3]));

		dNdY[i][0] = (1.0 / detJ[i]) * (-dXdEta[i] * dNdKsi[i][0] + dXdKsi[i] * dNdEta[i][0]);
		dNdY[i][1] = (1.0 / detJ[i]) * (-dXdEta[i] * dNdKsi[i][1] + dXdKsi[i] * dNdEta[i][1]);
		dNdY[i][2] = (1.0 / detJ[i]) * (-dXdEta[i] * dNdKsi[i][2] + dXdKsi[i] * dNdEta[i][2]);
		dNdY[i][3] = (1.0 / detJ[i]) * (-dXdEta[i] * dNdKsi[i][3] + dXdKsi[i] * dNdEta[i][3]);
	}
}

Jacobian::Jacobian() {};
Jacobian::~Jacobian() {};