#pragma once
#include "Element.h"
class Jacobian {
	public:
		double** N;
		double* ksi;
		double* eta;
		double** dNdKsi;
		double** dNdEta;
		double** dNdX;
		double** dNdY;
		double* dXdKsi;
		double* dXdEta;
		double* dYdKsi;
		double* dYdEta;
		double* detJ;
		double* x;
		double* y;

	public:
		Jacobian();
		Jacobian(Element element);
		~Jacobian();
};