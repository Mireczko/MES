#include "GlobalData.h"
#include <fstream>
#include <string>
#include <iostream>

GlobalData::GlobalData()
{
	std::fstream data;
	std::string filename = "data.txt";
	if (!filename.empty())
	{
		data.open(filename, std::ios::in);
		if (data.is_open())
		{
			std::string bufor;
			data >> bufor;
			initial_temperature = atof((bufor.replace(0, 20, "")).c_str());
			data >> bufor;
			simulation_time = atof((bufor.replace(0, 16, "")).c_str());
			data >> bufor;
			simulation_step_time = atof((bufor.replace(0, 21, "")).c_str());
			data >> bufor;
			ambient_temperature = atof((bufor.replace(0, 20, "")).c_str());
			data >> bufor;
			alfa = atof((bufor.replace(0, 5, "")).c_str());
			data >> bufor;
			H = atof((bufor.replace(0, 2, "")).c_str());
			data >> bufor;
			L = atof((bufor.replace(0, 2, "")).c_str());
			data >> bufor;
			nH = atoi((bufor.replace(0, 3, "")).c_str());
			data >> bufor;
			nL = atoi((bufor.replace(0, 3, "")).c_str());
			data >> bufor;
			specific_heat = atof((bufor.replace(0, 14, "")).c_str());
			data >> bufor;
			conductivity = atof((bufor.replace(0, 13, "")).c_str());
			data >> bufor;
			density = atof((bufor.replace(0, 8, "")).c_str());
			data.close();
		}
		else
		{
			std::cout << "Plik nie istnieje, lub nie ma dostêpu do pliku";
		}
	}
	else
	{
		std::cout << "Niepoprawna nazwa, lub brak pliku";
	}
}
GlobalData::~GlobalData(){}
GlobalData::GlobalData(int nH, int nL, double H, double L, double conductivity)
{
	this->nH = nH;
	this->nL = nL;
	this->H = H;
	this->L = L;
	this->conductivity = conductivity;
}

GlobalData* data = new GlobalData();