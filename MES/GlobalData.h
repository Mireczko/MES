#pragma once
class GlobalData
{
	public:
		double initial_temperature;
		double simulation_time;
		double simulation_step_time;
		double ambient_temperature;
		double alfa;
		double H;
		double L;
		int nH;
		int nL;
		double specific_heat;
		double conductivity;
		double density;
	public:
		GlobalData();
		~GlobalData();
		GlobalData(int, int, double, double, double);
};