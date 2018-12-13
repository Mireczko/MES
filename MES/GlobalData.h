#pragma once
class GlobalData
{
	public:
		int nH;
		int nL;
		double H;
		double L;
		double conductivity;
	public:
		GlobalData();
		~GlobalData();
		GlobalData(int, int, double, double, double);
};