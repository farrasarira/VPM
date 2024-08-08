#ifndef BASE_SAVE_DATA_H
#define BASE_SAVE_DATA_H

#include <iomanip> // std::setw, std::setfill 
#include <sstream> // std::stringstream
#include <fstream> // std::c_str()
#include <vector>

class base_save_data
{
	#define w20 std::setw(20) // spare width while saving data

public:
	void force_pen(int it, int nPi, double lambda, std::vector<double>& kai, std::vector<double>& uPi, std::vector<double>& vPi,
		std::vector<double>& uSi, std::vector<double>& vSi, std::vector<double>& sPi);

	void force_pen_3d(int it, int nPi, double lambda, std::vector<double>& kai, std::vector<double>& uPi, std::vector<double>& vPi,
		std::vector<double>& wPi, std::vector<double>& uSi, std::vector<double>& vSi,  std::vector<double>& wSi, std::vector<double>& sPi);

};

#endif