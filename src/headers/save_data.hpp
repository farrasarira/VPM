#ifndef SAVEDATA_H
#define SAVEDATA_H

#include <vector>
#include "particle.hpp"
#include "base_save_data.hpp"

// class base_poisson;
// class base_grid;
// class base_remeshing;
// class base_penalization;
class base_save_data;
class save_data
{
#define w10 std::setw(10) // spare width while saving data
#define w16 std::setw(16) // spare width while saving data
#define w20 std::setw(20) // spare width while saving data

#pragma region instances
	// base_poisson d_base_poisson;
	// base_grid d_base_grid;
	// base_remeshing d_base_remeshing;
	// penalization d_penalization;
	base_save_data d_base_save_data;
#pragma endregion

#pragma region internal_variables
	// internal variables
	// @param; lifetime = singleton
	// ...variables...

	// @param; lifetime = transient
	std::vector<double> _Isumx;
	std::vector<double> _Isumy;
	std::vector<double> _Isumz;
#pragma endregion

	// -- local functions
	void force_linear_impulse(int it, int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &zp,
							  const std::vector<double> &gpx, const std::vector<double> &gpy, const std::vector<double> &gpz, const std::vector<double> &sp);

	void grid_data(int np, std::vector<double> &xp, std::vector<double> &yp, std::vector<double> &sp, std::vector<double> &gpz,
				   std::vector<double> &up, std::vector<double> &vp, int &nvertn, std::vector<std::vector<double>> &vertex);
				   
	void energy_calc(int it, int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &zp,
                              const std::vector<double> &up, const std::vector<double> &vp, const std::vector<double> &wp, 
                              const std::vector<double> &gpx, const std::vector<double> &gpy, const std::vector<double> &gpz,
                              const std::vector<double> &sp);
public:
	// constructor
	save_data();
	// destructor
	~save_data();

	void particle_data_reading(int it, int &np, std::vector<double> &xp, std::vector<double> &yp, std::vector<double> &sp, std::vector<double> &gpz, std::vector<double> &up, std::vector<double> &vp);
	void particle_data_reading(int &nvertn, std::vector<std::vector<double>> &vertex);

	void output(int it, Particle &p, Body &b, double &cum_time);
};

#endif
