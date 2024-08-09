#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <vector>
#include "particle.hpp"
#include "LSMPSa.hpp"

////class base_dc;
class diffusion
{
	// ** creating instances
	// base_poisson d_base_poisson;
	// neighbor d_neighbor; // linked_list, direct_find
	LSMPSa lsmpsa_x;		 // tu calculate laplacian
	LSMPSa lsmpsa_y;
	LSMPSa lsmpsa_z;
	//// base_dc d_base_dc;

	// TODO: PSE scheme
	// ! change input [pair_i, pair_j]
	std::vector<double> pse(int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &sp, const std::vector<double> &gpz, const std::vector<int> &pair_i, const std::vector<int> &pair_j);

public:
	// void main_diffusion(Particle &p);
	void main_diffusion(Particle &p, std::vector<std::vector<double>> &dfdt);
	// void main_diffusion(Particle &p, const Particle &pBase, std::vector<double> &dfdt);
};



#endif