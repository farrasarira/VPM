#ifndef POISSON_H
#define POISSON_H

#include "particle.hpp"
#include "base_poisson_3d.hpp"

class LSMPSa;
class LSMPSb;
class neighbor;
class base_poisson;
class base_poisson_3d;
class base_remeshing;
class poisson_solver
{
	// -- creating instance
	// base_poisson d_base_poisson;
	base_poisson_3d d_base_poisson_3d;
	// neighbor d_neighbor;			// linked_list, direct_find
	// base_remeshing _base_remeshing; // inter-particle searching
	// LSMPSb _lsmpsb;

public:
	void poisson(Particle &p);
	// void poisson(Particle &p, Particle &pBase);

	// void advection(Particle &p);
	void advection(Particle &p, std::vector<double> &dfdt);
	// void advection(Particle &p, const Particle &pBase, std::vector<double> &dfdt);
};

#endif