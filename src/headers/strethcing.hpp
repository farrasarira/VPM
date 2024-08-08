#ifndef STRETCHING_H
#define STRETCHING_H

#include "particle.hpp"

class neighbor;
class LSMPSa;
class stretching
{
	// ** creating instances
	// neighbor d_neighbor; // linked_list, direct_find
	// LSMPSa lsmpsa_x;		 // tu calculate derivative
	// LSMPSa lsmpsa_y;
	// LSMPSa lsmpsa_z;

public:

	void main_stretching(Particle &p, std::vector<std::vector<double>> &dfdt);

};

#endif