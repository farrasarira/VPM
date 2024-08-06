#ifndef PARTICLE_H
#define PARTICKLE_H

#include <vector>

class Body
{
public:
	int num;
	std::vector<double> x;	// body coordinate in x
	std::vector<double> y;	// body coordinate in y
	std::vector<double> z; 
	std::vector<double> uT;   // u's translational velocity
	std::vector<double> vT;   // v's translational velocity
	std::vector<double> wT;
	std::vector<double> uR;   // u's rotational velocity
	std::vector<double> vR;   // v's rotational velocity
	std::vector<double> wR;
	std::vector<double> uDEF; // u's deformation velocity
	std::vector<double> vDEF; // u's deformation velocity
	std::vector<double> wDEF;
	//tambahan ical
	std::vector<double> avelocity; // kecepatan sudut dari body
};

class Particle
{
public:
	int num;
	std::vector<int> label;
	std::vector<double> x;					// particle coordinates in x
	std::vector<double> y;					// particle coordiantes in y
	std::vector<double> z;
	std::vector<double> s;					// core size
 	std::vector<double> u;					// velocities in x direction
	std::vector<double> v;					// velocities in y direction
	std::vector<double> w;
	std::vector<double> gz;					// strength
	std::vector<double> vorticity;
	std::vector<std::vector<int>> neighbor; // neighbor index
	std::vector<bool> isActive;		// indicates active particles [@param: transient particle]
	std::vector<int> isActiveIndex; // indicates active particles' index [@param: transient particle]
	std::vector<int> isboundary;
	std::vector<bool> inside; 
	std::vector<double> boundaryval;
	
};

#endif