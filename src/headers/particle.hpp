#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

class Body
{
public:
	int num;
	std::vector<double> x;	// body coordinate in x
	std::vector<double> y;	// body coordinate in y
	std::vector<double> z;	// body coordinate in z
	std::vector<double> uT;   // u's translational velocity
	std::vector<double> vT;   // v's translational velocity
	std::vector<double> wT;   // w's translational velocity
	std::vector<double> uR;   // u's rotational velocity
	std::vector<double> vR;   // v's rotational velocity
	std::vector<double> wR;   // w's rotational velocity
	std::vector<double> uDEF; // u's deformation velocity
	std::vector<double> vDEF; // v's deformation velocity
	std::vector<double> wDEF; // w's deformation velocity

	// public:
	// // constructor
	// Body();
	// // methods
	// // 	 getter
	// const std::vector<double>& get_x();
	// const std::vector<double>& get_y();
	// const std::vector<double>& get_uT();
	// const std::vector<double>& get_vT();
	// const std::vector<double>& get_uR();
	// const std::vector<double>& get_vR();
	// const std::vector<double>& get_uDEF();
	// const std::vector<double>& get_vDEF();
	// // 	setter
	// void set_x(const std::vector<double> &x);
	// void set_y(int size);
	// void set_uT(int size);
	// void set_vT(int size);
	// void set_uR(int size);
	// void set_vR(int size);
	// void set_uDEF(int size);
	// void set_vDEF(int size);
	// // destructor
	// ~Body();
};

class Particle
{
public:
	int num;
	std::vector<double> x;					// particle coordinates in x
	std::vector<double> y;					// particle coordiantes in y
	std::vector<double> z;					// particle coordiantes in z
	std::vector<double> s;					// core size (sigma)
	// std::vector<double> vol;				// volume
	std::vector<double> u;					// velocities in x direction
	std::vector<double> v;					// velocities in y direction
	std::vector<double> w;					// velocities in z direction
	std::vector<double> gx;					// strength in x direction
	std::vector<double> gy;					// strength in y direction
	std::vector<double> gz;					// strength in z direction
	std::vector<double> vort_x; 			// vorticity in x direction
	std::vector<double> vort_y; 			// vorticity in y direction
	std::vector<double> vort_z; 			// vorticity in z direction
	std::vector<std::vector<int>> neighbor; // neighbor index

	std::vector<bool> isActive;		// indicates active particles [@param: transient particle]
	std::vector<int> isActiveIndex; // indicates active particles' index [@param: transient particle]

	// public:
	// 	Particle();
	// 	~Particle();
};

class Cell
{
public:
	int num;
	std::vector<int> nleaf;
	std::vector<int> nchild;
	std::vector<int> parent;
	std::vector<double> xc;					// cell center coordinates in x
	std::vector<double> yc;					// cell center coordiantes in y
	std::vector<double> zc;					// cell center coordiantes in z
	std::vector<double> rc;					// radius cell
	std::vector<std::vector<double>> multipole_x;
	std::vector<std::vector<double>> multipole_y;
	std::vector<std::vector<double>> multipole_z;			
	std::vector<std::vector<int>> leaf;
	std::vector<std::vector<int>> child;

};

#endif