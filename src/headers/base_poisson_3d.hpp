#ifndef BASE_POISSON_H
#define BASE_POISSON_H

#include "particle.hpp"

class base_poisson_3d
{
public:
	// regularization function 2d
	void gaussian_function(double rij, double sij, double &gaussij);
	void regul_func_3d		(double rij, double sij, double &q);

	// -- Biot-Savart direct
	void biotsavart_direct_3d(Particle &pi, Particle &pj);
	// -- Biot-Savart FMM
	void biotsavart_treecode_3d(Particle &pi, Particle &pj);

private:
	// -- subroutine tree code 3d
	void built_tree(Particle &pi, Cell &root, const int n_crit,Cell &cells);
	void get_multipole(Particle &pi,const int p, Cell &cells, std::vector<int> &leaves,const int n_crit);
	void upward_sweep(Cell &cells);
	void eval_particle(Particle &pi, Cell &cells, const int n_crit, const double theta);

	void multipole_helper(std::vector<std::vector<double>> &multipole,double &dx,double &dy,double &dz,const int p, double &m);
	void add_child(int octant, int p, Cell &cells,const int n_crit);
	void split_cell(Particle &pi, int p, Cell &cells,const int n_crit);
	void M2M (int p, int c, Cell &cells);
	void M2M_helper(std::vector<std::vector<double>> &multipole,double &dx,double &dy,double &dz,int &p, int &c);
	void direct_sum(Particle &pi, Particle &pj);
	void evaluate(Particle &pi, int p, int idx, Cell &cells, const int n_crit, const double theta); 
};

#endif