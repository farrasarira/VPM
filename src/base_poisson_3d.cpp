// !!!!============================================================
// !!!!==== subroutine for biotsavart-treecode 3D =================
// !!!!============================================================
#include "base_poisson_3d.hpp"
// #include "poisson.hpp"
#include <iostream>
#include <math.h>
#include "parameters.hpp"
// !!!!============================================================
// !!!!============================================================

void base_poisson_3d::built_tree(Particle &pi, Cell &root, const int n_crit, Cell &cells)
{
	// """Construct a hierarchical octree to store the particles and return 
    // the tree (list) of cells.
    
    // Arguments:
    //     particles: the list of particles.
    //     root: the root cell.
    //     n_crit: maximum number of particles in a single cell.
    
    // Returns:
    //     cells: the list of cells.
    
    // """
    // # set root cell
	cells.xc.push_back(root.xc[0]);
	cells.yc.push_back(root.yc[0]);
	cells.zc.push_back(root.zc[0]);
	cells.rc.push_back(root.rc[0]);
	cells.num = 1;
	cells.parent.push_back(0);
	cells.nchild.push_back(0);
	cells.nleaf.push_back(0);
	
	for(int i=0;i<n_crit;i++){
		cells.leaf[i].push_back(0);
	}
	for(int i=0;i<8;i++){
		cells.child[i].push_back(0);
	}
	
	// # build tree
    int n = pi.num;
	for(int i = 0; i < n; i++){
		// # traverse from the root down to a leaf cell
		int curr = 0;
		while (cells.nleaf[curr] >= n_crit)
		{
			cells.nleaf[curr] += 1;
			int octant = (pi.x[i] > cells.xc[curr]) + ((pi.y[i] > cells.yc[curr]) << 1) + ((pi.z[i] > cells.zc[curr]) << 2);
			// # if there is no child cell in the particles octant, then create one
			if(!(cells.nchild[curr] & (1 << octant))){
				add_child(octant,curr,cells,n_crit);
			}
			curr = cells.child[octant][curr];		
		}
		// # allocate the particle in the leaf cell
		cells.leaf[cells.nleaf[curr]][curr] = i;
		cells.nleaf[curr] += 1;
		//  # check whether to split or not
		if(cells.nleaf[curr]>= n_crit)
		{
			split_cell(pi,curr,cells,n_crit);
		}
	}

}

void base_poisson_3d::get_multipole(Particle &pi,const int p, Cell &cells, std::vector<int> &leaves,const int n_crit)
{
	// """Calculate multipole arrays for all leaf cells under cell p. If leaf
    // number of cell p is equal or bigger than n_crit (non-leaf), traverse down
    // recursively. Otherwise (leaf), calculate the multipole arrays for leaf cell p.
    
    // Arguments:
    //     p: current cell's index.
    //     cells: the list of cells.
    //     leaves: the array of all leaf cells.
    //     n_crit: maximum number of particles in a leaf cell.     
    // """
   
    // # if the current cell p is not a leaf cell, then recursively traverse down
	if( cells.nleaf[p] >= n_crit )
	{
		for(int c = 0; c < 8; c++)
		{
			if(cells.nchild[p] & (1 << c))
			{
				get_multipole(pi,cells.child[c][p],cells,leaves,n_crit);
			}
		}
	}
	// # otherwise cell p is a leaf cell
	else
	{
		// # loop in leaf particles, do P2M
		for(int i = 0; i < cells.nleaf[p]; i++)
		{
			int l = cells.leaf[i][p];
			double dx = cells.xc[p] - pi.x[l];
			double dy = cells.yc[p] - pi.y[l];
			double dz = cells.zc[p] - pi.z[l];

			double mx,my,mz;
			mx = pi.gx[l];
			my = pi.gy[l];
			mz = pi.gz[l];
			
			// cells.multipole[0][p] += m;
			// cells.multipole[1][p] += m * dx;
			// cells.multipole[2][p] += m * dy;
			// cells.multipole[3][p] += m * dz;
			// cells.multipole[4][p] += m * dx * dx/2;
			// cells.multipole[5][p] += m * dy * dy/2;
			// cells.multipole[6][p] += m * dz * dz/2;
			// cells.multipole[7][p] += m * dx * dy/2;
			// cells.multipole[8][p] += m * dy * dz/2;
			// cells.multipole[9][p] += m * dz * dx/2;

			multipole_helper(cells.multipole_x,dx,dy,dz,p,mx);
			multipole_helper(cells.multipole_y,dx,dy,dz,p,my);
			multipole_helper(cells.multipole_z,dx,dy,dz,p,mz);
		}
		leaves.push_back(p);
	}
}

void base_poisson_3d::multipole_helper(std::vector<std::vector<double>> &multipole,double &dx,double &dy,double &dz,const int p, double &m)
{
	multipole[0][p] += m;
	multipole[1][p] += m * dx;
	multipole[2][p] += m * dy;
	multipole[3][p] += m * dz;
	multipole[4][p] += m * dx * dx/2;
	multipole[5][p] += m * dy * dy/2;
	multipole[6][p] += m * dz * dz/2;
	multipole[7][p] += m * dx * dy;
	multipole[8][p] += m * dy * dz;
	multipole[9][p] += m * dz * dx;
}

void base_poisson_3d::upward_sweep(Cell &cells)
{
    // """Traverse from leaves to root, in order to calculate multipoles of all the cells.
    
    // Arguments:
    //     cells: the list of cells.    
    // """
	for(int c = cells.num-1; c > 0; c--)
	{
		int p = cells.parent[c];
		M2M(p, c, cells);
	}
}

void base_poisson_3d::eval_particle(Particle &pi, Cell &cells, const int n_crit, const double theta)
{
    // """Evaluate the gravitational potential at all target points 
    
    // Arguments:
    //     particles: the list of particles.
    //     cells: the list of cells.
    //     n_crit: maximum number of particles in a single cell.
    //     theta: tolerance parameter.    
    // """
	for(int i=0; i<pi.num; i++)
	{
		evaluate(pi,0,i,cells,n_crit,theta);
	}
}

void base_poisson_3d::add_child(int octant, int p, Cell &cells, const int n_crit)
{
    // """Add a cell to the end of cells list as a child of p, initialize the
    // center and radius of the child cell c, and establish mutual reference
    // between child c and parent p.
    
    // Arguments:
    //     octant: reference to one of the eight divisions in three dimensions.
    //     p: parent cell index in cells list.
    //     cells: the list of cells.
    //     n_crit: maximum number of particles in a leaf cell.
    // """
	
	// Create a new cell
	double _xc,_yc,_zc,_rc;
	int c = cells.num;
	_rc = cells.rc[p]/2;
	_xc = cells.xc[p] + _rc * ((octant & 1) * 2 - 1);
	_yc = cells.yc[p] + _rc * ((octant & 2) - 1);
	_zc = cells.zc[p] + _rc * ((octant & 4) / 2 - 1);

	cells.xc.push_back(_xc);
	cells.yc.push_back(_yc);
	cells.zc.push_back(_zc);
	cells.rc.push_back(_rc);
	cells.num += 1;
	cells.parent.push_back(p);
	cells.nchild.push_back(0);
	cells.nleaf.push_back(0);
	for(int i=0;i<n_crit;i++){
		cells.leaf[i].push_back(0);
	}
	for(int i=0;i<8;i++){
		cells.child[i].push_back(0);
	}
	// Establish mutual reference in the cells list
	cells.child[octant][p] = c;
	cells.nchild[p] = (cells.nchild[p] | (1 << octant)); 
}
	
void base_poisson_3d::split_cell(Particle &pi, int p, Cell &cells,const int n_crit)
{
	// """Loop in parent p's leafs and reallocate the particles to subcells. 
    // If a subcell has not been created in that octant, create one using add_child. 
    // If the subcell c's leaf number exceeds n_crit, split the subcell c recursively.
    
    // Arguments: 
    //     particles: the list of particles.
    //     p: parent cell index in cells list.
    //     cells: the list of cells.
    //     n_crit: maximum number of particles in a leaf cell.
    // """
	
	// # loop in the particles stored in the parent cell that you want to split
	for(int i=0; i< cells.nleaf[p];i++){
		
		int l = cells.leaf[i][p];
		int octant = (pi.x[l] > cells.xc[p]) + ((pi.y[l] > cells.yc[p]) << 1) + ((pi.z[l] > cells.zc[p]) << 2); // finds the particle's octant
		
		// # if there is not a child cell in the particles octant, then create one
		if(!(cells.nchild[p] & (1 << octant))){
			add_child(octant,p,cells,n_crit);
		}

		// # reallocate the particle in the child cell
		int c = cells.child[octant][p];
		cells.leaf[cells.nleaf[c]][c] = l;
		cells.nleaf[c]+= 1;	
		// # check if the child reach n_crit
		if( cells.nleaf[c]>=n_crit){
			split_cell(pi,c,cells,n_crit);
		}
	}

}

void base_poisson_3d::M2M (int p, int c, Cell &cells)
{
    // """Calculate parent cell p's multipole array based on child cell c's 
    // multipoles.
    
    // Arguments:
    //     p: parent cell index in cells list.
    //     c: child cell index in cells list.
    //     cells: the list of cells.
    // """
	
	double dx,dy,dz;
	dx = cells.xc[p] - cells.xc[c];
	dy = cells.yc[p] - cells.yc[c];
	dz = cells.zc[p] - cells.zc[c];

	// cells.multipole[0][p] += cells.multipole[0][c];
	// cells.multipole[1][p] += cells.multipole[1][c] + cells.multipole[0][c] * dx;
	// cells.multipole[2][p] += cells.multipole[2][c] + cells.multipole[0][c] * dy;
	// cells.multipole[3][p] += cells.multipole[3][c] + cells.multipole[0][c] * dz;
	// cells.multipole[4][p] += cells.multipole[4][c] + cells.multipole[0][c] * dx * dx/2 + cells.multipole[1][c] * dx;
	// cells.multipole[5][p] += cells.multipole[5][c] + cells.multipole[0][c] * dy * dy/2 + cells.multipole[2][c] * dy;
	// cells.multipole[6][p] += cells.multipole[6][c] + cells.multipole[0][c] * dz * dz/2 + cells.multipole[3][c] * dz;
	// cells.multipole[7][p] += cells.multipole[7][c] + cells.multipole[0][c] * dx * dy/2 + cells.multipole[2][c] * dx/2 + cells.multipole[1][c] * dy/2;
	// cells.multipole[8][p] += cells.multipole[8][c] + cells.multipole[0][c] * dy * dz/2 + cells.multipole[3][c] * dy/2 + cells.multipole[2][c] * dz/2;
	// cells.multipole[9][p] += cells.multipole[9][c] + cells.multipole[0][c] * dz * dx/2 + cells.multipole[1][c] * dz/2 + cells.multipole[3][c] * dx/2;	
	
	M2M_helper(cells.multipole_x,dx,dy,dz,p,c);
	M2M_helper(cells.multipole_y,dx,dy,dz,p,c);
	M2M_helper(cells.multipole_z,dx,dy,dz,p,c);

}

void base_poisson_3d::M2M_helper(std::vector<std::vector<double>> &multipole,double &dx,double &dy,double &dz, int &p, int &c)
{
	multipole[0][p] += multipole[0][c];
	multipole[1][p] += multipole[1][c] + multipole[0][c] * dx;
	multipole[2][p] += multipole[2][c] + multipole[0][c] * dy;
	multipole[3][p] += multipole[3][c] + multipole[0][c] * dz;
	multipole[4][p] += multipole[4][c] + multipole[0][c] * dx * dx/2 + multipole[1][c] * dx;
	multipole[5][p] += multipole[5][c] + multipole[0][c] * dy * dy/2 + multipole[2][c] * dy;
	multipole[6][p] += multipole[6][c] + multipole[0][c] * dz * dz/2 + multipole[3][c] * dz;
	multipole[7][p] += multipole[7][c] + multipole[0][c] * dx * dy + multipole[2][c] * dx + multipole[1][c] * dy;
	multipole[8][p] += multipole[8][c] + multipole[0][c] * dy * dz + multipole[3][c] * dy + multipole[2][c] * dz;
	multipole[9][p] += multipole[9][c] + multipole[0][c] * dz * dx + multipole[1][c] * dz + multipole[3][c] * dx;
}

void base_poisson_3d::direct_sum(Particle &pi, Particle &pj)
{
	// internal variables
	double dx, dy, dz;
	double rij3, rij2, rij, sij, sij2;
	double gsig;

	// #pragma omp parallel for private(j)
	for (int i = 1; i <= pi.num; i++)
	{
		for (int j = 1; j <= pj.num; j++)
		{
			dx = pi.x[i - 1] - pj.x[j - 1];
			dy = pi.y[i - 1] - pj.y[j - 1];
			dz = pi.z[i - 1] - pj.z[j - 1];
			rij2 = std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2);
			rij = std::sqrt(rij2);
			rij3 = rij2 * rij;

			if (rij2 > 0.0e0)
			{
				sij2 = (std::pow(pi.s[i - 1], 2) + std::pow(pj.s[j - 1], 2)) * 0.25e0;
				sij = std::sqrt(sij2);
				regul_func_3d(rij, sij, gsig);

				pi.u[i - 1] = pi.u[i - 1] - (dy * pj.gz[j - 1] - dz * pj.gy[j - 1]) * gsig / rij3;
				pi.v[i - 1] = pi.v[i - 1] - (dz * pj.gx[j - 1] - dx * pj.gz[j - 1]) * gsig / rij3;
				pi.w[i - 1] = pi.w[i - 1] - (dx * pj.gy[j - 1] - dy * pj.gx[j - 1]) * gsig / rij3;
			}
		}
	}
}

void base_poisson_3d::evaluate(Particle &pi, int p, int idx, Cell &cells, const int n_crit, const double theta)
{
	//  """Evaluate the gravitational potential at a target point i, 
    // caused by source particles cell p. If nleaf of cell p is less 
    // than n_crit (leaf), use direct summation. Otherwise (non-leaf), loop
    // in p's child cells. If child cell c is in far-field of target particle i,
    // use multipole expansion. Otherwise (near-field), call the function
    // recursively.
    
    // Arguments:
    //     particles: the list of particles
    //     p: cell index in cells list
    //     i: target particle index
    //     cells: the list of cells
    //     n_crit: maximum number of leaves in a single cell
    //     theta: tolerance parameter    
    // """
   
    // # non-leaf cell
	if(cells.nleaf[p] >= n_crit)
	{
		// # loop in p's child cells (8 octants)
		for(int octant = 0; octant < 8; octant++)
		{
			if(cells.nchild[p] & (1 << octant))
			{
				int c = cells.child[octant][p];
				double dx,dy,dz,r,r2;
				dx = pi.x[idx] - cells.xc[c];
                dy = pi.y[idx] - cells.yc[c];
				dz = pi.z[idx] - cells.zc[c];
				r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
				r = std::sqrt(r2);

				// # near-field child cell
				if( cells.rc[c] > theta*r)
				{
					evaluate(pi,c,idx,cells,n_crit,theta);
				}
				// # far-field child cell
				else
				{
					double r3,r5,r7;
					r3 = r2*r;
					r5 = r3*r2;
					r7 = r5*r2;

					// # calculate the weight for each multipole
					std::vector<double> weight_x(10);
					std::vector<double> weight_y(10);
					std::vector<double> weight_z(10);

					// f = dx/r^3
					weight_x[0] = dx/r3;
					weight_x[1] = -(3*dx*dx/r5) + (1/r3);
					weight_x[2] = -(3*dx*dy/r5);
					weight_x[3] = -(3*dx*dz/r5);
					weight_x[4] = (15*dx*dx*dx/r7) - (9*dx/r5);
					weight_x[5] = (15*dx*dy*dy/r7) - (3*dx/r5);
					weight_x[6] = (15*dx*dz*dz/r7) - (3*dx/r5);
					weight_x[7] = (15*dx*dx*dy/r7) - (3*dy/r5);
					weight_x[8] = (15*dx*dy*dz/r7);
					weight_x[9] = (15*dx*dx*dz/r7) - (3*dz/r5);
					// f = dy/r^3
					weight_y[0] = dy/r3;
					weight_y[1] = -(3*dy*dx/r5);
					weight_y[2] = -(3*dy*dy/r5) + (1/r3);
					weight_y[3] = -(3*dy*dz/r5);
					weight_y[4] = (15*dy*dx*dx/r7) - (3*dy/r5);
					weight_y[5] = (15*dy*dy*dy/r7) - (9*dy/r5);
					weight_y[6] = (15*dy*dz*dz/r7) - (3*dy/r5);
					weight_y[7] = (15*dy*dx*dy/r7) - (3*dx/r5);
					weight_y[8] = (15*dy*dy*dz/r7) - (3*dz/r5);
					weight_y[9] = (15*dy*dx*dz/r7);
					// f = dz/r^3
					weight_z[0] = dz/r3;
					weight_z[1] = -(3*dz*dx/r5);
					weight_z[2] = -(3*dz*dy/r5);
					weight_z[3] = -(3*dz*dz/r5) + (1/r3);
					weight_z[4] = (15*dz*dx*dx/r7) - (3*dz/r5);
					weight_z[5] = (15*dz*dy*dy/r7) - (3*dz/r5);
					weight_z[6] = (15*dz*dz*dz/r7) - (9*dz/r5);
					weight_z[7] = (15*dz*dx*dy/r7);
					weight_z[8] = (15*dz*dy*dz/r7) - (3*dy/r5);
					weight_z[9] = (15*dz*dx*dz/r7) - (3*dx/r5);

					double K = -1/(4.0e0*Parameters::pi);
					for(size_t j = 0; j < 10; j++)
					{
						pi.u[idx] += K * (cells.multipole_z[j][c] * weight_y[j] - cells.multipole_y[j][c] * weight_z[j]);
						pi.v[idx] += K * (cells.multipole_x[j][c] * weight_z[j] - cells.multipole_z[j][c] * weight_x[j]);
						pi.w[idx] += K * (cells.multipole_y[j][c] * weight_x[j] - cells.multipole_x[j][c] * weight_y[j]);
					}
				}	
			}
		}
	}
	// # leaf cell
	else
	{
		// # loop in twig cell's particles
		for(int j = 0; j < cells.nleaf[p]; j++)
		{
			int is = cells.leaf[j][p];
			double dx,dy,dz,rij2,rij,rij3;
			dx = pi.x[idx] - pi.x[is];
			dy = pi.y[idx] - pi.y[is];
			dz = pi.z[idx] - pi.z[is];
			rij2 = std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2);
			rij = std::sqrt(rij2);
			rij3 = rij2 * rij;

			if (rij2 > 0.0e0)
			{
				double sij2,sij,gsig;
				sij2 = (std::pow(pi.s[idx], 2) + std::pow(pi.s[is], 2)) * 0.25e0;
				sij = std::sqrt(sij2);
				regul_func_3d(rij, sij, gsig);

				pi.u[idx] += - (dy * pi.gz[is] - dz * pi.gy[is]) * gsig / rij3;
				pi.v[idx] += - (dz * pi.gx[is] - dx * pi.gz[is]) * gsig / rij3;
				pi.w[idx] += - (dx * pi.gy[is] - dy * pi.gx[is]) * gsig / rij3;
			}
		}
	}
}
