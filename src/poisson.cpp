#include "poisson.hpp"

#include <iostream>

// =================================================================
//  TODO: solving POISSON EQUATION for Particles by using Biot_Savart and technique tree-code to accelerate computation speed
//  input: PARTICLE info: vorticity-strength, position
//  output: calculate velocity 
// =================================================================

void poisson_solver::poisson(Particle &p)
{
	// -- reset velocity
	p.u.clear();
	p.u.resize(p.num);
	p.v.clear();
	p.v.resize(p.num);
	p.w.clear();
	p.w.resize(p.num);

	// -- internal variables
	const static int n_inter = 1; // ! change into class static variables
	std::vector<int> _index;
	Particle _particleActive; // store current active particles
	Particle _particle;		  // store current active particles box (to increase robustness)
	Particle _particleDense;  // store particle with core size as smallest from current active particles
	_particleActive.num = 0;
	_particle.num = 0;
	_particleDense.num = 0;

	// TODO: store current active particles
	for (int i = 0; i < p.num; i++)
	{
		if (p.isActive[i] == true)
		{
			_particle.x.emplace_back(p.x[i]);
			_particle.y.emplace_back(p.y[i]);
			_particle.z.emplace_back(p.z[i]);
			_particle.s.emplace_back(p.s[i]);
			// _particle.s.push_back(p.s[i] * Pars::er);
			// _particle.gz.push_back(p.gz[i] * (Pars::er * std::pow(p.s[i],2)) * Pars::er * std::pow(p.s[i] / Pars::sigma, 2));  //Convert Vorticity into Strength
			_particle.gx.push_back(p.gx[i]);
			_particle.gy.push_back(p.gy[i]);
			_particle.gz.push_back(p.gz[i]); 
			//_particle.gx.emplace_back(p.gx[i] * std::pow(p.s[i] / Pars::sigma, 3));
			//_particle.gy.emplace_back(p.gy[i]* std::pow(p.s[i] / Pars::sigma, 3));
			//_particle.gz.emplace_back(p.gz[i]* std::pow(p.s[i] / Pars::sigma, 3));
			// _particle.gz.push_back(p.gz[i]);
			_particle.u.emplace_back(0.0e0);
			_particle.v.emplace_back(0.0e0);
			_particle.w.emplace_back(0.0e0);
			_index.push_back(i);
			_particle.neighbor.push_back(p.neighbor[i]);
			_particle.num++;
		}
	}

	// TODO: calculate velocity using treecode 
	// d_base_poisson_3d.biotsavart_treecode_3d(_particle,_particle);
	
	// // TODO: calculate velocity using direct sum
	d_base_poisson_3d.biotsavart_direct_3d(_particle,_particle);

	// update base particle velocity
	for (size_t i = 0; i < _particle.num; i++)
	{
		int _isActiveIndex = _index[i];
		p.u[_isActiveIndex] = _particle.u[i];
		p.v[_isActiveIndex] = _particle.v[i];
		p.w[_isActiveIndex] = _particle.w[i];		
	}
	// !STEP 4: end
}
