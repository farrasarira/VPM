#include "parameters.hpp"
#include "particle.hpp"
#include "initialization.hpp"
#include "poisson.hpp"
#include "stretching.hpp"	
#include "diffusion.hpp"
#include "save_data.hpp"
#include <time.h>
#include <stdio.h>
#include <algorithm>


#include <iostream>
#include <numeric>
#include <fstream>      // std::ofstrea

int main(int argc, char const *argv[]){

    Particle particle; // create simulation particle object
	Body body;		   // create obstacle object
	initialization initialization_step;
	poisson_solver advection_step;
	stretching stretching_step;
	diffusion diffusion_step;
	save_data save_step;

	// TODO: Initial Particles Generation
	std::cout << "Domain initialization ..." << std::endl;
	initialization_step.init_vortex_ring(particle);
	// initialization_step.init_2D_Test_Domain(particle);
	std::cout << "Domain initialization (DONE)" << std::endl;
	
	int nt_start = 0;
	double cum_time = 0.0; // actual cumulative time
	std::vector<double> _cumulativeTime;
	clock_t delta_t; // actual time difference

    for(int it = nt_start; it < Parameters::nt; ++it)
    {
		printf("+--------------- iter no. %d -------------------+\n", (int)it);

		// TODO: Poisson: solving Rotational Velocity, Stretching
		advection_step.poisson(particle);

		// TODO: Diffusion and Stretching Sub-step
		std::vector<std::vector<double>> _dfdtStr(3,std::vector<double>(particle.num));
		stretching_step.main_stretching(particle, _dfdtStr); 
		std::vector<std::vector<double>> _dfdtDiff(3,std::vector<double>(particle.num));
		diffusion_step.main_diffusion(particle, _dfdtDiff); 

		// Update Particle (Euler Scheme)
		for (size_t i = 0; i < particle.num; i++)
		{
        	particle.x[i] += Parameters::dt * (particle.u[i]);
        	particle.y[i] += Parameters::dt * (particle.v[i]);
        	particle.z[i] += Parameters::dt * (particle.w[i]);		
			// std::cout << particle.vort_x[i] << " | " << (_dfdtStr[0][i]) << std::endl;
			if (_dfdtStr[0][i] < 2 || _dfdtStr[1][i] < 2 || _dfdtStr[2][i] < 2 || _dfdtDiff[0][i] < 2 || _dfdtDiff[1][i] < 2 || _dfdtDiff[2][i] < 2){
				particle.vort_x[i] += Parameters::dt * (_dfdtDiff[0][i]+_dfdtStr[0][i]); // 
				particle.vort_y[i] += Parameters::dt * (_dfdtDiff[1][i]+_dfdtStr[1][i]); // 
				particle.vort_z[i] += Parameters::dt * (_dfdtDiff[2][i]+_dfdtStr[2][i]); // 
				particle.gx[i] = particle.vort_x[i] * particle.s[i] * particle.s[i] * particle.s[i];
				particle.gy[i] = particle.vort_y[i] * particle.s[i] * particle.s[i] * particle.s[i];
				particle.gz[i] = particle.vort_z[i] * particle.s[i] * particle.s[i] * particle.s[i];
			}
		}

		// TODO: Saving data
		// delta_t = clock() - delta_t;							   // actual time [s]
		// _cumulativeTime.push_back(static_cast<double>(delta_t) / static_cast<double>(CLOCKS_PER_SEC));
		// cum_time = accumulate(_cumulativeTime.begin(), _cumulativeTime.end(), 0.0e0);
		// cum_time = cum_time + (double)delta_t / CLOCKS_PER_SEC; // cumulative time

		int _numberOfActiceParticle = 0;
		for (int i = 0; i < particle.num; i++)
		{
			if (particle.isActive[i])
			{
				_numberOfActiceParticle += 1;
			}
		}

		std::ofstream outs;
		if (it == 0)
		{
			outs.open("output/particle_number_time.dat");
			outs << it << "," << it * Parameters::dt << "," << cum_time << "," << _numberOfActiceParticle << "\n";
			outs.close();
		}
		else if (it >= 1)
		{
			outs.open("output/particle_number_time.dat", std::ofstream::out | std::ofstream::app);
			outs << it << "," << it * Parameters::dt << "," << cum_time << "," << _numberOfActiceParticle << "\n";
			outs.close();
		}
		std::cout << "\nnumber of particles: " << _numberOfActiceParticle << std::endl;
		// std::cout << "\ncomputation time: " << cum_time << std::endl;
		save_step.output(it, particle, body, cum_time);

		// if(it == nt_start)
		// 	break;



    }

}