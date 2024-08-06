#include "particle.hpp"
#include "initialization.hpp"
#include "parameters.hpp"
#include <iostream>


int main(int argc, char const *argv[]){
    int nt_start = 0;
	double cum_time = 0.0; // actual cumulative time

    Particle particle; // create simulation particle object

	initialization initialization_step;

    for(size_t it = nt_start; it < Parameters::nt; ++it)
    {
		printf("+--------------- iter no. %d -------------------+\n", (int)it);
        
        if(it == 0)
        {
            for (size_t ii = 0; ii < particle.x.size();ii++)
            {
				if (particle.x[ii] <= -Parameters::xdom - 10 * Parameters::sigma || particle.x[ii]  >= ((Parameters::lxdom-Parameters::xdom) + 6 * Parameters::sigma)
					 || particle.y[ii]  <= (-Parameters::lydom/2) - 10 * Parameters::sigma || particle.y[ii]  >= ((Parameters::lydom/2) + 6 * Parameters::sigma ))
            	{
                	particle.isboundary.push_back(1);
					particle.boundaryval.push_back(0.0);
            	}
            	else
            	{
               		particle.isboundary.push_back(0);
					particle.boundaryval.push_back(0.0);
            	}
			}
        }

        if (it == 0){
			particle.vorticity.resize(particle.num, 0.0e0);
		}else{
			for (size_t i = 0; i < particle.num; i++)
			{
				particle.vorticity[i] = particle.gz[i] / pow(particle.s[i],2);
			}
		}

        if (it == 0){
			for (int j = 0; j < particle.x.size(); j++){
				particle.label.push_back(j);
			} 
		}

        // TODO: Poisson: solving Rotational Velocity, Stretching
		advection_step.poisson(particle, it);	

		penalization_step.get_penalization(particle, body, it);

        // TODO: Convection Sub-step
		std::vector<double> _dfdtConv(particle.num);
		advection_step.advection(particle, _dfdtConv); // ! later: do 2nd order scheme

        //mulai step 2
		// TODO: Diffusion Sub-step
		std::vector<double> _dfdtDiff(particle.num);
		diffusion_step.main_diffusion(particle, _dfdtDiff); // ! later: do 2nd order scheme
		
		// time integration (1st order; Diffusion)
		for (size_t i = 0; i < particle.num; i++)
		{
			particle.gz[i] += Pars::dt * (_dfdtDiff[i]); 
		}
		
		//TODO: Move body (deformation, translation, rotation)
		geom_step.moving_body(it, body);
		//selesai step 2 dan ulangi

        // TODO: Saving data
		//delta_t = clock() - delta_t;							   // actual time [s]
		std::chrono::duration<double> elapsed_time_ms = (std::chrono::system_clock::now() - t_start);
		//_cumulativeTime.push_back(static_cast<double>(delta_t) / static_cast<double>(CLOCKS_PER_SEC));
		_cumulativeTime.push_back(elapsed_time_ms.count());
		cum_time = accumulate(_cumulativeTime.begin(), _cumulativeTime.end(), 0.0e0);
		// cum_time = cum_time + (double)delta_t / CLOCKS_PER_SEC; // cumulative time
		/*
		int _numberOfActiceParticle = 0;
		for (size_t i = 0; i < particle.num; i++)
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
			outs << it << "," << it * Pars::dt << "," << cum_time << "," << _numberOfActiceParticle << "\n";
			outs.close();
		}
		else if (it >= 1)
		{
			outs.open("output/particle_number_time.dat", std::ofstream::out | std::ofstream::app);
			outs << it << "," << it * Pars::dt << "," << cum_time << "," << _numberOfActiceParticle << "\n";
			outs.close();
		}
		std::cout << "\nnumber of particles: " << _numberOfActiceParticle << std::endl;*/
		std::cout << "\ncomputation time: " << cum_time << std::endl;
		//if (it % Pars::nt_sf == 0)
		//{
		save_step.output(it, particle, body, cum_time);



    }

}