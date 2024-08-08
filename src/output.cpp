#include "save_data.hpp"
#include "parameters.hpp"
#include <vector>
#include <iostream>

save_data::save_data()
{
	// ---initializing internal variables inside 'save_data' class
	// into their respective size and/or value
}

save_data::~save_data()
{
	// ---deallocating memory from internal variables
}

void save_data::output(int it, Particle &p, Body &b, double &cum_time)
{
	// -- accessing struct data, variables that are not commented are an input only
	int &np = p.num;
	std::vector<double> &xp = p.x;
	std::vector<double> &yp = p.y;
	std::vector<double> &zp = p.z;
	std::vector<double> &gpx = p.vort_x;
	std::vector<double> &gpy = p.vort_y;
	std::vector<double> &gpz = p.vort_z;
	std::vector<double> &sp = p.s;
	std::vector<double> &up = p.u;
	std::vector<double> &vp = p.v;
	std::vector<double> &wp = p.w;

	// internal variables
	std::string name1, name2, name3, name4;

	printf("Saving data ...\n");

	
    // TODO : Save Energy Calc
    energy_calc(it, np, xp, yp, zp, up, vp, wp, gpx, gpy, gpz, sp); // ! Should calculate all domain and delete variable gamtrim  


	// TODO: Saving common data and extra data
	if (true)
	{
		// defining output variables
		std::ofstream ofs;
		std::stringstream ss;

		cum_time = 0.0e0; // Reset for cumulative time saving file
		// Giving name and opening the files which want to save.
		// Giving "number of name" = "time step"
		if (true)
		{
			ss << std::setw(5) << std::setfill('0') << it; //because c++ indexing is started from 0
		}
		else
		{
			ss << std::setw(5) << std::setfill('0') << it - 1; //because c++ indexing is started from 0
		}

		std::string nData = ss.str();
		// If not incase of "later saving grid data", (if not do the file can be rewired)
		if (true) // save only particles, or save particles and grid at same time
		{
			std::cout << "masuk pakk" << std::endl;
			// Save common data (information of free particles):
			name1.append("output/x_vor_vel_particle_");
			name1.append(nData);
			name1.append(".csv");
			ofs.open(name1.c_str());
			ofs << "" << "xp" << "," << "yp" << "," << "zp" << "," << "vortx" << "," << "vorty"<< "," << "vortz" << "," << "up" << "," << "vp"<< "," << "wp" << "," << "sp\n";
			for (int i = 0; i < np; i++)
			{
				if (p.isActive[i] == true)
				{
					ofs << "" << xp[i]
						<< "," << yp[i]
						<< "," << zp[i]
						<< "," << gpx[i]
						<< "," << gpy[i]
						<< "," << gpz[i]
						<< "," << up[i] 
						<< "," << vp[i] 
						<< "," << wp[i]
						<< "," << sp[i]
						<< "\n";
				}
			}
			ofs.close();
		}

	}

} // end of function
