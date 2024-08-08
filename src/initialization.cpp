#include "initialization.hpp"
#include "parameters.hpp"
#include "particle.hpp"

#include <math.h>
#include <iostream>

initialization::initialization(/* args */)
{
}

initialization::~initialization()
{
}


void initialization::init_domain(Particle &p)
{
    const static double lxdom = Parameters::lxdom;   // ! should be put on global paremeters
    const static double lydom = Parameters::lydom;   // ! should be put on global paremeters
    const static double lzdom = Parameters::lzdom;   // ! should be put on global paremeters
    
    int _nx = std::ceil(lxdom / (Parameters::sigma * 1.0)) + 1;
    int _ny = std::ceil(lydom / Parameters::sigma) + 1;
    int _nz = std::ceil(lzdom / Parameters::sigma) + 1;
 
    Particle _particle;
    int _nparticle = 0;
    // TODO: generate particles
    for (size_t k = 0; k < _nz; k++)
    {
        for (size_t i = 0; i < _nx; i++)
        {
            for (size_t j = 0; j < _ny; j++)
            {
                double _x = -Parameters::xdom +  static_cast<double>(i) * Parameters::sigma * 1.0;
                double _y = -lydom / 2 + static_cast<double>(j) * Parameters::sigma;
                double _z = -lzdom / 2 + static_cast<double>(k) * Parameters::sigma;
                // double _x = -lxdom / 2 + static_cast<double>(i) * Parameters::D_0;
                // double _y = -lydom / 2 + static_cast<double>(j) * Parameters::D_0;
                // double _z = -lzdom / 2 + static_cast<double>(k) * Parameters::D_0;
                _particle.x.push_back(_x);
                _particle.y.push_back(_y);
                _particle.z.push_back(_z);
                _particle.s.push_back(Parameters::sigma);
                _particle.gx.push_back(0.0e0);
                _particle.gy.push_back(0.0e0);
                _particle.gz.push_back(0.0e0); 
                _particle.vort_x.push_back(0.0e0);
                _particle.vort_y.push_back(0.0e0);
                _particle.vort_z.push_back(0.0e0);                 
                _particle.u.push_back(0.0e0);
                _particle.v.push_back(0.0e0);
                _particle.w.push_back(0.0e0);
                _nparticle++;
            }
        }
    }    
    _particle.num = _nparticle;
    _particle.isActive.resize(_nparticle, false);

    // // * set active particles
    // int _num = 0;
    // for (size_t i = 0; i < _nparticle; i++)
    // {
    //     if (_particle.x[i]>= (-lxdom/4) && _particle.x[i] <= (lxdom/4) &&
    //         _particle.y[i]>= (-lydom/4) && _particle.y[i] <= (lydom/4) &&
    //         _particle.z[i]>= (-lzdom/4) && _particle.z[i] <= (lzdom/4)
    //     )
    //     {
    //         _particle.isActive[i] = true;
    //         _num++;
    //     }
    // }

    // TODO: generate neighbor
    // neighbor searching
    // ---direct find
    // if (Parameters::neighbor_opt == 0)
    // {
      _particle.neighbor = d_neighbor.direct_find(_particle.num, _particle.s, _particle.x, _particle.y, _particle.z, Parameters::r_scale);
    // }
    // // ---linked list
    // else if (Parameters::neighbor_opt == 1)
    // {
        // _particle.neighbor = d_neighbor.link_list(_particle.num, _particle.s, _particle.x, _particle.y, _particle.z, Parameters::r_scale);
    // }

    // cout << "initial particle number: " << _num << " : " << _nparticle << endl;

    // * assign into parameter
    p = _particle;

    /*std::ofstream ofs;
    ofs.open("output/init_particle.csv");
		ofs << "" << "xp" << "," << "yp" << "," << "zp"<< "," << "gpx" << "," <<"gpy" << "," <<"gpz" << "," << "up" << "," << "vp" << "," << "wp" << "," << "sp\n";
		for (int i = 0; i < p.x.size(); i++){
			ofs << "" << p.x[i]<< "," << p.y[i]<< "," << p.z[i] << "," << p.gx[i]<< "," <<p.gy[i]<< "," <<p.gz[i]<< "," << p.u[i]<< "," << p.v[i]<< "," << p.w[i]<< "," << p.s[i]<< "\n";
		}
	ofs.close();*/

    // std::ofstream ofs1;
    // ofs1.open("output/cek_neighbor.csv");
	// 	int np = p.x.size();
    //     int midp = std::ceil(np/2);
    //     ofs1 << "" << "xp" << "," << "yp" << "," << "zp\n";    
    //     ofs1 << p.x[midp] <<"," <<p.y[midp]<< ","<<p.z[midp]<<"\n";
	// 	for (int i = 0; i < p.neighbor[midp].size(); i++){
    //         int index_midp = p.neighbor[midp][i];
	// 		ofs1 << "" << p.x[index_midp]<< "," << p.y[index_midp]<<","<< p.z[index_midp]<< "\n";
	// 	}
	// ofs1.close();

}

void initialization::init_vorticity(Particle &p)
{
    // define vortex ring 1
    double R1 = Parameters::ly;
    double Gamma1 = 1.0e0;
    double sigma1 = 0.05 * Parameters::ly;
    double sigma12 = std::pow(sigma1,2);
    double *xcen_pen1 = new double[3];
    xcen_pen1[0] = 0.50 * Parameters::lx;
    xcen_pen1[1] = 0.0e0;
    xcen_pen1[2] = 0.50 * Parameters::lz;
     
    // define vortex ring 2
    double R2 = Parameters::ly;
    double Gamma2 = 1.0e0;
    double sigma2 = 0.05 * Parameters::ly;
    double sigma22 = std::pow(sigma2,2);
    double *xcen_pen2 = new double[3];
    xcen_pen2[0] = -(0.50 * Parameters::lx);
    xcen_pen2[1] = 0.0e0;
    xcen_pen2[2] = -(0.50 * Parameters::lz);
    
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        
      if(p.z[i] >= (xcen_pen2[2]-1.0e0) && p.z[i] <= (xcen_pen1[2]+1.0e0))
        {
        double xpi,ypi,zpi; // Coordinates relative to center of ring
        double phi;
        double xci,yci,zci; // Coordinates core center relative to center of ring
        double xn,yn,zn,gxn,gyn,gzn; // Coordinates transformation
        double gtheta, theta, rci, rci2,rcxy2,rcxy;
        double gpx1,gpy1,gpz1;
        double gpx2,gpy2,gpz2;
        double Vp = std::pow(p.s[i],3);


        // vortex ring 1 
        if(p.z[i] >= 0.0e0){
        
        xpi = p.x[i] - xcen_pen1[0];
        ypi = p.y[i] - xcen_pen1[1];
        zpi = p.z[i] - xcen_pen1[2];
        phi = atan2(ypi,xpi);
        xci = R1 * cos(phi);
        yci = R1 * sin(phi);
        zci = 0.0e0;
        rci2 = std::pow((xpi-xci),2) + std::pow((ypi-yci),2) + std::pow((zpi-zci),2);
        rci = std::sqrt(rci2);
        
        gtheta = (Gamma1/(2.0e0*Parameters::pi*sigma12))*exp(-rci2/(2.0e0*sigma12));
        
        theta = phi;
        gpx1 = gtheta * sin(theta);
        gpy1 = - gtheta * cos(theta);
        gpz1 = 0.0e0;

        // xn = xpi * cos(phi) - ypi * sin(phi);
        // yn = xpi * sin(phi) + ypi * cos(phi);
        // zn = zpi;
        // theta = atan2(zn,xn);
        
        // gxn = -gtheta * sin(theta);
        // gzn = gtheta * cos(theta);
        // gyn = 0.0e0;

        // gpx1 =  gxn * cos(phi) + gyn * sin(phi);
        // gpy1 =  - gxn * sin(phi) + gyn * cos(phi);
        // gpz1 =  gzn;

        // rcxy2 = std::pow((xpi-xci),2) + std::pow((ypi-yci),2);
        // rcxy = std::sqrt(rcxy2);
        // theta = atan2(rcxy,(zpi-zci));
        // if(rci != 0.0e0){
        // theta = acos((zpi-zci)/rci);

        // gpx1 = gtheta * cos(theta) * cos(phi);
        // gpy1 = gtheta * cos(theta) * sin(phi);
        // gpz1 = gtheta * sin(theta); 

        // }    
        // else{
        // gtheta = (Gamma1/(2.0e0*Parameters::pi*sigma12))*exp(-rci2/(2.0e0*sigma12));
        // gpx1 = 0.0e0;
        // gpy1 = 0.0e0;
        // gpz1 = gtheta;         
        // }

        // p.gx[i] += gpx1;
        // p.gy[i] += gpy1;
        // p.gz[i] += gpz1;
        p.gx[i] = gpx1 * Vp;
        p.gy[i] = gpy1 * Vp;
        p.gz[i] = gpz1 * Vp;
        
        }
        else{
 
        // vortex ring 2
        xpi = p.x[i] - xcen_pen2[0];
        ypi = p.y[i] - xcen_pen2[1];
        zpi = p.z[i] - xcen_pen2[2];
        phi = atan2(ypi,xpi);
        xci = R2 * cos(phi);
        yci = R2 * sin(phi);
        zci = 0.0e0;
        rci2 = std::pow((xpi-xci),2) + std::pow((ypi-yci),2) + std::pow((zpi-zci),2);
        rci = std::sqrt(rci2);
        
        gtheta = (Gamma2/(2.0e0*Parameters::pi*sigma22))*exp(-rci2/(2.0e0*sigma22));

        theta = phi;
        gpx2 = - gtheta * sin(theta);
        gpy2 = gtheta * cos(theta);
        gpz2 = 0.0e0;

        // xn = xpi * cos(phi) - ypi * sin(phi);
        // yn = xpi * sin(phi) + ypi * cos(phi);
        // zn = zpi;
        // theta = atan2(zn,xn);

        // gxn = -gtheta * sin(theta);
        // gzn = gtheta * cos(theta);
        // gyn = 0.0e0;

        // gpx2 =  gxn * cos(phi) + gyn * sin(phi);
        // gpy2 =  - gxn * sin(phi) + gyn * cos(phi);
        // gpz2 =  gzn;

        // rcxy2 = std::pow((xpi-xci),2) + std::pow((ypi-yci),2);
        // rcxy = std::sqrt(rcxy2);
        // theta = atan2(rcxy,(zpi-zci));
        // if(rci != 0.0e0){
        // theta = acos((zpi-zci)/rci);

        // gpx2 = gtheta * cos(theta) * cos(phi);
        // gpy2 = gtheta * cos(theta) * sin(phi);
        // gpz2 = gtheta * sin(theta); 
        // }    
        // else{
        // gtheta = (Gamma1/(2.0e0*Parameters::pi*sigma12))*exp(-rci2/(2.0e0*sigma12));
        // gpx1 = 0.0e0;
        // gpy1 = 0.0e0;
        // gpz1 = gtheta;         
        // }
        
        // p.gx[i] += gpx2;
        // p.gy[i] += gpy2;
        // p.gz[i] += gpz2;
        p.gx[i] = gpx2 * Vp;
        p.gy[i] = gpy2 * Vp;
        p.gz[i] = gpz2 * Vp;         
 
        }
      }
    }   

    // // * prevent truncated error
    // double maxx_temp = 0.0e0;    double maxy_temp = 0.0e0;   double maxz_temp = 0.0e0;
    // double gpx_max = maxx_temp;  double gpy_max = maxy_temp; double gpz_max = maxz_temp;
    
    // for (int i = 0; i < p.num; i++)
    // {
    //     gpx_max = maxx_temp > std::abs(p.gx[i]) ? maxx_temp : std::abs(p.gx[i]);
    //     maxx_temp = gpx_max;
    //     gpy_max = maxy_temp > std::abs(p.gy[i]) ? maxy_temp : std::abs(p.gy[i]);
    //     maxy_temp = gpy_max;
    //     gpz_max = maxy_temp > std::abs(p.gz[i]) ? maxz_temp : std::abs(p.gz[i]);
    //     maxz_temp = gpz_max;
    // }

    // for (size_t i = 0; i < p.num; i++)
    // {
    //     if (p.gx[i] < 1.0e-4 * gpx_max){  
    //         p.gx[i] = 0.0e0; // prevent truncated error   
    //     }
    //     if (p.gy[i] < 1.0e-4 * gpy_max){  
    //         p.gy[i] = 0.0e0; // prevent truncated error   
    //     }
    //     if (p.gz[i] < 1.0e-4 * gpz_max){  
    //         p.gz[i] = 0.0e0; // prevent truncated error   
    //     }
    // }
    /*
    std::ofstream ofs;
    ofs.open("output/init_vorticity.csv");
		ofs << "" << "xp" << "," << "yp" << "," << "zp"<< "," << "gpx" << "," <<"gpy" << "," <<"gpz" << "," << "up" << "," << "vp" << "," << "wp" << "," << "sp\n";
		for (int i = 0; i < p.x.size(); i++){
			ofs << "" << p.x[i]<< "," << p.y[i]<< "," << p.z[i] << "," << p.gx[i]<< "," <<p.gy[i]<< "," <<p.gz[i]<< "," << p.u[i]<< "," << p.v[i]<< "," << p.w[i]<< "," << p.s[i]<< "\n";
		}
	ofs.close();*/

}

void initialization::initializeVortexRing(Particle &particle, int num_particles, double ring_radius, double ring_thickness, double core_size, double vortex_strength)
{
    // Clear existing data
    particle.num = num_particles;
    particle.x.resize(num_particles);
    particle.y.resize(num_particles);
    particle.z.resize(num_particles);
    particle.s.resize(num_particles);
    particle.u.resize(num_particles);
    particle.v.resize(num_particles);
    particle.w.resize(num_particles);
    particle.gx.resize(num_particles);
    particle.gy.resize(num_particles);
    particle.gz.resize(num_particles);
    particle.vort_x.resize(num_particles);
    particle.vort_y.resize(num_particles);
    particle.vort_z.resize(num_particles);

    double angle_step = 2.0 * M_PI / num_particles; // Angular step to distribute particles around the ring
    for (int i = 0; i < num_particles; ++i)
    {
        double angle = i * angle_step;
        
        // Particle coordinates
        particle.x[i] = ring_radius * cos(angle);
        particle.y[i] = ring_radius * sin(angle);
        particle.z[i] = 0.0; // Assume the vortex ring lies in the xy-plane

        // Core size
        particle.s[i] = core_size;

        // Initialize velocities (assuming a simple vortex ring)
        double u_vortex = -vortex_strength * sin(angle);
        double v_vortex = vortex_strength * cos(angle);
        double w_vortex = 0.0; // No z-component for a planar vortex ring

        particle.u[i] = u_vortex;
        particle.v[i] = v_vortex;
        particle.w[i] = w_vortex;

        // Strengths
        particle.gx[i] = vortex_strength * cos(angle);
        particle.gy[i] = vortex_strength * sin(angle);
        particle.gz[i] = 0.0; // No z-component for a planar vortex ring

        // Vorticity (consistent with a simple vortex ring)
        particle.vort_x[i] = 0.0;
        particle.vort_y[i] = 0.0;
        particle.vort_z[i] = vortex_strength; // Vorticity in z-direction for a planar ring
    }

    std::cout << "Initialized vortex ring with " << num_particles << " particles." << std::endl;
}