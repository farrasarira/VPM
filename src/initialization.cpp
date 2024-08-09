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

void initialization::init_vortex_ring(Particle &p)
{    
    const int nT = 200;
    const int nc = 3;

    const double Gamma = 1.0; // circulation
    const double Rad = 1.0;
    const double rad = 0.05*Rad;
    const double rl = rad / (2.0*nc+1.0);
    const double cell_size = sqrt(Parameters::pi)*rl;
    const double deltaT = 2*Parameters::pi / nT;

    const double x_center = 0.0;
    const double y_center = 0.0;
    const double z_center = 0.0;

    Particle _particle;
    int _nparticle = 0;

    for (int iter_T = 0; iter_T < nT; ++iter_T){
        double Theta = iter_T * deltaT;
        for (int iter_c = 0; iter_c <= nc; ++iter_c){
            double xpos = 0.0;
            double ypos = 0.0;
            double zpos = 0.0;
            double vort_mag = 0.0;
            double xvort = 0.0;
            double yvort = 0.0;
            double zvort = 0.0;

            if (iter_c == 0){
                xpos = x_center + Rad*cos(Theta);
                ypos = y_center + Rad*sin(Theta);
                zpos = z_center;

                vort_mag = Gamma/(2.0*Parameters::pi*rad*rad);

                xvort = vort_mag * sin(Theta);
                yvort = vort_mag * -1*cos(Theta);
                zvort = 0.0;

                _particle.x.push_back(xpos);
                _particle.y.push_back(ypos);
                _particle.z.push_back(zpos);
                _particle.s.push_back(cell_size);
                _particle.gx.push_back(xvort*cell_size*cell_size*cell_size);
                _particle.gy.push_back(yvort*cell_size*cell_size*cell_size);
                _particle.gz.push_back(zvort*cell_size*cell_size*cell_size); 
                _particle.vort_x.push_back(xvort);
                _particle.vort_y.push_back(yvort);
                _particle.vort_z.push_back(zvort);                 
                _particle.u.push_back(0.0e0);
                _particle.v.push_back(0.0e0);
                _particle.w.push_back(0.0e0);
                _nparticle++;
            }
            else{
                double rc = rl * (1.0+12.0*iter_c*iter_c)/(6.0*iter_c);
                int n_phi = iter_c * 4;
                double delta_phi = 2.0*Parameters::pi / n_phi;
                for (int iter_phi = 1; iter_phi <= n_phi; ++iter_phi){
                    double phi = iter_phi * delta_phi;

                    xpos = x_center + Rad*cos(Theta) + rc*cos(Theta)*cos(phi);
                    ypos = y_center + Rad*sin(Theta) + rc*sin(Theta)*cos(phi);
                    zpos = z_center + rc*sin(phi);

                    vort_mag = Gamma/(2.0*Parameters::pi*rad*rad) * exp(-1.0*rc*rc/(2.0*rad*rad));

                    xvort = vort_mag * sin(Theta);
                    yvort = vort_mag * -1*cos(Theta);
                    zvort = 0.0;

                    _particle.x.push_back(xpos);
                    _particle.y.push_back(ypos);
                    _particle.z.push_back(zpos);
                    _particle.s.push_back(cell_size);
                    _particle.gx.push_back(xvort*cell_size*cell_size*cell_size);
                    _particle.gy.push_back(yvort*cell_size*cell_size*cell_size);
                    _particle.gz.push_back(zvort*cell_size*cell_size*cell_size); 
                    _particle.vort_x.push_back(xvort);
                    _particle.vort_y.push_back(yvort);
                    _particle.vort_z.push_back(zvort);                 
                    _particle.u.push_back(0.0e0);
                    _particle.v.push_back(0.0e0);
                    _particle.w.push_back(0.0e0);
                    _nparticle++;
                }
            }

        }
    }
    _particle.num = _nparticle;
    _particle.isActive.resize(_nparticle, true);

    _particle.neighbor = d_neighbor.direct_find(_particle.num, _particle.s, _particle.x, _particle.y, _particle.z, Parameters::r_scale);

    p = _particle;
}

void initialization::init_2D_Test_Domain(Particle &p)
{    
    const double lx = 1.0;
    const double ly = 1.0;

    const int nx = 100;
    const int ny = 100;

    const double blob_size = 1.0/nx;

    Particle _particle;
    int _nparticle = 0;

    for (int iter_x = 0; iter_x <= nx; ++iter_x){
        for (int iter_y = 0; iter_y <= ny; ++iter_y){
            double xpos = blob_size*iter_x;
            double ypos = blob_size*iter_y;
            double zpos = 0.0;
            double vort_mag = 0.0;
            double xvort = 0.0;
            double yvort = 0.0;
            double zvort = 0.0;
            double ux = 0.1*pow(xpos,2)+0.2*pow(ypos,2)+0.3*xpos*ypos+0.4*xpos+0.5*ypos+0.6;
            double uy = 0.0;
            double uz = 0.0;

            _particle.x.push_back(xpos);
            _particle.y.push_back(ypos);
            _particle.z.push_back(zpos);
            _particle.s.push_back(blob_size);
            _particle.gx.push_back(xvort*blob_size*blob_size*blob_size);
            _particle.gy.push_back(yvort*blob_size*blob_size*blob_size);
            _particle.gz.push_back(zvort*blob_size*blob_size*blob_size); 
            _particle.vort_x.push_back(xvort);
            _particle.vort_y.push_back(yvort);
            _particle.vort_z.push_back(zvort);                 
            _particle.u.push_back(ux);
            _particle.v.push_back(uy);
            _particle.w.push_back(uz);
            _nparticle++;
            

        }
    }
    _particle.num = _nparticle;
    _particle.isActive.resize(_nparticle, true);

    _particle.neighbor = d_neighbor.direct_find(_particle.num, _particle.s, _particle.x, _particle.y, _particle.z, Parameters::r_scale);

    p = _particle;
}
