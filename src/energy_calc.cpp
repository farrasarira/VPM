#include "parameters.hpp"
#include "save_data.hpp"

#include <math.h>

void save_data::energy_calc(int it, int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &zp,
                              const std::vector<double> &up, const std::vector<double> &vp, const std::vector<double> &wp, 
                              const std::vector<double> &gpx, const std::vector<double> &gpy, const std::vector<double> &gpz,
                              const std::vector<double> &sp)
{
    // internal variables
    double Ek, Es; // Ek = kinetic energy and Es = enstrophy
    Ek = 0.0e0;
    Es = 0.0e0;

    for (int i = 0; i < np; i++)
    {
        double Vp1 = std::pow(sp[i],3) * 1.0;
        // Kinetic Energy
        double Ut = std::pow(up[i],2) + std::pow(vp[i],2) + std::pow(wp[i],2); 
        Ek += 0.5 * Ut * Vp1;

        // Enstrophy
        double Vp2 = 1.0e0;
        // double Vp2 = std::pow(sp[i] / Pars::sigma, 3);
        double gx,gy,gz;
        gx = gpx[i]/Vp1; gy = gpy[i]/Vp1; gz = gpz[i]/Vp1;
        double gp = std::pow(gx,2) + std::pow(gy,2) + std::pow(gz,2);  
        Es += gp * Vp1 * Vp2;
    }

    // saving data
    std::ofstream ofs;
    if (it == 0)
    {
        ofs.open("output/energy_calc.dat");
        ofs << w20 << "time" << w20 << "kinetic energy" << w20 << "enstrophy\n";
        ofs << w20 << it * Parameters::dt << w20 << Ek << w20 << Es << "\n";
        ofs.close();
    }
    else if (it >= 1)
    {
        ofs.open("output/energy_calc.dat", std::ofstream::out | std::ofstream::app);
        ofs << w20 << it * Parameters::dt << w20 << Ek << w20 << Es << "\n";
        ofs.close();
    }
} // end of function
