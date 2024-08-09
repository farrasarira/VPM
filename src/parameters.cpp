#include "parameters.hpp"
#include <math.h>

namespace Parameters {

    extern const double pi = M_PI;
    extern const double lxdom = 5.0;
    extern const double lydom = 5.0;
    extern const double lzdom = 5.0;
    extern const double xdom = 1.5;

    extern const double Re = 250.0e0; // REYNOLDs NUMBER
    extern const double vis = 1.0e0 / Re;    
    extern const double sigma = 0.1;                                  // Coresize, default: 0.0025, dt~(1/2)*sigma - (1/3)*sigma


    /** @brief parallel_comp */
    extern const double simulation_time = 10.0e0;
    extern const double comtime_sf = 1200.0e0 * 6000.0e0 * 6000.0e0; // % Save file frequently after how many [second]=[hour]*60*60, only in case of running is stopped suddenly, but nt_sf take long days
    extern const double dt = 0.0005;      // default:0.001, dt <= phi_s*sigma^2/vis (Ploumhans [2000]) OR dt = phi_s*sigma^2*Courant^2/vis, where 0 < Courant <= 1
    extern const int nt = std::ceil(simulation_time / dt); // number of time step
    extern const int nt_sf = 50;                           // save data per ( nt_sf ) time step, for saving storage MEMORY and for case of running is stopped suddenly

    extern const double ly = 1.0;              // 3D Unbounded : vortex ring radius, 3D flat plate: length, 3D cylinder disc: diameter
    extern const double lx = 1.0;              // 3D Unbounded : vortex ring radius, 3D cylinder disc: thickness
    extern const double lz = 1.0; // 3D Unbounded : vortex ring thickness, 3D cylinder disc: diameter

    extern const int r_scale = 3 * 1.0; 
    extern const double er = 1.0;       // elliptical ratio ex:ey


    extern const int icutoff = 2;    // 0.singular ;  1.super (high-oder) algebraic  ; 2. Gaussian  ; 3 super Gaussian


};