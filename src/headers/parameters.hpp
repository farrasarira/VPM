#ifndef PARAMETERS_H
#define PARAMETERS_H



/** @brief Defines certain parameters used in the VPM */

namespace Parameters {

    extern const double pi;
    extern const double lxdom;
    extern const double lydom;
    extern const double lzdom;
    extern const double xdom;

    extern const double sigma;   // Coresize, default: 0.0025, dt~(1/2)*sigma - (1/3)*sigma

    extern const double Re; // REYNOLDs NUMBER
    extern const double vis;    

    /** @brief parallel_comp */
    extern const double simulation_time;
    extern const double comtime_sf; // % Save file frequently after how many [second]=[hour]*60*60, only in case of running is stopped suddenly, but nt_sf take long days
    extern const double dt;         // default:0.001, dt <= phi_s*sigma^2/vis (Ploumhans [2000]) OR dt = phi_s*sigma^2*Courant^2/vis, where 0 < Courant <= 1
    extern const int nt;            // number of time step
    extern const int nt_sf;         // save data per ( nt_sf ) time step, for saving storage MEMORY and for case of running is stopped suddenly

    extern const double ly;         // 3D Unbounded : vortex ring radius, 3D flat plate: length, 3D cylinder disc: diameter
    extern const double lx;         // 3D Unbounded : vortex ring radius, 3D cylinder disc: thickness
    extern const double lz;         // 3D Unbounded : vortex ring thickness, 3D cylinder disc: diameter

    extern const int r_scale; 
    extern const double er;       // elliptical ratio ex:ey


    extern const int icutoff;    // 0.singular ;  1.super (high-oder) algebraic  ; 2. Gaussian  ; 3 super Gaussian


};

#endif