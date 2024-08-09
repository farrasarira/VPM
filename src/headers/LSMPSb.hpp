#ifndef LSMPSB_H
#define LSMPSB_H

#include <vector>
#include "../Eigen/Dense"

class LSMPSb
{
private:
//// #define MAT_SIZE 6 // size of matrix
#define R_e 3.1                    // effective radius ratio
    static const int MAT_SIZE = 10; // size of matrix

    std::vector<double> _d00;    // d^
    std::vector<double> _ddx;    // d{}/dx
    std::vector<double> _ddy;    // d{}/dy
    std::vector<double> _ddz;    // d{}/dz
    std::vector<double> _d2d2x;  // d^2{}/d^2x
    std::vector<double> _d2d2y;  // d^2{}/d^2y
    std::vector<double> _d2d2z;  // d^2{}/d^2z
    std::vector<double> _d2dxdy; // d^2{}/dxdy
    std::vector<double> _d2dydz; // d^2{}/dydz
    std::vector<double> _d2dzdx; // d^2{}/dzdx

    std::vector<double> get_p(const double &x, const double &y, const double &z, const double &s);
    double weight_function(const double &rij, const double &Rij);

    void calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, 
                         const std::vector<double> &s, const std::vector<double> &f,
                         const std::vector<double> &xi, const std::vector<double> &yi, const std::vector<double> &zi,
                         const std::vector<double> &si, const std::vector<double> &fi,
                         std::vector<std::vector<int>> &neighborlistInter, std::vector<std::vector<int>> &neighborlistIntra);

public:
    void set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, 
                   const std::vector<double> &s, const std::vector<double> &f,
                   const std::vector<double> &xi, const std::vector<double> &yi, const std::vector<double> &zi,
                   const std::vector<double> &si, const std::vector<double> &fi,
                   std::vector<std::vector<int>> &neighborlistInter, std::vector<std::vector<int>> &neighborlistIntra);

    std::vector<double> get_d00();
    std::vector<double> get_ddx();
    std::vector<double> get_ddy();
    std::vector<double> get_ddz();
    std::vector<double> get_d2d2x();
    std::vector<double> get_d2d2y();
    std::vector<double> get_d2d2z();
    std::vector<double> get_d2dxdy();
    std::vector<double> get_d2dydz();
    std::vector<double> get_d2dzdx();

};

#endif