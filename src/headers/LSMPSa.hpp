#ifndef LSMPSA_H
#define LSMPSA_H

#include <vector>
#include "../Eigen/Dense"

//// template <typename T>
class LSMPSa
{
private:
//// #define MAT_SIZE 5 // size of matrix
#define R_e 3.1                    // effective radius ratio
    static const int MAT_SIZE = 9; // size of matrix

    std::vector<double> _ddx;    // d{}/dx
    std::vector<double> _ddy;    // d{}/dy
    std::vector<double> _ddz;    // d{}/dy
    std::vector<double> _d2d2x;  // d^2{}/d^2x
    std::vector<double> _d2d2y;  // d^2{}/d^2y
    std::vector<double> _d2d2z;  // d^2{}/d^2y
    std::vector<double> _d2dxdy; // d^2{}/dxdy
    std::vector<double> _d2dydz; // d^2{}/dydz
    std::vector<double> _d2dzdx; // d^2{}/dzdx

    std::vector<double> get_p(const double &x, const double &y, const double &z, const double &s);
    double weight_function(const double &rij, const double &Rij);

    void calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, const std::vector<double> &s,
                         const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist);
    void calculate_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &zSource, 
                         const std::vector<double> &sSource, const std::vector<double> &fSource,
                         const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &zCollocation,
                         const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                         std::vector<std::vector<int>> &neighborlist);

public:
    // LSMPS(/* args */);
    // ~LSMPS();

    void set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, const std::vector<double> &s,
                   const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist);
    void set_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &zSource, 
                   const std::vector<double> &sSource, const std::vector<double> &fSource,
                   const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &zCollocation,
                   const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                   std::vector<std::vector<int>> &neighborlist);

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