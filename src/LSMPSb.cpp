#include "LSMPSb.hpp"
#include "parameters.hpp"

#pragma region public methods
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void LSMPSb::set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, 
                       const std::vector<double> &s, const std::vector<double> &f,
                       const std::vector<double> &xi, const std::vector<double> &yi, const std::vector<double> &zi, 
                       const std::vector<double> &si, const std::vector<double> &fi,
                       std::vector<std::vector<int>> &neighborlistInter, std::vector<std::vector<int>> &neighborlistIntra)
{
    // clear local variables
    this->_d00.clear();
    this->_ddx.clear();
    this->_ddy.clear();
    this->_ddz.clear();    
    this->_d2d2x.clear();
    this->_d2d2y.clear();
    this->_d2d2z.clear();
    this->_d2dxdy.clear();
    this->_d2dydz.clear();
    this->_d2dzdx.clear();

    // resize local variabels
    int nparticle = x.size();
    this->_d00.resize(nparticle);
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_ddz.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2d2y.resize(nparticle);
    this->_d2d2z.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2dydz.resize(nparticle);
    this->_d2dzdx.resize(nparticle);

    // perform calculation
    this->calculate_LSMPS(x, y, z, s, f, xi, yi, zi, si, fi, neighborlistInter, neighborlistIntra);

}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d00()
{
    return this->_d00;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_ddx()
{
    return this->_ddx;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_ddy()
{
    return this->_ddy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_ddz()
{
    return this->_ddz;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d2d2x()
{
    return this->_d2d2x;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d2d2y()
{
    return this->_d2d2y;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d2d2z()
{
    return this->_d2d2z;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d2dxdy()
{
    return this->_d2dxdy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d2dydz()
{
    return this->_d2dydz;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSb::get_d2dzdx()
{
    return this->_d2dzdx;
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
#pragma endregion

#pragma region private methods
// ===========================================================================
// ===========================================================================
void LSMPSb::calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, 
                             const std::vector<double> &s, const std::vector<double> &f,
                             const std::vector<double> &xi, const std::vector<double> &yi, const std::vector<double> &zi, 
                             const std::vector<double> &si, const std::vector<double> &fi,
                             std::vector<std::vector<int>> &neighborlistInter, std::vector<std::vector<int>> &neighborlistIntra)
{
    int nparticle = x.size();
    std::vector<double> averageDiameter(nparticle);

    for (size_t i = 0; i < nparticle; i++)
    {
        Eigen::MatrixXd Hbar = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::MatrixXd Mbar = Eigen::MatrixXd::Zero(MAT_SIZE, MAT_SIZE);
        Eigen::VectorXd bbar = Eigen::VectorXd::Zero(MAT_SIZE);

        std::vector<int> _neighborInterIndex = neighborlistInter[i];
        std::vector<int> _neighborIntraIndex = neighborlistIntra[i];

        if (_neighborInterIndex.empty())
        {
            // assign to private variables
            this->_d00[i] = 0;
            this->_ddx[i] = 0;
            this->_ddy[i] = 0;
            this->_ddz[i] = 0;            
            this->_d2d2x[i] = 0;
            this->_d2d2y[i] = 0;
            this->_d2d2z[i] = 0;
            this->_d2dxdy[i] = 0;
            this->_d2dydz[i] = 0;
            this->_d2dzdx[i] = 0;
        }
        else
        {
            // TODO: calculate average diameter
            for (size_t j = 0; j < _neighborInterIndex.size(); j++)
            {
                int idx = _neighborInterIndex[j];
                averageDiameter[i] += si[idx];
            }
            // if (_neighborInterIndex.size() > 0)
            // {
                averageDiameter[i] = averageDiameter[i] / (double)_neighborInterIndex.size();
            // }
            // else
            // {
            //     averageDiameter[i] = s[i];
            // }

            Hbar(0, 0) = 1;
            Hbar(1, 1) = std::pow(averageDiameter[i], -1);
            Hbar(2, 2) = std::pow(averageDiameter[i], -1);
            Hbar(3, 3) = std::pow(averageDiameter[i], -1);
            Hbar(4, 4) = std::pow(averageDiameter[i], -2) * 2;
            Hbar(5, 5) = std::pow(averageDiameter[i], -2) * 2;
            Hbar(6, 6) = std::pow(averageDiameter[i], -2) * 2;
            Hbar(7, 7) = std::pow(averageDiameter[i], -2);
            Hbar(8, 8) = std::pow(averageDiameter[i], -2);
            Hbar(9, 9) = std::pow(averageDiameter[i], -2);

            // TODO: perform Least Square MPS
            for (size_t j = 0; j < _neighborInterIndex.size(); j++)
            {
                int idxi = i;
                int idxj = _neighborInterIndex[j];
                double _xi = x[idxi];
                double _yi = y[idxi];
                double _zi = z[idxi];
                double _xj = xi[idxj];
                double _yj = yi[idxj];
                double _zj = zi[idxj];
                double _xij = _xj - _xi;
                double _yij = _yj - _yi;
                double _zij = _zj - _zi;
                double _fij = fi[idxj] - 0.0e0;

                double _rij = std::sqrt(std::pow(_xij/Parameters::er, 2) + std::pow(_yij, 2)+ std::pow(_zij, 2));
                                
                double _Ri = /*R_e*/ 2.1 * averageDiameter[idxi];
                double _Rj = /*R_e*/ 2.1 * si[idxj];
                double _Rij = (_Ri + _Rj) * 0.5;

                std::vector<double> _p1 = this->get_p(_xij, _yij, _zij, averageDiameter[i]);
                std::vector<double> _p2 = this->get_p(_xij, _yij, _zij, averageDiameter[i]);
                double _wij = this->weight_function(_rij, _Rij);

                for (size_t k1 = 0; k1 < MAT_SIZE; k1++)
                {
                    for (size_t k2 = 0; k2 < MAT_SIZE; k2++)
                    {
                        // generate tensor product between p
                        Mbar(k1, k2) = Mbar(k1, k2) + (_wij * _p1[k1] * _p2[k2]);
                    }
                    // generate mooment matrix
                    bbar(k1) = bbar(k1) + (_wij * _p1[k1] * _fij);
                }
            }
            // TODO: to increase robustness <- considering intra-particle neighbor
            if (std::abs(Mbar.determinant()) < 1.0e-6) // ! if tend to be NON-invertible matrix
            {
                for (size_t j = 0; j < _neighborIntraIndex.size(); j++)
                {
                    int idxi = i;
                    int idxj = _neighborIntraIndex[j];
                    double _xi = x[idxi];
                    double _yi = y[idxi];
                    double _zi = z[idxi];
                    double _xj = x[idxj];
                    double _yj = y[idxj];
                    double _zj = z[idxj];
                    double _xij = _xj - _xi;
                    double _yij = _yj - _yi;
                    double _zij = _zj - _zi;
                    double _fij = 0.0e0; // ! no interaction of the function

                    double _rij = std::sqrt(std::pow(_xij/Parameters::er, 2) + std::pow(_yij, 2)+ std::pow(_zij, 2));
                    
                    double _Ri = /*R_e*/ 2.1 * averageDiameter[idxi];
                    double _Rj = /*R_e*/ 2.1 * s[idxj];
                    double _Rij = (_Ri + _Rj) * 0.5;

                    std::vector<double> _p1 = this->get_p(_xij, _yij, _zij, averageDiameter[i]);
                    std::vector<double> _p2 = this->get_p(_xij, _yij, _zij, averageDiameter[i]);
                    double _wij = this->weight_function(_rij, _Rij);

                    for (size_t k1 = 0; k1 < MAT_SIZE; k1++)
                    {
                        for (size_t k2 = 0; k2 < MAT_SIZE; k2++)
                        {
                            // generate tensor product between p
                            Mbar(k1, k2) = Mbar(k1, k2) + (_wij * _p1[k1] * _p2[k2]);
                        }
                        // generate mooment matrix
                        bbar(k1) = bbar(k1) + (_wij * _p1[k1] * _fij);
                    }
                }
            }

            // ? solve Least Square
            // LLT<MatrixXd> llt;
            // llt.compute(Mbar);
            // Eigen::VectorXd MbarInv_Bbar = llt.solve(bbar);
            Eigen::VectorXd MbarInv_Bbar = Mbar.fullPivHouseholderQr().solve(bbar);
            //Eigen::VectorXd MbarInv_Bbar = Mbar.bdcSvd(ComputeThinU | ComputeThinV).solve(bbar); // (MAT_SIZE x 1)
            Eigen::VectorXd Dx = Hbar * MbarInv_Bbar; // (MAT_SIZE x 1)

            // assign to private variables
            this->_d00[i] = Dx[0];
            this->_ddx[i] = Dx[1];
            this->_ddy[i] = Dx[2];
            this->_ddz[i] = Dx[3];
            this->_d2d2x[i] = Dx[4];
            this->_d2d2y[i] = Dx[5];
            this->_d2d2z[i] = Dx[6];
            this->_d2dxdy[i] = Dx[7]; 
            this->_d2dydz[i] = Dx[8];
            this->_d2dzdx[i] = Dx[9];           
        }
    }
}
// ===========================================================================
// ===========================================================================
std::vector<double> LSMPSb::get_p(const double &xij, const double &yij, const double &zij, const double &si)
{
    std::vector<double> _p(MAT_SIZE);

    double _xij = xij / si;
    double _yij = yij / si;
    double _zij = zij / si;

    _p[0] = 1;
    _p[1] = _xij;
    _p[2] = _yij;
    _p[3] = _zij;
    _p[4] = _xij * _xij;
    _p[5] = _yij * _yij;
    _p[6] = _zij * _zij;
    _p[7] = _xij * _yij;
    _p[8] = _yij * _zij;
    _p[9] = _zij * _xij;    
    
    return _p;
}
// ===========================================================================
// ===========================================================================
double LSMPSb::weight_function(const double &rij, const double &Rij)
{
    double _wij;
    if (rij <= Rij)
    {
        _wij = std::pow(1 - (rij / Rij), 2);
        // _wij = std::pow((rij / Rij) - 1, 2);
    }
    else
    {
        _wij = 0.0e0;
    }

    return _wij;
}
// ===========================================================================
// ===========================================================================
#pragma endregion