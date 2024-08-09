#include "LSMPSa.hpp"
#include "parameters.hpp"

#pragma region public methods
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void LSMPSa::set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, 
                       const std::vector<double> &s,const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist)
{
    // clear local variables
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
    this->calculate_LSMPS(x, y, z, s, f, neighborlist);
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void LSMPSa::set_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &zSource, 
                       const std::vector<double> &sSource, const std::vector<double> &fSource,
                       const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &zCollocation,
                       const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                       std::vector<std::vector<int>> &neighborlist)
{
    // clear local variables
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
    int nparticle = xSource.size();
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
    this->calculate_LSMPS(xSource, ySource, zSource, sSource, fSource, xCollocation, yCollocation, zCollocation, sCollocation, fCollocation, neighborlist);
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_ddx()
{
    return this->_ddx;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_ddy()
{
    return this->_ddy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_ddz()
{
    return this->_ddz;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2d2x()
{
    return this->_d2d2x;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2d2y()
{
    return this->_d2d2y;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2d2z()
{
    return this->_d2d2z;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2dxdy()
{
    return this->_d2dxdy;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2dydz()
{
    return this->_d2dydz;
}
// ---------------------------------------------------------------------------
std::vector<double> LSMPSa::get_d2dzdx()
{
    return this->_d2dzdx;
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
#pragma endregion

#pragma region private methods
// ===========================================================================
// ===========================================================================
void LSMPSa::calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, 
                             const std::vector<double> &s, const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist)
{
    int nparticle = x.size();
    for (int i = 0; i < nparticle; i++)
    {
        Eigen::MatrixXd Hi = Eigen::MatrixXd::Zero(MAT_SIZE + 10, MAT_SIZE + 10);
        Eigen::MatrixXd Mi = Eigen::MatrixXd::Zero(MAT_SIZE + 10, MAT_SIZE + 10);
        Eigen::VectorXd bi = Eigen::VectorXd::Zero(MAT_SIZE + 10);

        Hi(0, 0) = std::pow(s[i], -1);
        Hi(1, 1) = std::pow(s[i], -1);
        Hi(2, 2) = std::pow(s[i], -1);
        Hi(3, 3) = std::pow(s[i], -2) * 2;
        Hi(4, 4) = std::pow(s[i], -2) * 2;
        Hi(5, 5) = std::pow(s[i], -2) * 2;
        Hi(6, 6) = std::pow(s[i], -2);
        Hi(7, 7) = std::pow(s[i], -2);
        Hi(8, 8) = std::pow(s[i], -2);

        Hi(9, 9) = std::pow(s[i], -3) * 6;
        Hi(10, 10) = std::pow(s[i], -3) * 6;
        Hi(11, 11) = std::pow(s[i], -3) * 6;
        Hi(12, 12) = std::pow(s[i], -3) * 1;
        Hi(13, 13) = std::pow(s[i], -3) * 2;
        Hi(14, 14) = std::pow(s[i], -3) * 2;
        Hi(15, 15) = std::pow(s[i], -3) * 2;
        Hi(16, 16) = std::pow(s[i], -3) * 2;
        Hi(17, 17) = std::pow(s[i], -3) * 2;
        Hi(18, 18) = std::pow(s[i], -3) * 2;

        std::vector<int> _neighborIndex = neighborlist[i];
        for (size_t j = 0; j < _neighborIndex.size(); j++)
        {
            int idxi = i;
            int idxj = _neighborIndex[j];
            double _xi = x[idxi];
            double _yi = y[idxi];
            double _zi = z[idxi];
            double _xj = x[idxj];
            double _yj = y[idxj];
            double _zj = z[idxj];
            double _xij = _xj - _xi;
            double _yij = _yj - _yi;
            double _zij = _zj - _zi;
            double _fij = f[idxj] - f[idxi];
            
            double _rij = std::sqrt(std::pow(_xij/Parameters::er, 2) + std::pow(_yij, 2) + std::pow(_zij, 2));
            
            double _Ri = R_e * s[idxi];
            double _Rj = R_e * s[idxj];
            double _Rij = (_Ri + _Rj) * 0.5;

            std::vector<double> _p1 = this->get_p(_xij, _yij, _zij, s[idxi]);
            std::vector<double> _p2 = this->get_p(_xij, _yij, _zij, s[idxi]);
            double _wij = this->weight_function(_rij, _Rij);

            for (size_t k1 = 0; k1 < MAT_SIZE + 10; k1++)
            {
                for (size_t k2 = 0; k2 < MAT_SIZE + 10; k2++)
                {
                    // generate tensor product between p
                    Mi(k1, k2) = Mi(k1, k2) + (_wij * _p1[k1] * _p2[k2]);
                }
                // generate mooment matrix
                bi(k1) = bi(k1) + (_wij * _p1[k1] * _fij);
            }
        }

        // solve Least Square
        Eigen::VectorXd MiInv_Bi = Mi.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bi); // (MAT_SIZE x 1)
        Eigen::VectorXd Dx = Hi * MiInv_Bi;                                          // (MAT_SIZE x 1)

        // assign to private variables
        this->_ddx[i] = Dx[0];
        this->_ddy[i] = Dx[1];
        this->_ddz[i] = Dx[2];
        this->_d2d2x[i] = Dx[3];
        this->_d2d2y[i] = Dx[4];
        this->_d2d2z[i] = Dx[5];
        this->_d2dxdy[i] = Dx[6];
        this->_d2dydz[i] = Dx[7];
        this->_d2dzdx[i] = Dx[8];
    }
}
// ===========================================================================
// ===========================================================================
std::vector<double> LSMPSa::get_p(const double &xij, const double &yij, const double &zij, const double &si)
{
    std::vector<double> _p(MAT_SIZE + 10);

    double _xij = xij / si;
    double _yij = yij / si;
    double _zij = zij / si;

    _p[0] = _xij;
    _p[1] = _yij;
    _p[2] = _zij;    
    _p[3] = _xij * _xij;
    _p[4] = _yij * _yij;
    _p[5] = _zij * _zij;
    _p[6] = _xij * _yij;
    _p[7] = _yij * _zij;
    _p[8] = _zij * _xij;
    
    _p[9] = _xij * _xij * _xij;
    _p[10] = _yij * _yij * _yij;
    _p[11] = _zij * _zij * _zij;
    _p[12] = _xij * _yij * _zij;
    _p[13] = _xij * _xij * _yij;
    _p[14] = _xij * _xij * _zij;
    _p[15] = _yij * _yij * _xij;
    _p[16] = _yij * _yij * _zij;
    _p[17] = _zij * _zij * _xij;
    _p[18] = _zij * _zij * _yij;
    return _p;
}
// ===========================================================================
// ===========================================================================
double LSMPSa::weight_function(const double &rij, const double &Rij)
{
    double _wij;
    if (rij <= Rij)
    {
        _wij = std::pow(1 - (rij / Rij), 2);
    }
    else
    {
        _wij = 0.0e0;
    }

    return _wij;
}
// ===========================================================================
// ===========================================================================
void LSMPSa::calculate_LSMPS(const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &zSource, 
                             const std::vector<double> &sSource, const std::vector<double> &fSource,
                             const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &zCollocation, 
                             const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                             std::vector<std::vector<int>> &neighborlist)
{
    int nparticle = xSource.size();
    for (int i = 0; i < nparticle; i++)
    {
        Eigen::MatrixXd Hi = Eigen::MatrixXd::Zero(MAT_SIZE + 10, MAT_SIZE + 10);
        Eigen::MatrixXd Mi = Eigen::MatrixXd::Zero(MAT_SIZE + 10, MAT_SIZE + 10);
        Eigen::VectorXd bi = Eigen::VectorXd::Zero(MAT_SIZE + 10);

        Hi(0, 0) = std::pow(sSource[i], -1);
        Hi(1, 1) = std::pow(sSource[i], -1);
        Hi(2, 2) = std::pow(sSource[i], -1);
        Hi(3, 3) = std::pow(sSource[i], -2) * 2;
        Hi(4, 4) = std::pow(sSource[i], -2) * 2;
        Hi(5, 5) = std::pow(sSource[i], -2) * 2;
        Hi(6, 6) = std::pow(sSource[i], -2);
        Hi(7, 7) = std::pow(sSource[i], -2);
        Hi(8, 8) = std::pow(sSource[i], -2);

        Hi(9, 9) = std::pow(sSource[i], -3) * 6;
        Hi(10, 10) = std::pow(sSource[i], -3) * 6;
        Hi(11, 11) = std::pow(sSource[i], -3) * 6;
        Hi(12, 12) = std::pow(sSource[i], -3) * 1;
        Hi(13, 13) = std::pow(sSource[i], -3) * 2;
        Hi(14, 14) = std::pow(sSource[i], -3) * 2;
        Hi(15, 15) = std::pow(sSource[i], -3) * 2;
        Hi(16, 16) = std::pow(sSource[i], -3) * 2;
        Hi(17, 17) = std::pow(sSource[i], -3) * 2;
        Hi(18, 18) = std::pow(sSource[i], -3) * 2;

        std::vector<int> _neighborIndex = neighborlist[i];
        for (size_t j = 0; j < _neighborIndex.size(); j++)
        {
            int idxi = i;
            int idxj = _neighborIndex[j];
            double _xi = xSource[idxi];
            double _yi = ySource[idxi];
            double _zi = zSource[idxi];
            double _xj = xCollocation[idxj];
            double _yj = yCollocation[idxj];
            double _zj = zCollocation[idxj];
            double _xij = _xj - _xi;
            double _yij = _yj - _yi;
            double _zij = _zj - _zi;
            double _fij = fCollocation[idxj] - fSource[idxi];

            double _rij = std::sqrt(std::pow(_xij/Parameters::er, 2) + std::pow(_yij, 2)+ std::pow(_zij, 2));
            
            double _Ri = R_e * sSource[idxi];
            double _Rj = R_e * sCollocation[idxj];
            double _Rij = (_Ri + _Rj) * 0.5;

            std::vector<double> _p1 = this->get_p(_xij, _yij, _zij, sSource[idxi]);
            std::vector<double> _p2 = this->get_p(_xij, _yij, _zij, sSource[idxi]);
            double _wij = this->weight_function(_rij, _Rij);

            for (size_t k1 = 0; k1 < MAT_SIZE+10; k1++)
            {
                for (size_t k2 = 0; k2 < MAT_SIZE +10; k2++)
                {
                    // generate tensor product between p
                    Mi(k1, k2) = Mi(k1, k2) + (_wij * _p1[k1] * _p2[k2]);
                }
                // generate mooment matrix
                bi(k1) = bi(k1) + (_wij * _p1[k1] * _fij);
            }
        }

        // solve Least Square
        Eigen::VectorXd MiInv_Bi = Mi.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bi); // (MAT_SIZE x 1)
        Eigen::VectorXd Dx = Hi * MiInv_Bi;                                          // (MAT_SIZE x 1)

        // assign to private variables
        this->_ddx[i] = Dx[0];
        this->_ddy[i] = Dx[1];
        this->_ddz[i] = Dx[2];
        this->_d2d2x[i] = Dx[3];
        this->_d2d2y[i] = Dx[4];
        this->_d2d2z[i] = Dx[5];
        this->_d2dxdy[i] = Dx[6];
        this->_d2dydz[i] = Dx[7];
        this->_d2dzdx[i] = Dx[8];

    }
}
// ===========================================================================
// ===========================================================================
#pragma endregion