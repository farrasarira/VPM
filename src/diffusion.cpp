#include "diffusion.hpp"
#include "parameters.hpp"
#include <time.h>
#include <stdio.h>
#include <algorithm>
#include "LSMPSa.hpp"

void diffusion::main_diffusion(Particle &p, std::vector<std::vector<double>> &dfdt)
{
    // * internal variables
    const static int neighbor_scale = Parameters::r_scale; // ! neighbor_scale. use from global parameters [DONE]
    Particle _particle;                              // active particle and its buffer region
    std::vector<int> _index;                         // to store 'active particle + buffer zone' index

    clock_t t; // defines clock
    printf("\nCalculating diffusion...");
    t = clock();

    // TODO: assign active particles
    for (size_t i = 0; i < p.num; i++)
    {
        if (p.isActive[i] == true)
        {
            std::vector<int> _neighbor = p.neighbor[i];
            for (size_t j = 0; j < _neighbor.size(); j++)
            {
                _index.push_back(_neighbor[j]); // assign its neighbor
            }
            _index.push_back(i); // assign its own value
        }
    }
    // store only unique ID
    std::sort(_index.begin(), _index.end());
    std::vector<int>::iterator it;
    it = std::unique(_index.begin(), _index.end());
    _index.resize(std::distance(_index.begin(), it));

    for (it = _index.begin(); it != _index.end(); it++)
    {
        _particle.x.push_back(p.x[*it]);
        _particle.y.push_back(p.y[*it]);
        _particle.z.push_back(p.z[*it]);
        _particle.s.push_back(p.s[*it]);
        _particle.vort_x.push_back(p.vort_x[*it]);
        _particle.vort_y.push_back(p.vort_y[*it]);
        _particle.vort_z.push_back(p.vort_z[*it]);
        _particle.gx.push_back(p.gx[*it]);
        _particle.gy.push_back(p.gy[*it]);
        _particle.gz.push_back(p.gz[*it]);
        _particle.u.push_back(p.u[*it]);
        _particle.v.push_back(p.v[*it]);
        _particle.w.push_back(p.w[*it]);
        _particle.neighbor.push_back(p.neighbor[*it]);
    }
    _particle.num = _index.size();

    // TODO : Diffusion eq.
    // ! USE LSMPS type A - diffusion
    // lsmpsa.set_LSMPS(p.x, p.y, p.s, p.gz, p.neighbor);
    lsmpsa_x.set_LSMPS(_particle.x, _particle.y, _particle.z, _particle.s, _particle.vort_x, p.x, p.y, p.z, p.s, p.vort_x, _particle.neighbor);
    lsmpsa_y.set_LSMPS(_particle.x, _particle.y, _particle.z, _particle.s, _particle.vort_y, p.x, p.y, p.z, p.s, p.vort_y, _particle.neighbor);
    lsmpsa_z.set_LSMPS(_particle.x, _particle.y, _particle.z, _particle.s, _particle.vort_z, p.x, p.y, p.z, p.s, p.vort_z, _particle.neighbor);
    // x - dir
    std::vector<double> _d2fxd2x = lsmpsa_x.get_d2d2x();
    std::vector<double> _d2fxd2y = lsmpsa_x.get_d2d2y();
    std::vector<double> _d2fxd2z = lsmpsa_x.get_d2d2z();
    // y - dir
    std::vector<double> _d2fyd2x = lsmpsa_y.get_d2d2x();
    std::vector<double> _d2fyd2y = lsmpsa_y.get_d2d2y();
    std::vector<double> _d2fyd2z = lsmpsa_y.get_d2d2z();
    // z - dir
    std::vector<double> _d2fzd2x = lsmpsa_z.get_d2d2x();
    std::vector<double> _d2fzd2y = lsmpsa_z.get_d2d2y();
    std::vector<double> _d2fzd2z = lsmpsa_z.get_d2d2z();

    dfdt.resize(3,std::vector<double>(p.num, 0.0e0)); // output
    for (int i = 0; i < _particle.num; i++)
    {
        dfdt[0][_index[i]] = Parameters::vis * (_d2fxd2x[i] + _d2fxd2y[i] + _d2fxd2z[i]);
        dfdt[1][_index[i]] = Parameters::vis * (_d2fyd2x[i] + _d2fyd2y[i] + _d2fyd2z[i]);
        dfdt[2][_index[i]] = Parameters::vis * (_d2fzd2x[i] + _d2fzd2y[i] + _d2fzd2z[i]);
    }
    t = clock() - t;
    printf(" [%f s].\n", ((float)t) / CLOCKS_PER_SEC);

    _d2fxd2x.clear();
    _d2fxd2y.clear();
    _d2fxd2z.clear();
    _d2fyd2x.clear();
    _d2fyd2y.clear();
    _d2fyd2z.clear();
    _d2fzd2x.clear();
    _d2fzd2y.clear();
    _d2fzd2z.clear();
}
