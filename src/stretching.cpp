#include "stretching.hpp"
#include "parameters.hpp"
#include <time.h>
#include <stdio.h>
#include <algorithm>

void stretching::main_stretching(Particle &p, std::vector<std::vector<double>> &dfdt)
{
    // * internal variables
    const static int neighbor_scale = Parameters::r_scale;  // ! neighbor_scale. use from global parameters [DONE]
    Particle _particle;                                     // active particle and its buffer region
    std::vector<int> _index;                                // to store 'active particle + buffer zone' index

    clock_t t; // defines clock
    printf("\nCalculating stretching...");
    t = clock();

    // TODO: assign active particles
    for (int i = 0; i < p.num; i++)
    {
        if (p.isActive[i] == true)
        {
            std::vector<int> _neighbor = p.neighbor[i];
            for (int j = 0; j < _neighbor.size(); j++)
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
        _particle.x.emplace_back(p.x[*it]);
        _particle.y.emplace_back(p.y[*it]);
        _particle.z.emplace_back(p.z[*it]);
        _particle.s.emplace_back(p.s[*it]);
        _particle.vort_x.emplace_back(p.vort_x[*it]);
        _particle.vort_y.emplace_back(p.vort_y[*it]);
        _particle.vort_z.emplace_back(p.vort_z[*it]);
        _particle.gx.emplace_back(p.gx[*it]);
        _particle.gy.emplace_back(p.gy[*it]);
        _particle.gz.emplace_back(p.gz[*it]);
        _particle.u.emplace_back(p.u[*it]);
        _particle.v.emplace_back(p.v[*it]);
        _particle.w.emplace_back(p.w[*it]);
        _particle.neighbor.emplace_back(p.neighbor[*it]);
    }
    _particle.num = _index.size();

    // TODO : Stretching term
    // ! USE LSMPS type A - diffusion
    // lsmpsa.set_LSMPS(p.x, p.y, p.s, p.gz, p.neighbor);
    lsmpsa_x.set_LSMPS(_particle.x, _particle.y, _particle.z, _particle.s, _particle.u, p.x, p.y, p.z, p.s, p.u, _particle.neighbor);
    lsmpsa_y.set_LSMPS(_particle.x, _particle.y, _particle.z, _particle.s, _particle.v, p.x, p.y, p.z, p.s, p.v, _particle.neighbor);
    lsmpsa_z.set_LSMPS(_particle.x, _particle.y, _particle.z, _particle.s, _particle.w, p.x, p.y, p.z, p.s, p.w, _particle.neighbor);
    // x - dir
    std::vector<double> _duxdx = lsmpsa_x.get_ddx();
    std::vector<double> _duxdy = lsmpsa_x.get_ddy();
    std::vector<double> _duxdz = lsmpsa_x.get_ddz();
    // y - dir
    std::vector<double> _duydx = lsmpsa_y.get_ddx();
    std::vector<double> _duydy = lsmpsa_y.get_ddy();
    std::vector<double> _duydz = lsmpsa_y.get_ddz();
    // z - dir
    std::vector<double> _duzdx = lsmpsa_z.get_ddx();
    std::vector<double> _duzdy = lsmpsa_z.get_ddy();
    std::vector<double> _duzdz = lsmpsa_z.get_ddz();

    dfdt.resize(3,std::vector<double>(p.num, 0.0e0)); // output
    for (int i = 0; i < _particle.num; i++)
    {
        dfdt[0][_index[i]] = (p.vort_x[_index[i]] * _duxdx[i]) + (p.vort_y[_index[i]] * _duxdy[i]) + (p.vort_z[_index[i]] * _duxdz[i]);
        dfdt[1][_index[i]] = (p.vort_x[_index[i]] * _duydx[i]) + (p.vort_y[_index[i]] * _duydy[i]) + (p.vort_z[_index[i]] * _duydz[i]);
        dfdt[2][_index[i]] = (p.vort_x[_index[i]] * _duzdx[i]) + (p.vort_y[_index[i]] * _duzdy[i]) + (p.vort_z[_index[i]] * _duzdz[i]);
    }
    t = clock() - t;
    printf(" [%f s].\n", ((float)t) / CLOCKS_PER_SEC);

    _duxdx.clear();
    _duxdy.clear();
    _duxdz.clear();
    _duydx.clear();
    _duydy.clear();
    _duydz.clear();
    _duzdx.clear();
    _duzdy.clear();
    _duzdz.clear();        
}
