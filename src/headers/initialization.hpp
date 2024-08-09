#ifndef INIT_H
#define INIT_H

#include "particle.hpp"
#include "neighbor.hpp"

class initialization
{
private:
    // base_grid d_base_grid;
    neighbor d_neighbor;
    /* data */
public:
    initialization(/* args */);
    ~initialization();

    //3D Particles initialization
    void init_domain(Particle &p);
    void init_vortex_ring(Particle &p);


};

#endif