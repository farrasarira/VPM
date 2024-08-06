#ifndef INIT_H
#define INIT_H

#include "particle.hpp"

class initialization
{
private:
    base_grid d_base_grid;
    neighbor d_neighbor;
    /* data */
public:
    initialization(/* args */);
    ~initialization();
    //2D Particles initialization
    void init_particle(const Body &b, Particle &p);
    void init_particle_test(Particle &p);
    void init_domain(Particle &p);

    //3D Particles initialization
    void init_domain_3d(Particle &p);
};

#endif