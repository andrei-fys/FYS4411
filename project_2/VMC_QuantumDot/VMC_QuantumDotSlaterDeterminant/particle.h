#ifndef PARTICLE_H
#define PARTICLE_H
#include "vec3.h"

class Particle
{
public:
    Particle();
    vec3 position;
    vec3 positionNew;
    double spin;
};

#endif // PARTICLE_H
