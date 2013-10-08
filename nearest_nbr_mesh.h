#ifndef NEAREST_NBR_MESH_H
#define NEAREST_NBR_MESH_H 1

#include <vector>
#include "particle.h"


void calculate_nearest_particle_u(std::vector<ParticleData>& v, 
				  const float boxsize,
				  const int nc, float vmesh[]);


#endif
