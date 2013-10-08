#ifndef HALO_FILE_H
#define HALO_FILE_H 1

#include <vector>
#include "particle.h"

void read_mock_text(const char filename[], std::vector<ParticleData>& v);
void read_fof_text(const char filename[], std::vector<ParticleData>& v, 
		   const float m, const float logMmin, const float logmMax);
void read_fof_binary(const char filename[], std::vector<ParticleData>& v, const int nfof_min);

#endif
