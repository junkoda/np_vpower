//
// Computes Nearest Particle line-of-sight velocity mesh
//

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "particle.h"
#include "binary_tree.h"
#include "nbr_finder.h"
#include "fft_mesh.h"
#include "nearest_nbr_mesh.h"

using namespace std;


void calculate_nearest_particle_u(vector<ParticleData>& v, const float boxsize,
				  const int nc, float * const vmesh)
{
  // Finds nearest neighbor for each grid point, and assign its velocity
  // to the grid

  const int np= v.size();
  const int quota = 8;
  nbr_finder::set_boxsize(boxsize);

  BinaryTree* const tree= new BinaryTree[np];
  index_t ntree= tree->construct(&v.front(), np, quota);
  cerr << ntree << " trees used. searching nearest neighbors...\n";
  KthValue* knbrs= new KthValue(1); // 1= nearest particle

  const int ncz= 2*(nc/2+1);

  float x[3];
  for(int ix=0; ix<nc; ++ix) {
    x[0]= boxsize*(ix + 0.5f)/nc;
    for(int iy=0; iy<nc; ++iy) {
      x[1]= boxsize*(iy + 0.5f)/nc;
      for(int iz=0; iz<nc; ++iz) {
	x[2]= boxsize*(iz + 0.5f)/nc;
	
	index_t i= nbr_finder::for_neighbors_k(tree, x, knbrs);
	if(! (0 <= i && i<np))
	  cerr << "Error: nearest nbr index i=" << i << endl;
	assert(0<=i && i<np);
	
	vmesh[(ix*nc + iy)*ncz + iz]= 0.01f*v[i].v[2]; 
                                            // velocity in units 100 km/s
      }
    }
  }

  double vrms= 0.0;
  for(int ix=0; ix<nc; ++ix) {
    for(int iy=0; iy<nc; ++iy) {
      for(int iz=0; iz<nc; ++iz) {
	float u= vmesh[ncz*(nc*ix + iy) + iz];
	vrms += u*u;
      }
    }
  }
  vrms /= nc*nc*nc;

  printf("# vrms %e\n", sqrt(vrms));
}
