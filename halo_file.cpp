//
// Reads halo data from file
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "halo_file.h"

using namespace std;

void read_mock_text(const char filename[], vector<ParticleData>& v)
{
  // Ascii text file
  // Column 1-3: x[3]  halo positions  [1/h Mpc]
  // Column 4-6: v[3]  halo velocities [km/s]

  ParticleData p;

  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open fof file: " << filename << endl;
    throw 1;
  }

  char buf[256];
  while(fgets(buf, 255, fp)) {
    int ret= sscanf(buf, "%e %e %e %e %e %e", 
	     p.x, p.x+1, p.x+2, p.v, p.v+1, p.v+2);
    if(ret == 6)
      v.push_back(p);
  }
}

void read_fof_text(const char filename[], vector<ParticleData>& v, 
		   const float m, const float logMmin, const float logMmax)
{
  // Ascii text file with mass information
  // Column 1:   Number of halo particles
  // Column 2-4: x[3]  halo positions  [1/h Mpc]
  // Column 5-7: v[3]  halo velocities [km/s]

  ParticleData p;
  const float M_min= pow(10.0f, logMmin);
  const float M_max= pow(10.0f, logMmax);

  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open fof file: " << filename << endl;
    throw 1;
  }

  int nfof;
  char buf[256];
  while(fgets(buf, 255, fp)) {
    int ret= sscanf(buf, "%d %e %e %e %e %e %e", 
	     &nfof, p.x, p.x+1, p.x+2, p.v, p.v+1, p.v+2);
    if(ret == 7 && M_min <= m*nfof && m*nfof < M_max)
      v.push_back(p);
  }
}

void read_fof_binary(const char filename[], vector<ParticleData>& v, const int nfof_min)
{
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open fof file: " << filename << endl;
    throw "FOF File Error";
  }

  int id;
  float ll;
  int nhalo;
  float boxsize;

  fread(&id, sizeof(int), 1, fp);
  assert(id == 604);
  fread(&ll, sizeof(float), 1, fp);
  fread(&boxsize, sizeof(float), 1, fp);
  fread(&nhalo, sizeof(int), 1, fp); assert(nhalo);

  int nfof;
  ParticleData p;

  for(int j=0; j<nhalo; ++j) {
    fread(&nfof, sizeof(int), 1, fp);
    fread(p.x, sizeof(float), 3, fp);
    fread(p.v, sizeof(float), 3, fp);
    if(nfof >= nfof_min)
      v.push_back(p);
  }
  
  int nhalo_check= 0;
  int ret= fread(&nhalo_check, sizeof(int), 1, fp); assert(ret == 1);
  assert(nhalo_check == nhalo);

  fclose(fp);
}
