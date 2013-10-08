//
// A utility of histogram
//
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "histogram.h"

using namespace std;

Histogram::Histogram(const float xmin_, const float xmax_, const float dx) :
  xmin(xmin_), xmax(xmax_), dxinv(1.0f/dx), count(0)
{
  nbin= (int) floor((xmax - xmin)/dx);
  an= (int*) calloc(sizeof(int), nbin); assert(an);
  ax= (double*) calloc(sizeof(double), nbin); assert(ax);
  ay= (double*) calloc(sizeof(double), nbin); assert(ay);
}

Histogram::~Histogram()
{
  free(an); free(ax); free(ay);
}
