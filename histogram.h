#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1

#include <cmath>

class Histogram {
 public:
  Histogram(const float xmin, const float xmax, const float dx);
  ~Histogram();
  void fill(const float x_, const float y_) {
    int i= (int) floor((x_-xmin)*dxinv);
    if(0 <= i && i < nbin) {
      an[i]++;
      ax[i] += x_;
      ay[i] += y_;
      count++;
    }
  }
  int numbin() const { return nbin; }
  int   n(const int i) const { return an[i]; }
  float pdf(const int i) const { return (double)an[i]/count*dxinv; } 
  float x(const int i) const { return ax[i]/an[i]; }
  float y(const int i) const { return ay[i]/an[i]; }
  float x_left(const int i) const { return xmin + i/dxinv; }
  float x_right(const int i) const { return xmin + (i+1)/dxinv; }
  float x_mid(const int i) const { return xmin + (i+0.5)/dxinv; }
 private:
  const float xmin, xmax, dxinv;
  int nbin;
  unsigned long long count;
  int* an;
  double *ax, *ay;
};

#endif
