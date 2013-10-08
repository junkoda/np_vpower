//
// Calculates angle-averaged power spectra from given Fourier modes,
// delta_k and u_k
// Shotnoise and alias correction performed assuming,
// CIC density + nearest neighbor velocity

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "histogram.h"
#include "power_spectra.h"

using namespace std;

struct FourierMesh {
  float *delta_k, *u_k;
  int nc, np;
  float boxsize;
};

static inline float w(const float tx)
{
  return tx == 0.0f ? 1.0f : sin(tx)/tx;
}

void calc_power_spectrum_sa(const int nc, const float boxsize, const float delta_k[], const float u_k[], const float nbar, const float neff, const float dk, float kmax)
{
  // Power spectrum subtracting shot noise and aliasing correction
  // P(k) ~ k^neff
  const float knq= (M_PI*nc)/boxsize;
  const float knq_inv= boxsize/(M_PI*nc);
  const float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  const float sin_fac= 0.5*boxsize/nc;

  const float fac= 2.0*M_PI/boxsize;
  const int ncz= nc/2+1;

  const float kp= 2.0*M_PI*nc/boxsize;
  const float nbar_inv= nbar > 0.0f ? 1.0/nbar : 0.0f;

  const int na=2;

  printf("# dk kmax %e %e\n", dk, kmax);
  printf("# k n-mode Pgg Pgu Puu Pgg_raw Pgu_raw Puu_raw\n");

  if(kmax == 0.0f) kmax= nc*M_PI/boxsize; //Nyquist frequency

  Histogram Pgg(0.0, kmax, dk);
  Histogram Pgu(0.0, kmax, dk);
  Histogram Puu(0.0, kmax, dk);

  Histogram Pgg_raw(0.0, kmax, dk);
  Histogram Pgu_raw(0.0, kmax, dk);
  Histogram Puu_raw(0.0, kmax, dk);

  for(int ix=0; ix<nc; ++ix) {
    float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    float sintx= sin(sin_fac*kx);

    float c1x= 1.0f - 2.0f/3.0f*sintx*sintx;

    for(int iy=0; iy<nc; ++iy) {
      float ky= iy <= nc/2 ? fac*iy : fac*(iy-nc);

      float sinty= sin(sin_fac*ky);
      float c1y= 1.0f - 2.0f/3.0f*sinty*sinty;

      // Avoid double counting in the kz=0 plane
      // delta_[k] and delta[-k] have same information by the reality condition
      int iz0 = !(kx > 0.0f || (kx == 0.0f && ky > 0.0f));

      for(int iz=iz0; iz<nc/2+1; ++iz) {
	float kz= fac*iz;
	float sintz= sin(sin_fac*kz);
	float c1z= 1.0f - 2.0f/3.0f*sintz*sintz;

	float k= sqrt(kx*kx + ky*ky + kz*kz);
	float shot_noise= c1x*c1y*c1z*nbar_inv; // C1 function in Jing 2005

	if(k <= knq) {

	// C2 function in Jing 2005
	float c2gg= 0.0f, c2gu= 0.0f;
	for(int ax=-na; ax<na; ++ax) {
	  float kax= kx+ax*kp;
	  for(int ay=-na; ay<na; ++ay) {
	    float kay= ky+ay*kp;
	    for(int az=-na; az<na; ++az) {
	      float kaz= kz+az*kp;
	      float ka= sqrt(kax*kax + kay*kay + kaz*kaz);
	      //float kk= ka/k;

	      float w1= w(sin_fac*(kx+ax*kp))*
		        w(sin_fac*(ky+ay*kp))*
                        w(sin_fac*(kz+az*kp));
	      float w2= w1*w1;
	      float w4= w2*w2;
	      float Pkfac= pow(ka/k, neff);;
	      c2gg += w4*Pkfac;
	      //c2gu += w4*w1*Pkfac/kk; // debug***. empirical w5
	      //c2uu += Pkfac/(kk*kk);
	      //fprintf(stderr, "%d %d %d %e\n", ax, ay, az, Pkfac);
	    }
	  }
	}

	// Empirical c2uu factor (understood partially for the value 45)
	float c2uu = 1.0 + 45.0f*pow(0.5f*k*knq_inv, 2.0f-neff);

	// empirical c2gu factor
	{
	  float w1= w(sin_fac*k);
	  float w2= w1*w1;
	  float w4= w2*w2;
	  c2gu = w1*w4*(1.0 - 0.27*pow(k*knq_inv, 5.0f));
	}

	assert(c2gg > 0.0f && c2gu > 0.0f && c2uu > 0.0f);

	int index= ncz*(nc*ix + iy) + iz;
	float d_re= delta_k[2*index  ];
	float d_im= delta_k[2*index+1];

	float u_re= u_k[2*index  ];
	float u_im= u_k[2*index+1];

	Pgg.fill(k, (pk_fac*(d_re*d_re + d_im*d_im) - shot_noise)/c2gg);
	Pgu.fill(k,  pk_fac*(d_re*u_im - d_im*u_re)/c2gu);
	Puu.fill(k,  pk_fac*(u_re*u_re + u_im*u_im)/c2uu);

	Pgg_raw.fill(k, pk_fac*(d_re*d_re + d_im*d_im));
	Pgu_raw.fill(k, pk_fac*(d_re*u_im - d_im*u_re));
	Puu_raw.fill(k, pk_fac*(u_re*u_re + u_im*u_im));

	}
      }
    }
  }

  printf("# Column 1: k\n");
  printf("# Column 2: nmode\n");
  printf("# Column 3: Pgg(k)\n");
  printf("# Column 4: Pgu(k)\n");
  printf("# Column 5: Puu(k)\n");
  printf("# Column 6: Pgg(k) without shotnoise/alias correction\n");
  printf("# Column 7: Pgu(k) without alias correction\n");
  printf("# Column 8: Puu(k) without alias correction\n")

  for(int j=0; j<Pgg.numbin(); ++j) {
    if(Pgg.n(j) > 0)
      printf("%e %d %e %e %e %e %e %e\n", 
	     Pgg.x(j), Pgg.n(j), Pgg.y(j), Pgu.y(j), Puu.y(j),
	     Pgg_raw.y(j), Pgu_raw.y(j), Puu_raw.y(j)
	     );
  }
}


void calc_momentum_power_sa(const int nc, const float boxsize, const float delta_k[], const float u_k[], const float nbar, const float neff, const float dk, float kmax)
{
  // Power spectrum subtracting shot noise and aliasing correction
  // P(k) ~ k^neff
  const float knq= (M_PI*nc)/boxsize;
  const float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  const float sin_fac= 0.5*boxsize/nc;

  const float fac= 2.0*M_PI/boxsize;
  const int ncz= nc/2+1;

  const float kp= 2.0*M_PI*nc/boxsize;
  const float nbar_inv= nbar > 0.0f ? 1.0/nbar : 0.0f;

  const int na=2;

  printf("# dk kmax %e %e\n", dk, kmax);
  printf("# k n-mode Pgg Pgu Puu Pgg_raw Pgu_raw Puu_raw\n");

  if(kmax == 0.0f) kmax= nc*M_PI/boxsize; //Nyquist frequency

  Histogram Pgg(0.0, kmax, dk);
  Histogram Pgu(0.0, kmax, dk);
  Histogram Puu(0.0, kmax, dk);

  Histogram Pgg_raw(0.0, kmax, dk);
  Histogram Pgu_raw(0.0, kmax, dk);
  Histogram Puu_raw(0.0, kmax, dk);

  for(int ix=0; ix<nc; ++ix) {
    float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    float sintx= sin(sin_fac*kx);

    float c1x= 1.0f - 2.0f/3.0f*sintx*sintx;

    for(int iy=0; iy<nc; ++iy) {
      float ky= iy <= nc/2 ? fac*iy : fac*(iy-nc);

      float sinty= sin(sin_fac*ky);
      float c1y= 1.0f - 2.0f/3.0f*sinty*sinty;

      int iz0 = !(kx > 0.0f || (kx == 0.0f && ky > 0.0f));

      for(int iz=iz0; iz<nc/2+1; ++iz) {
	float kz= fac*iz;
	float sintz= sin(sin_fac*kz);
	float c1z= 1.0f - 2.0f/3.0f*sintz*sintz;

	float k= sqrt(kx*kx + ky*ky + kz*kz);
	float shot_noise= c1x*c1y*c1z*nbar_inv; // C1 function in Jing 2005

	if(k <= knq) {

	// C2 function in Jing 2005
	  float c2gg= 0.0f, c2gu= 0.0f;
	for(int ax=-na; ax<na; ++ax) {
	  float kax= kx+ax*kp;
	  for(int ay=-na; ay<na; ++ay) {
	    float kay= ky+ay*kp;
	    for(int az=-na; az<na; ++az) {
	      float kaz= kz+az*kp;
	      float ka= sqrt(kax*kax + kay*kay + kaz*kaz);
	      float kk= ka/k;

	      float w1= w(sin_fac*(kx+ax*kp))*
		        w(sin_fac*(ky+ay*kp))*
                        w(sin_fac*(kz+az*kp));
	      float w2= w1*w1;
	      float w4= w2*w2;
	      float Pkfac= pow(ka/k, neff);;
	      c2gg += w4*Pkfac;
	      c2gu += w4*Pkfac/kk;
	      c2gu += w4*Pkfac/(kk*kk);
	    }
	  }
	}

	int index= ncz*(nc*ix + iy) + iz;
	float d_re= delta_k[2*index  ];
	float d_im= delta_k[2*index+1];

	float u_re= u_k[2*index  ];
	float u_im= u_k[2*index+1];

	Pgg.fill(k, (pk_fac*(d_re*d_re + d_im*d_im) - shot_noise)/c2gg);
	Pgu.fill(k,  pk_fac*(d_re*u_im - d_im*u_re)/c2gg);
	Puu.fill(k,  pk_fac*(u_re*u_re + u_im*u_im)/c2gg);

	Pgg_raw.fill(k, pk_fac*(d_re*d_re + d_im*d_im));
	Pgu_raw.fill(k, pk_fac*(d_re*u_im - d_im*u_re));
	Puu_raw.fill(k, pk_fac*(u_re*u_re + u_im*u_im));
	}
      }
    }
  }

  printf("# momentum power spectrum\n");
  printf("# Column 1: k\n");
  printf("# Column 2: nmode\n");
  printf("# Column 3: Pgg(k)\n");
  printf("# Column 4: Pgp(k)\n");
  printf("# Column 5: Ppp(k)\n");
  printf("# Column 6: Pgg(k) without shotnoise/alias correction\n");
  printf("# Column 7: Pgp(k) without alias correction\n");
  printf("# Column 8: Ppp(k) without alias correction\n")

  for(int j=0; j<Pgg.numbin(); ++j) {
    if(Pgg.n(j) > 0)
      printf("%e %d %e %e %e %e %e %e\n", 
	     Pgg.x(j), Pgg.n(j), Pgg.y(j), Pgu.y(j), Puu.y(j),
	     Pgg_raw.y(j), Pgu_raw.y(j), Puu_raw.y(j)
	     );
  }
}
