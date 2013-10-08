#ifndef POWER_SPECTRA_H
#define POWER_SPECTRA_H 1

void calc_power_spectrum_sa(const int nc, const float boxsize, const float delta_k[], const float u_k[], const float nbar, const float neff, const float dk, float kmax);

void calc_momentum_power_sa(const int nc, const float boxsize, const float delta_k[], const float u_k[], const float nbar, const float neff, const float dk, float kmax);

#endif
