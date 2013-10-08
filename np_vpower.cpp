//
// Main code for the np_vpower
//
// This code calculates matter (or halo) and line-of-sight 
// velocity auto- and cross-power spectrum 
// in real or redshift space (--redshift-space)
//
// Also calculates matter (or halo) and momentum power spectra (--momentum)
//
// $ np_vpower [options] <particle/halo filename>
//
#include <iostream>
#include <vector>

#include "particle.h"
#include "halo_file.h"
#include "fft_mesh.h"
#include "nearest_nbr_mesh.h"
#include "gadget_file2.h"
#include "power_spectra.h"
#include "transformation.h"
#include "density_mesh.h"

#ifdef USE_BOOST_PROGRAM_OPTIONS
  #include <boost/program_options.hpp>
  using namespace boost::program_options;
#endif

using namespace std;

enum InputFileType {
  gadget_binary, fof_text
};

int main(int argc, char* argv[])
{
#ifdef USE_BOOST_PROGRAM_OPTIONS
  //
  // Command line options and default values
  //
  options_description opt("vpower_spectrum5 [options] <filename>");
  opt.add_options()
    ("help,h", "display this help")
    ("fof-text", "Input is a FoF text file")
    ("gadget-binary", "Input is a Gadget binary file")
    ("nc", value<int>()->default_value(128), "number of density mesh per dim")
    ("boxsize", value<float>()->default_value(1000.0f), 
                                             "boxsize (for subfind case only)")
    ("redshift-space,r", "redshift-space distortion in z direction")
    ("z", value<float>()->default_value(0.0f), "redshift for --redshift-space")
    ("omegam", value<float>()->default_value(0.273, "0.273"), "omega_m necessary for --redshift-space")
    ("logMmin", value<double>()->default_value(1,"1"), 
                 "log10 Minimum halo mass")
    ("logMmax", value<double>()->default_value(20,"20"), 
                 "log10 Maximum halo mass")
    ("m", value<float>()->default_value(0.0f), "particle mass for FoF file")
    ("dk", value<float>()->default_value(0.01f, "0.01"), "output k bin widtth")
    ("kmax", value<float>()->default_value(0.0f, "0"), "output kmax (default kNq)")
    ("neff", value<float>()->default_value(-1.6f, "-1.6"), "Power spectrum beyond kNq")
    ("momentum-power", "Momemtum power spectra, instead of velocity")
    ;
  
  positional_options_description p;
  p.add("filename", -1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);
  
  if(vm.count("help") || ! vm.count("filename")) {
    cout << opt << "\n"; 
    return 0;
  }

  const bool redshift_space_distortions= vm.count("redshift-space");
  const bool momentum_power_spectra= vm.count("momentum-power");
  const int nc= vm["nc"].as<int>(); assert(nc > 0);
  float boxsize = vm["boxsize"].as<float>(); 
  float z= vm["z"].as<float>();
  float omega_m= vm["omegam"].as<float>();

  string filename= vm["fof-text"].as<string>();

  const float m= vm["m"].as<float>(); assert(m > 0.0f);
  const float logMmin= vm["logMmin"].as<float>();
  const float logMmax= vm["logMmax"].as<float>();

  const float neff= vm["neff"].as<float>();
  const float dk= vm["dk"].as<float>();
  const float kmax= vm["kmax"].as<float>();

  InputFileType input_file_type;
  if(vm.count("gadget-binary"))
    input_file_type= gadget_binary;
  else if(vm.count("fof-text"))
    input_file_type= fof_text;
  else {
    cerr <<"Error: input file type not specified --gadget-binary --fof-text?\n";
    return 1;
  }
#else
  //
  // Hard coded parameters without USE_BOOST_PROGRAM_OPTIONS
  //
  InputFileType input_file_type= gadget_binary;
  const int nc= 128;                 // number of grid cells per dimension
  const bool redshift_space_distortions= false;
  const bool momentum_power_spectra= false;

  // Parameters for FoF 
  //   (these parameters are automatically set for gadget-binary)
  float boxsize = 1000.0f;
  float m= 0;
  const float logMmin= 1.0;
  const float logMmax= 20.0;
  float z= 0.0f;                     // for FoF in redshift space only
  float omega_m= 0.273f;             // for FoF in redshift space only

  // File
  string filename= string(argv[argc-1]);

  // Details
  const float neff= -1.6f;
  const float dk= 0.01f;
  const float kmax= 0.0f;            // kmax=kNq if it is 0.
#endif

  //
  // Read particles
  //
  vector<ParticleData> v;
  float nbar= 0.0f;

  switch(input_file_type) {
   case gadget_binary: {
    cout << "# gadget-binary: " << filename << endl;

    gadget_file<particle_data_sph_all, ParticleData> gf;
    gf.set_velocity_conversion(true);

    gf.set_cdm(&v, 1);
    if(!gf.read(filename.c_str())) {
      cerr << "Unable to read gadget file: " << filename << endl;
      throw filename;
    }

    gadget_header1* header= gf.get_header();
    omega_m= header->omega0;
    z= header->redshift;
    boxsize= header->boxsize;
    // nbar= 0.0f, which means shotnoise is not subtracted for 
    // dark matter particles
    }
    break;

   case fof_text:
    cout << "# fof-text: " << filename << endl;
    if(m > 0.0f)
      printf("# Halo mass range log10(M): %4.2e - %4.2e\n", logMmin, logMmax); 
    else
      printf("# Reading all haloes (because m=0)\n");

    read_fof_text(filename.c_str(), v, m, logMmin, logMmax);

    nbar= v.size() / (boxsize*boxsize*boxsize);
    break;
  }

  if(v.empty()) {
    cerr << "Error: Zero particles\n";
    return 1;
  }

  // Redshift-space distortions (line of sight is along the z axis)
  if(redshift_space_distortions) {
    redshift_space_distortion(v, z, omega_m);
    printf("# Redshift space: z=%.4f omega_m=%.4f\n", z, omega_m);
  }
  //
  // Mesh
  //
  assert(boxsize > 0.0f);

  FFTmesh dmesh(nc), vmesh(nc);

  calculate_cic_density_mesh(v, boxsize, nc, dmesh.data());

  if(! momentum_power_spectra) {
    // Calculate velocity field
    calculate_nearest_particle_u(v, boxsize, nc, vmesh.data());

    // FFT
    cerr << "FFT...\n";
    dmesh.fft(); vmesh.fft();
    
    // Power spectrum calculation
    calc_power_spectrum_sa(nc, boxsize, dmesh.data(), vmesh.data(),
			   nbar, neff, dk, kmax);
  }
  else {
    // Calculate momentu field
    calculate_cic_momentum_mesh(v, boxsize, nc, vmesh.data());

    // FFT
    cerr << "FFT...\n";
    dmesh.fft(); vmesh.fft();

    // Power spectrum calculation
    const float neff= -1.6f;
    const float dk= vm["dk"].as<float>();
    const float kmax= vm["kmax"].as<float>();

    calc_momentum_power_sa(nc, boxsize, dmesh.data(), vmesh.data(),
			   nbar, neff, dk, kmax);
  }

  return 0;
}

