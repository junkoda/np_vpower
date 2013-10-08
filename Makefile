# 1.  If you don't want/have the Boost command line options, set
#     USE_BOOST_PROGRAM_OPTIONS = 0.
# 2.  Review CXX, BOOST_DIR, FFTW3_DIR and modify them if necessary.
#     Single-precision FFTW3 library is required.
# 3a. make
# 3b. make momentum_power (optional; if you want matter momentum power spectra).

EXEC = np_vpower 

# Basic setups
CXX       ?= g++
BOOST_DIR ?= /opt/local
FFTW3_DIR ?= #/Users/junkoda/opt/gcc/fftw3/fftw-3.3.3  # for example

#
# Compile options
# USE_BOOST_PROGRAM_OPTIONS = 0
#   Program parameters are hard coded in np_vpower.cpp
# USE_BOOST_PROGRAM_OPTIONS = 1
#   Allows you to change parameter in command line rather than hard-coded values
#   Boost library is required
#
USE_BOOST_PROGRAM_OPTIONS = 1

DIR_PATH = $(BOOST_DIR) $(FFTW3_DIR)

ifeq ($(CXX), g++)
  # Optional warnings to detect potential bugs
  WOPT   ?= -Wall -W -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wmissing-prototypes -Wredundant-decls
endif

CXXFLAGS := -O2 $(WOPT)
CXXFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS     := -lm
LIBS     += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

ifeq ($(USE_BOOST_PROGRAM_OPTIONS),1)
  CXXFLAGS += -DUSE_BOOST_PROGRAM_OPTIONS
  LIBS     += -lboost_program_options
              # Could be -lboost_program_options-mt for MacPorts
endif


all: $(EXEC)

#
# Velocity Power Spectrum with the Nearest Particle Method
#
OBJS1 := np_vpower.o
OBJS1 += binary_tree.o nbr_finder.o fft_mesh.o
OBJS1 += nearest_nbr_mesh.o density_mesh.o power_spectra.o histogram.o
OBJS1 += halo_file.o transformation.o 

LIBS1  = $(LIBS) -lfftw3f

np_vpower: $(OBJS1)
	$(CXX) $(OBJS1) $(LIBS1) -o $@

# g++ -MM *.cpp
binary_tree.o: binary_tree.cpp binary_tree.h particle.h gadget_file2.h
density_mesh.o: density_mesh.cpp density_mesh.h particle.h gadget_file2.h
fft_mesh.o: fft_mesh.cpp fft_mesh.h
halo_file.o: halo_file.cpp halo_file.h particle.h gadget_file2.h
histogram.o: histogram.cpp histogram.h
nbr_finder.o: nbr_finder.cpp nbr_finder.h binary_tree.h particle.h \
  gadget_file2.h kth_value.h
nearest_nbr_mesh.o: nearest_nbr_mesh.cpp particle.h gadget_file2.h \
  binary_tree.h nbr_finder.h kth_value.h fft_mesh.h nearest_nbr_mesh.h
np_vpower.o: np_vpower.cpp particle.h gadget_file2.h halo_file.h \
  fft_mesh.h nearest_nbr_mesh.h power_spectra.h transformation.h \
  density_mesh.h Makefile
power_spectra.o: power_spectra.cpp histogram.h power_spectra.h
transformation.o: transformation.cpp transformation.h particle.h \
  gadget_file2.h


.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1) $(OBJS2)
