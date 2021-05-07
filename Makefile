# Makefile for pix_fft

# specify the location of the pure-data sources
PDDIR = /Users/fd/Documents/pure-data

lib.name = pix_fft

cflags = -std=c++11
# specify the location of FFTW header file
cflags += -Iinclude -I/usr/local/include
# specify the location of the GEM library
cflags += -I/Users/fd/Documents/Gem/src

# specify the location and name of the FFTW library
ldlibs = -L/usr/local/lib -lfftw3

class.sources = src/pix_fft.cpp

datafiles = help/pix_fft-help.pd README.txt LICENSE.txt

# define min version for macos to avoid multiple warnings
# - compiling FFTW on 10.14 generates these warnings 
define forDarwin
	cflags += -mmacosx-version-min=10.14
endef

# make a multi-object library
# make-lib-executable = yes

# provide path to pdlibbuilder dir
PDLIBBUILDERDIR=./pd-lib-builder
include $(PDLIBBUILDERDIR)/Makefile.pdlibbuilder