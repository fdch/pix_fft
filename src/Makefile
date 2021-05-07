# change to your local directories!
PD_DIR = /Users/fd/Documents/pure-data/src
GEM_DIR = /Users/fd/Documents/Gem/src

# build flags
INCLUDES = -I$(PD_DIR) -I$(GEM_DIR)
CPPFLAGS = -g -Wall -O2 -fPIC -Os -std=c++11 -pedantic
CPPFLAGS += -funroll-loops -ffast-math -mmmx -fpascal-strings -DDARWIN
LDFLAGS = -dynamiclib -undefined dynamic_lookup -L/usr/local/lib 
LIBS = -ldl -lz -lpthread -lm -lfftw3
EXTENSION = pd_darwin

SOURCES = pix_fft

all:
	g++ $(INCLUDES) $(CPPFLAGS) -o $(SOURCES).o -c $(SOURCES).cpp
	g++ -o $(SOURCES).$(EXTENSION) $(LDFLAGS) $(SOURCES).o $(LIBS)

clean:
	rm -f $(SOURCES)*.o
	rm -f $(SOURCES)*.$(EXTENSION)
