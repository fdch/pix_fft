////////////////////////////////////////////////////////
//
// pix_fft
//
// Calculates the Forward Fourier Transform using FFTW
// 
//
// fdch.github.io/tv
// camarafede@gmail.com
// Fede Camara Halac 2017
//
//
//
// GEM - Graphics Environment for Multimedia
//
// zmoelnig@iem.kug.ac.at
//
// Implementation file
//
//    Copyright (c) 1997-1998 Mark Danks.
//    Copyright (c) Günther Geiger.
//    Copyright (c) 2001-2011 IOhannes m zmölnig. forum::für::umläute. IEM. zmoelnig@iem.at
//    For information on usage and redistribution, and for a DISCLAIMER OF ALL
//    WARRANTIES, see the file, "GEM.LICENSE.TERMS" in this distribution.
//
/////////////////////////////////////////////////////////
#include "pix_fft.h"
#include "Utils/Functions.h"//for CLAMP
#include <cmath>
#define PLANFLAG FFTW_ESTIMATE
#define DEF 64

CPPEXTERN_NEW_WITH_ONE_ARG(pix_fft, t_floatarg, A_DEFFLOAT);
/////////////////////////////////////////////////////////
//
// pix_fft
//
/////////////////////////////////////////////////////////
// Constructor
//
/////////////////////////////////////////////////////////
pix_fft :: pix_fft(t_floatarg n):
  m_dataOut(0),
  m_size(0), 
  m_xsize(0), 
  m_ysize(0), 
  m_insize(0), 
  m_enable(false),
  m_convolve(false),
  m_display(1),
  m_squelch(10),
  m_norm(0.0001)
{
  n=n<=0?DEF:n;
  m_xsize=n;
  m_ysize=n;
  m_insize=n*n;
  m_size=n*(n/2+1);
  reallocAll(n,n);
  inlet_new(this->x_obj, &this->x_obj->ob_pd, gensym("list"),
            gensym("data"));

  m_dataOut = outlet_new(this->x_obj, &s_list);
}
/////////////////////////////////////////////////////////
// Destructor
//
/////////////////////////////////////////////////////////
pix_fft :: ~pix_fft()
{
  outlet_free(m_dataOut);
  if(!m_size)return;
  else{
    delete [] fftwIn;
    delete [] fftwOutR;
    delete [] q1;
    delete [] q2;
    delete [] q3;
    delete [] q4;
    delete [] m_buffer;
    delete [] m_mag2;
    delete [] m_mag;
    fftw_free(fftwOut);
    fftw_free(fftwInC);
    fftw_destroy_plan(fftwPlanF);
    fftw_destroy_plan(fftwPlanB);
  }
}
/////////////////////////////////////////////////////////
// Utility functions
//
/////////////////////////////////////////////////////////
void pix_fft :: reallocAll(int n, int m)
{
  m_enable=false;
  // Get new sizes
  m_xsize = n;
  m_ysize = m;
  m_insize = n*m; //actual size of image
  m_size = n*(m/2+1); //FFTW output size
  // Allocate arrays
  m_mag    = new double [m_size];
  m_mag2   = new double [m_size];
  fftwIn   = new double [m_insize];
  fftwOutR = new double [m_insize];
  q1 = new unsigned char [m_insize/4];
  q2 = new unsigned char [m_insize/4];
  q3 = new unsigned char [m_insize/4];
  q4 = new unsigned char [m_insize/4];
  m_buffer = new t_atom[m_size];
  fftwOut = (fftw_complex *)fftw_alloc_complex(m_size);
  fftwInC = (fftw_complex *)fftw_alloc_complex(m_insize);
  fftwPlanF = fftw_plan_dft_r2c_2d(n, m, fftwIn,  fftwOut,  PLANFLAG);
  fftwPlanB = fftw_plan_dft_c2r_2d(n, m, fftwInC, fftwOutR, PLANFLAG);
  // Notify and enable computing
  post("m_insize=%d, m_size=%d", m_insize, m_size);
  m_enable=true;
}

void pix_fft :: copyRect(unsigned char*s,unsigned char *t,bool dir,bool Yoff, bool Xoff)
{
// copyRect: 
// *s  *t  direction(1=s->t, 0=t->s) Y  X (offsets)
  unsigned char *src = s;//safe local copies
  unsigned char *tar = t;
  int cols = m_xsize/2;//n
  int rows = m_ysize/2;//m
  int Xoffset = Xoff?cols:0;//add n/2 to start at mid column
  int Yoffset = Yoff?m_insize/2:0;//n*m/2 to start at mid row
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++) {
      int step=i*rows*2+j+Xoffset+Yoffset;
      if (dir) tar[step] = *src++;
      else *src++ = tar[step];
    }
}
void pix_fft :: shiftFFT(unsigned char*data)
{
  //original
  copyRect(q1,data, 0, 0, 0);
  copyRect(q2,data, 0, 0, 1);
  copyRect(q3,data, 0, 1, 1);
  copyRect(q4,data, 0, 1, 0);
  //swapped
  copyRect(q1,data, 1, 1, 1);
  copyRect(q2,data, 1, 1, 0);
  copyRect(q3,data, 1, 0, 0);
  copyRect(q4,data, 1, 0, 1);
}

double pix_fft :: rsqrt(double x)
{
    double xhalf = 0.5 * x;
    long i = *(long*)&x;
    i = 0x5f3759df - (i >> 1);
    x = *(double*)&i;  x = x*(1.5f-(xhalf*x*x));
    return x;
}

/////////////////////////////////////////////////////////
// Process image (grey space only)
//
/////////////////////////////////////////////////////////
void pix_fft :: processGrayImage(imageStruct &image)
{
  // Pointer to the pixels (unsigned char 0-255)
  unsigned char *pixels = image.data;
  int rows = image.ysize;
  int cols = image.xsize;
  if(!m_enable)return;
  // Check if sizes match and reallocate.
  if(m_insize!=rows*cols) {
    reallocAll(cols, rows);
  } else {
    //input to FFTW
    long i, j, k=0, step;
    double re, im, mag, rmag, ctrl;
    for (i=0; i<m_insize; i++) fftwIn[i] = (double)pixels[i]/255.;

    fftw_execute(fftwPlanF);

    switch(m_display) {
      case 0:
        // post("List output only");
        for(i=0;i<m_ysize;i++)
          for(j=0;j<m_xsize/2+1;j++) {
            re = fftwOut[k][0];
            im = fftwOut[k][1];
            mag = sqrt(re*re + im*im);
            SETFLOAT(&m_buffer[i], mag);
          }
        outlet_list(m_dataOut, gensym("list"), m_size, m_buffer);
        break;
      case 1:
      default:
        // post("Magnitude display to screen");
        // calculate magnitude
        for(i=0;i<m_ysize;i++)
          for(j=0;j<m_xsize/2+1;j++) {
            re = fftwOut[k][0];
            im = fftwOut[k][1];
            mag = sqrt(re*re + im*im);
            pixels[i*m_ysize+j] = CLAMP(255.0 * log10(1. + mag) * m_norm);
            k++;
        }
        //copy the non-computed symmetry back to the data
        for(i=0;i<m_ysize;i++)
          for(j=m_xsize/2+1;j<m_xsize;j++) {
            step=i*m_ysize+j;
            pixels[step] = pixels[m_insize-step];
          }
        // shift zero-th frequency to center
        shiftFFT(pixels);
      break;
      case 2:
        // timbre-stamp 
        if(m_convolve) {
          // post("Convolving");
          for(i=0;i<m_ysize;i++) {
            for(j=0;j<m_xsize/2+1;j++) {
              re = fftwOut[k][0];
              im = fftwOut[k][1];
              mag = rsqrt(re*re + im*im + 1e-20);
              rmag = mag<0.0?0.0:0.01*m_squelch*m_squelch;
              ctrl = rmag*m_mag2[k]*m_norm;
              fftwInC[k][0] = re * ctrl;
              fftwInC[k][1] = im * ctrl;
              k++;
            }
          }
          // post("executing IFFT");
          fftw_execute(fftwPlanB);
          // post("Filling pixels with ifft");
          for(i=0;i<m_insize;i++) {
              // pixels[i] = CLAMP(255*fftwOutR[i]/m_insize);
              pixels[i] = CLAMP(255*log10(1e-20+fftwOutR[i]/m_size));
          }

        }
      break;
      case 3:
          double mask;
          // post("Convolving");
          for(i=0;i<m_ysize;i++) {
            for(j=0;j<m_xsize/2+1;j++) {
              re = fftwOut[k][0];
              im = fftwOut[k][1];
              mag = re*re + im*im;
              mask = mag - m_mag2[k] * m_squelch;
              mask = mask<0.0?0.0:mask;
              ctrl = sqrt(mask / (mag + 1e-20)) * m_norm;
              fftwInC[k][0] = re * ctrl;
              fftwInC[k][1] = im * ctrl;
              k++;
            }
          }
          // post("executing IFFT");
          fftw_execute(fftwPlanB);
          // post("Filling pixels with ifft");
          for(i=0;i<m_insize;i++) {
              // pixels[i] = CLAMP(255*fftwOutR[i]/m_insize);
              pixels[i] = CLAMP(255*log10(1e-20+fftwOutR[i]/m_size));
          }
      break;
    }
  }
}
/////////////////////////////////////////////////////////
// DATAMess
//
/////////////////////////////////////////////////////////
void pix_fft :: DATAMess(t_symbol *s, int argc, t_atom *argv)
{
  if(argc!=m_size) {
    post("FFT data mismatch.");
    post("%d != %d", argc, m_size);
    m_convolve=false;
    return;
  } else {
    for(int i=0;i<argc;i++) {
      m_mag2[i] = (double)(atom_getfloat(&argv[i]));
    }
    m_convolve=true;
  }


}
void pix_fft :: displayMess(int f)
{
  m_display=f;
}

void pix_fft :: squelchMess(float f)
{
  m_squelch=f;
}

void pix_fft :: normMess(float f)
{
  m_norm=f;
}
/////////////////////////////////////////////////////////
// static member function
//
/////////////////////////////////////////////////////////
void pix_fft :: obj_setupCallback(t_class *classPtr) 
{
  // CPPEXTERN_MSG0(classPtr, "bang", trigger);
  CPPEXTERN_MSG1(classPtr, "display", displayMess, int);
  CPPEXTERN_MSG1(classPtr, "squelch", squelchMess, float);
  CPPEXTERN_MSG1(classPtr, "norm", normMess, float);
  CPPEXTERN_MSG (classPtr, "data", DATAMess);
}