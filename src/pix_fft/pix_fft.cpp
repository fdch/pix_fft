////////////////////////////////////////////////////////
//
// pix_fft
//
// Effects using the Fourier Transform (FFTW)
// 
//
// fdch.github.io
//
// Fede Camara Halac 2021
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
#include <complex>
#define PLANFLAG FFTW_ESTIMATE
#define DEF 64
#define _USE_MATH_DEFINES
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
  m_magOut(0),m_angOut(0),
  m_enable(false), m_convolve(false),
  m_size(0), m_xsize(0), m_ysize(0), m_insize(0),
  m_filter(0), m_display(0), m_shift(1), m_output(0),
  m_filter_type(0), m_output_type(0), m_convolve_type(0),
  m_filter_size(4.0),
  m_magscale(1.0), m_magclip(0.0), m_pshift(1.0),
  m_squelch(10.0), m_norm(0.0), m_pow(0.5)
{
  n=n<=0?DEF:n;
  m_xsize=n;
  m_ysize=n;
  m_insize=n*n;
  m_size=n*(n/2+1);
  reallocAll(n,n);
  inlet_new(this->x_obj, 
            &this->x_obj->ob_pd, 
            gensym("list"),
            gensym("xaxis"));
  inlet_new(this->x_obj, 
            &this->x_obj->ob_pd, 
            gensym("list"),
            gensym("yaxis"));

  m_magOut = outlet_new(this->x_obj, &s_list);
  m_angOut = outlet_new(this->x_obj, &s_list);
}
/////////////////////////////////////////////////////////
// Destructor
//
/////////////////////////////////////////////////////////
pix_fft :: ~pix_fft()
{
  outlet_free(m_magOut);
  outlet_free(m_angOut);
  if(!m_size)return;
  else{
    delete [] fftwIn;
    delete [] fftwOutR;
    delete [] q1;
    delete [] q2;
    delete [] q3;
    delete [] q4;
    delete [] m_magBuf;
    delete [] m_angBuf;
    delete [] m_mag2;
    delete [] m_mag;
    delete [] m_angle2;
    delete [] m_angle;
    fftw_free(fftwOut);
    fftw_free(fftwInC);
    fftw_destroy_plan(fftwPlanF);
    fftw_destroy_plan(fftwPlanB);
  }
}
/////////////////////////////////////////////////////////
// (Re)allocate arrays
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
  m_angle  = new double [m_size];
  m_angle2 = new double [m_size];
  fftwIn   = new double [m_insize];
  fftwOutR = new double [m_insize];
  q1 = new unsigned char [m_insize/4];
  q2 = new unsigned char [m_insize/4];
  q3 = new unsigned char [m_insize/4];
  q4 = new unsigned char [m_insize/4];
  m_magBuf = new t_atom[m_size];
  m_angBuf = new t_atom[m_size];
  fftwOut = (fftw_complex *)fftw_alloc_complex(m_size);
  fftwInC = (fftw_complex *)fftw_alloc_complex(m_insize);
  fftwPlanF = fftw_plan_dft_r2c_2d(n, m, fftwIn,  fftwOut,  PLANFLAG);
  fftwPlanB = fftw_plan_dft_c2r_2d(n, m, fftwInC, fftwOutR, PLANFLAG);
  for(int j=0;j<m_size;j++) m_mag2[j] = 0.0, m_angle2[j] = 0.0;
  // Notify and enable computing
  post("m_insize=%d, m_size=%d", m_insize, m_size);
  m_enable=true;
}

/////////////////////////////////////////////////////////
// Helper routine to copy rectangles from one place to another
//
/////////////////////////////////////////////////////////
void pix_fft :: copyRect(unsigned char*s,unsigned char *t,bool dir,bool Yoff, bool Xoff)
{
// copyRect: 
// *s  *t  direction(1=s->t, 0=t->s) Y  X (offsets)
  unsigned char *src = &s[0];
  unsigned char *tar = &t[0];
  int xsize = m_xsize/2;//n
  int ysize = m_ysize/2;//m
  int Xoffset = Xoff?xsize:0;//add n/2 to start at mid column
  int Yoffset = Yoff?m_insize/2:0;//n*m/2 to start at mid row
  for(int i=0;i<ysize;i++)
    for(int j=0;j<xsize;j++) {
      int step=i*ysize*2+j+Xoffset+Yoffset;
      if (dir) tar[step] = *src++;
      else *src++ = tar[step];
    }
}

/////////////////////////////////////////////////////////
// Shift Zeroth frequency to center
//
/////////////////////////////////////////////////////////
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

void pix_fft :: copySymmetry(unsigned char *data)
{
  //copy the non-computed symmetry back to the data
  long step, i=0, j;
  for(;i<m_ysize;i++) {
    for(j=m_xsize/2+1;j<m_xsize;j++) {
      step = i * m_ysize + j;
      data[step] = data[m_insize-step];
    }
  }
}


/////////////////////////////////////////////////////////
// Fast Inverse square root (found somewhere online)
//
/////////////////////////////////////////////////////////
double pix_fft :: rsqrt(double x)
{
    double xhalf = 0.5 * x;
    long i = *(long*)&x;
    i = 0x5f3759df - (i >> 1);
    x = *(double*)&i;  x = x*(1.5f-(xhalf*x*x));
    return x;
}

bool pix_fft :: insideBox(long x,long y, float s) 
{
  int top    = std::floor (m_ysize / 2 - s);
  int bottom = std::ceil  (m_ysize / 2 + s);
  int left   = std::floor (m_xsize / 2 - s);
  int right  = std::ceil  (m_xsize / 2 + s);
  return ((x > left) && (x < right)) && ((y > top) && (y < bottom));
}

void pix_fft :: display(unsigned char *data, double *src)
{
  long i,j;
  for(i=0;i<m_ysize;i++) {
    for(j=0;j<m_xsize/2+1;j++) {
      data[i * m_ysize + j] = CLAMP(255 * std::log10(1.0 + *src++) );
    }
  }
  //copy the non-computed symmetry back to the data
  copySymmetry(data);
  // shift zero-th frequency to center
  if(m_shift) shiftFFT(data);
}

/////////////////////////////////////////////////////////
// Process image (grey space only)
//
/////////////////////////////////////////////////////////
void pix_fft :: processGrayImage(imageStruct &image)
{
  // Pointer to the pixels (unsigned char 0-255)
  unsigned char *pixels = image.data;
  int ysize = image.ysize;
  int xsize = image.xsize;
  if(!m_enable)return;
  // Check if sizes match and reallocate.
  if(m_insize!=ysize*xsize) {
    reallocAll(xsize, ysize);
  } else {
    //ROI
    int left=0;
    int top=0;
    int bottom=m_ysize;
    int right=m_xsize/2+1;
    //input to FFTW
    long i, j, k=0;
    long mx = m_xsize / 2, my = m_ysize / 2;
    double mag, angle, re, im, ctrl, rsize = 1.0 / m_insize;
    // double magmax=0.0, angmax=0.0, angmin=1e20, rmagmax;
    // double re, im, mag, rmag, ctrl;
    // fill the fft input array
    for (i=0; i<m_insize; i++) {
      fftwIn[i] = (double)pixels[image.csize * i]/255.;
    }
    // calculate the forward fft
    fftw_execute(fftwPlanF);

    // populate magnitude and angle arrays
    for(i=top;i<bottom;i++) {
      for(j=left;j<right;j++) {
        // get real and imaginary values from the fftw output
        re = fftwOut[k][0], im = fftwOut[k][1];
        // calculate magnitude and store it in atom buf
        mag = std::pow(re * re + im * im, m_pow);
        SETFLOAT(&m_magBuf[k], m_output_type?re:mag);
        // get magnitude max for normalization
        // magmax = mag>magmax?mag:magmax;
        // clip, scale, and/or power the magnitude
        m_mag[k] = std::fmax(m_magclip, mag) * m_magscale;
        // calculate the angle and store it in atom buf
        angle = std::atan2(im, re);
        SETFLOAT(&m_angBuf[k], m_output_type?im:angle);
        // get angle max and min for normalization
        // angmax = angle>angmax?angle:angmax;
        // angmin = angle<angmin?angle:angmin;
        // shift the angle by some value
        m_angle[k] = angle * m_pshift;
        k++;
      }
    }
    
    if(m_output>0) {
      outlet_list(m_magOut, gensym("list"), m_size, m_magBuf);
      outlet_list(m_angOut, gensym("list"), m_size, m_angBuf);
    } 
    // rmagmax = m_norm / magmax;
    // inv_pow = 1.0 / m_pow;

    // if (m_norm >= 1.0) {
    //   // post("rmax: %f\nmax:%f",rmagmax,magmax);
    //   for(i=0;i<m_size;i++) {
    //     m_mag[i] *= rmagmax;
    //   }
    // }
    switch(m_display) {
      case 0:
      default:
        //display original
      break;
      case 1:
        //display magnitude
        display(pixels, m_mag);
      break;
      case 2:
        //display phase
        display(pixels, m_angle);
      break;
      case 3:
        //reconstruction
        switch(m_filter) {
          case 0:
          default:
            //display backwards fft
            for(k=0,i=top;i<bottom;i++) {
              for(j=left;j<right;j++) {
                fftwInC[k][0] = std::cos(m_angle[k]) * m_mag[k];
                fftwInC[k][1] = std::sin(m_angle[k]) * m_mag[k];
                k++;
              }
            }
          break;
          case 1:
            //display reconstructed phase only
            for(k=0,i=top;i<bottom;i++) {
              for(j=left;j<right;j++) {
                fftwInC[k][0] = std::cos(m_angle[k]) * m_magscale;
                fftwInC[k][1] = std::sin(m_angle[k]) * m_magscale;
                k++;
              }
            }
          break;
          case 2:
            // lopass filter
            switch(m_filter_type) {
              case 0:
              default:
                // square filter
                // apply a square black filter at the center of fft
                for(k=0,i=top;i<bottom;i++) {
                  for(j=left;j<right;j++) {
                    if(insideBox(j,i,m_filter_size/2.0)){
                      // pixels[image.csize * i*m_ysize+j] = 0;
                      fftwInC[k][0] = fftwInC[k][1] = 0.0;
                    } else {
                      // pixels[image.csize * i*m_ysize+j] = 80;
                      fftwInC[k][0] = std::cos(m_angle[k]) * m_mag[k];
                      fftwInC[k][1] = std::sin(m_angle[k]) * m_mag[k];
                    }
                    k++;
                  }
                }
              break;
              case 1:
                //circle filter
                // apply a circle filter at the center of fft
                for(k=0,i=top;i<bottom;i++) {
                  for(j=left;j<right;j++) {
                    
                    re = std::cos(m_angle[k]) * m_mag[k];
                    im = std::sin(m_angle[k]) * m_mag[k];

                    if (insideBox(j,i,m_filter_size))
                      if (m_filter_size>std::sqrt((i-my)*(i-my)+(j-mx)*(j-mx))) 
                        re = 0.0, im = 0.0;
                    
                    fftwInC[k][0] = re;
                    fftwInC[k][1] = im;
                    k++;
                  }
                }
              break;  
            }
          break;
          case 3: //hipass filter
              switch(m_filter_type) {
              case 0:
              default:
                // square filter
                // apply a square black filter at the center of fft
                for(k=0,i=top;i<bottom;i++) {
                  for(j=left;j<right;j++) {
                    if(!insideBox(j,i,m_filter_size/2.0)){
                      // pixels[image.csize * i*m_ysize+j] = 0;
                      fftwInC[k][0] = fftwInC[k][1] = 0.0;
                    } else {
                      // pixels[image.csize * i*m_ysize+j] = 80;
                      fftwInC[k][0] = std::cos(m_angle[k]) * m_mag[k];
                      fftwInC[k][1] = std::sin(m_angle[k]) * m_mag[k];
                    }
                    k++;
                  }
                }
              break;
              case 1:
                //circle filter
                // apply a circle filter at the center of fft
                for(k=0,i=top;i<bottom;i++) {
                  for(j=left;j<right;j++) {
                    
                    re = 0.0;
                    im = 0.0;

                    if (insideBox(j,i,m_filter_size)) {
                      if (m_filter_size>std::sqrt((i-my)*(i-my)+(j-mx)*(j-mx))) {
                        re = std::cos(m_angle[k]) * m_mag[k];
                        im = std::sin(m_angle[k]) * m_mag[k];
                      }
                    }
                    fftwInC[k][0] = re;
                    fftwInC[k][1] = im;
                    k++;
                  }
                }
              break;  
            }
          break; // end hipass case 3
          case 4:
            // mask (noise removal in audio world)
            for(k=0,i=top;i<bottom;i++) {
              for(j=left;j<right;j++) {
                mag = m_mag[k] * m_mag2[k];
                // ctrl = mag - m_mag2[k];
                // ctrl = m_mag[k]>0.0?0.5:0.0;
                // ctrl = std::sqrt(ctrl / (mag + 1e20));
                fftwInC[k][0] = std::cos(m_angle[k]) * mag;
                fftwInC[k][1] = std::sin(m_angle[k]) * mag;
                k++;
              }
            }
          break;
          // case 4: //bandpass filter
          // break;
          // case 5: //bandreject filter
          // break;
        }
        // copySymmetry(pixels);
        fftw_execute(fftwPlanB);
        for(i=top;i<m_insize;i++) pixels[i] = CLAMP(255 * fftwOutR[i] * rsize);

      break; //end m_display case 3
      case 4:
        // reconstruction with inputs
        if(!m_convolve) break;

        switch(m_convolve_type) {
          case 0:
          default:
            //reconstruct only the inlets (mag phase, aka output_type 0)
            for(k=0,i=top;i<bottom;i++) {
              for(j=left;j<right;j++) {
                angle = m_angle2[k] * m_pshift;
                mag = std::fmax(m_magclip, m_mag2[k]) * m_magscale;
                fftwInC[k][0] = std::cos(angle) * mag;
                fftwInC[k][1] = std::sin(angle) * mag;
                k++;
              }
            }
          break;
          case 1:
            //reconstruct only the inlets (real imag, aka output_type 1)
            for(k=0,i=top;i<bottom;i++) {
              for(j=left;j<right;j++) {
                fftwInC[k][0] = m_mag2[k];
                fftwInC[k][1] = m_angle2[k];
                k++;
              }
            }
          break;
          case 2:
            // timbre-stamp in audio world
            for(k=0,i=top;i<bottom;i++) {
              for(j=left;j<right;j++) {
                // reciprocal modulus of pixels
                ctrl = rsqrt(m_mag[k]+1e20);
                // clip it from 0 to the squelch parameter
                ctrl = ctrl<0.0?0.0:ctrl>m_squelch?m_squelch:ctrl;
                // and multiply by modulus of control amplitude
                ctrl *= m_mag2[k];
                fftwInC[k][0] = std::cos(m_angle[k]) * ctrl;
                fftwInC[k][1] = std::sin(m_angle[k]) * ctrl;
                k++;
              }
            }
          break;
        } // end m_convolve_type switch
      
        // computer backwards fft

        fftw_execute(fftwPlanB);
        for(i=top;i<m_insize;i++) pixels[i] = CLAMP(255 * fftwOutR[i] * rsize);
      
      break;//end m_display case 4
    }// end m_display switch
  }// end else statement
}
/////////////////////////////////////////////////////////
// xaxisMess
//
/////////////////////////////////////////////////////////
void pix_fft :: xaxisMess(t_symbol *s, int argc, t_atom *argv)
{
  if(argc!=m_size) {
    post("%s function says: FFT data mismatch.", s->s_name);
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
/////////////////////////////////////////////////////////
// yaxisMess
//
/////////////////////////////////////////////////////////
void pix_fft :: yaxisMess(t_symbol *s, int argc, t_atom *argv)
{
  if(argc!=m_size) {
    post("%s function says: FFT data mismatch.", s->s_name);
    post("%d != %d", argc, m_size);
    m_convolve=false;
    return;
  } else {
    for(int i=0;i<argc;i++) {
      m_angle2[i] = (double)(atom_getfloat(&argv[i]));
    }
    m_convolve=true;
  }
}
void pix_fft :: displayMess(int f)
{
  m_display=f;
}
void pix_fft :: filterMess(int f)
{
  m_filter=f;
}
void pix_fft :: shiftMess(int f)
{
  m_shift=f;
}
void pix_fft :: outputMess(int f)
{
  m_output=f;
}
void pix_fft :: filterTypeMess(int f)
{
  m_filter_type=f;
}
void pix_fft :: filterSizeMess(float f)
{
  m_filter_size = f<0.0?0.0:f>m_xsize?m_xsize:f>m_ysize?m_ysize:f;
  // post("%.2f",m_filter_size);
}
void pix_fft :: squelchMess(float f)
{
  m_squelch=0.01*f*f;
}
void pix_fft :: magScaleMess(float f)
{
  m_magscale=f;
}
void pix_fft :: magClipMess(float f)
{
  m_magclip=f;
}
void pix_fft :: phaseMess(float f)
{
  m_pshift=f;
}
void pix_fft :: powerMess(float f)
{
  m_pow=f;
}
void pix_fft :: normMess(float f)
{
  m_norm=f;
}
void pix_fft :: outputTypeMess(int f)
{
  m_output_type=f;
}
void pix_fft :: convolveTypeMess(int f)
{
  m_convolve_type=f;
}
void pix_fft :: bangMess(void)
{
  if(!m_enable) return;
  if(!m_output) {
    //only output oneshot if frame output is disabled
    outlet_list(m_magOut, gensym("list"), m_size, m_magBuf);
    outlet_list(m_angOut, gensym("list"), m_size, m_angBuf);
  }
}
void pix_fft :: maskMess(t_symbol *s, int argc, t_atom *argv)
{
  if(!m_enable) return;
  m_convolve=false;
  float f=0.0, g=100000.0;
  if(argc==2){
    f = atom_getfloat(&argv[0]);
    g = atom_getfloat(&argv[1]);
  } else if (argc==1) {
    f = atom_getfloat(&argv[0]);
  } else if(argc>2)  {
    post("%s: ignoring extra args.", s->s_name); 
  } else {
    post("using default min and max values:%.2f,%.2f",f,g);
  }
  for(int i=0;i<m_size;i++) {
      m_mag2[i]   = m_mag[i]<f?0.0:m_mag[i]>g?g:m_mag[i];
      m_angle2[i] = m_angle[i];
  }
  m_convolve=true;
}
/////////////////////////////////////////////////////////
// static member functions
//
/////////////////////////////////////////////////////////
void pix_fft :: obj_setupCallback(t_class *classPtr) 
{
  CPPEXTERN_MSG0(classPtr, "bang", bangMess);
  CPPEXTERN_MSG1(classPtr, "display", displayMess, int);
  CPPEXTERN_MSG1(classPtr, "shift", shiftMess, int);
  CPPEXTERN_MSG1(classPtr, "filter", filterMess, int);
  CPPEXTERN_MSG1(classPtr, "output", outputMess, int);
  CPPEXTERN_MSG1(classPtr, "filter_size", filterSizeMess, float);
  CPPEXTERN_MSG1(classPtr, "convolve_type", convolveTypeMess, int);
  CPPEXTERN_MSG1(classPtr, "filter_type", filterTypeMess, int);
  CPPEXTERN_MSG1(classPtr, "output_type", outputTypeMess, int);
  CPPEXTERN_MSG1(classPtr, "squelch", squelchMess, float);
  CPPEXTERN_MSG1(classPtr, "magscale", magScaleMess, float);
  CPPEXTERN_MSG1(classPtr, "magclip", magClipMess, float);
  CPPEXTERN_MSG1(classPtr, "phase", phaseMess, float);
  CPPEXTERN_MSG1(classPtr, "power", powerMess, float);
  CPPEXTERN_MSG1(classPtr, "norm", normMess, float);
  CPPEXTERN_MSG (classPtr, "mask", maskMess);
  CPPEXTERN_MSG (classPtr, "xaxis", xaxisMess);
  CPPEXTERN_MSG (classPtr, "yaxis", yaxisMess);
}