#ifndef _INCLUDE__GEM_PIXES_pix_fft_H_
#define _INCLUDE__GEM_PIXES_pix_fft_H_

#include "Base/GemPixObj.h"
#include "fftw3.h"


class GEM_EXTERN pix_fft : public GemPixObj
{
  CPPEXTERN_HEADER(pix_fft, GemPixObj)

  public:

  //////////
  // Constructor
  pix_fft(t_floatarg n);

  protected:

  //////////
  // Destructor
  virtual ~pix_fft(void);

  //////////
  // process the image
  virtual void  processGrayImage(imageStruct &image);

  //////////
  // the outlet
  t_outlet  *m_magOut, *m_angOut; 
  t_atom  *m_magBuf, *m_angBuf;
  
  //////////
  // the fftw arrays
  fftw_complex  *fftwOut, *fftwInC;
  fftw_plan  fftwPlanF, fftwPlanB;

  //////////
  // the quadrants
  unsigned char  *q1, *q2, *q3, *q4;

  //////////
  // the data arrays
  double *fftwIn, *fftwOutR;
  double *m_mag, *m_angle;
  double *m_mag2, *m_angle2;
  
  //////////
  // variables
  bool  m_enable, m_convolve;
  int   m_size, m_xsize, m_ysize, m_insize;
  int   m_filter, m_display, m_shift, m_output;
  int   m_filter_type, m_output_type, m_convolve_type;
  float m_filter_size;
  float m_magscale, m_magclip, m_pshift;
  float m_squelch, m_norm, m_pow;
  
  //////////
  // internal fuunctions 
  double  rsqrt(double x);
  void  reallocAll(int n, int m);
  bool  insideBox(long x,long y,float s);
  void  copyRect(unsigned char*s,unsigned char*t,bool dir,bool Yoff,bool Xoff);
  void  shiftFFT(unsigned char*data);
  void  display(unsigned char*data, double *src);
  void  copySymmetry(unsigned char*data);
  //////////
  // callbacks
  void  bangMess(void);
  void  xaxisMess(t_symbol*s, int argc, t_atom *argv);
  void  yaxisMess(t_symbol*s, int argc, t_atom *argv);
  void  displayMess(int);
  void  shiftMess(int);
  void  filterMess(int);
  void  filterTypeMess(int);
  void  convolveTypeMess(int);
  void  outputTypeMess(int);
  void  outputMess(int);
  void  filterSizeMess(float);
  void  squelchMess(float);
  void  magScaleMess(float);
  void  magClipMess(float);
  void  phaseMess(float);
  void  powerMess(float);
  void  normMess(float);
  void  maskMess(t_symbol*s, int argc, t_atom *argv);
};
#endif
