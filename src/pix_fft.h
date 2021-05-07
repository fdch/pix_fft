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

  virtual void    processGrayImage(imageStruct &image);
  void                    displayMess(int);
  void                    squelchMess(float);
  void                    normMess(float);

  //////////
  // The real and imaginary outlets
  t_outlet            *m_dataOut;
  
  t_atom               *m_buffer;

  unsigned char        *q1,*q2,*q3,*q4;
  
  double               *fftwIn,   *fftwOutR, *m_mag, *m_mag2;
  fftw_complex         *fftwOut,  *fftwInC;
  fftw_plan             fftwPlanF, fftwPlanB;

  int             m_size, m_xsize, m_ysize, m_insize;
  bool            m_enable, m_convolve;
  int             m_display;
  float           m_squelch, m_norm;
  
  void             reallocAll(int n, int m);
  void             copyRect(unsigned char*s,unsigned char*t, bool dir, bool Yoff, bool Xoff);
  void             shiftFFT(unsigned char*data);
  void             displayCallback(void *s, t_float val);
  //////////
  // Pass the data
  void            DATAMess(t_symbol*s, int argc, t_atom *argv);
  double          rsqrt(double x);

};
#endif
