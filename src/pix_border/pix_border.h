#ifndef _INCLUDE__GEM_PIXES_pix_border_H_
#define _INCLUDE__GEM_PIXES_pix_border_H_

#include "Base/GemPixObj.h"


class GEM_EXTERN pix_border : public GemPixObj
{
  CPPEXTERN_HEADER(pix_border, GemPixObj)

  public:

  //////////
  // Constructor
  pix_border(t_floatarg n);

  protected:

  //////////
  // Destructor
  virtual ~pix_border(void);

  virtual void    processGrayImage(imageStruct &image);
  //////////
  // The real and imaginary outlets
  t_outlet            *m_dataOut;

  t_atom        *m_buffer;
  
  unsigned char *m_data;
  
  int   m_top, m_bottom, m_left, m_right, m_size;
  float m_thresh, m_bwidth;
  int m_contours;

  
  void  threshMess(float f);
  void  headMess(float f);
  void  leftMess(float f);
  void  bottomMess(float f);
  void  borderWidthMess(float f);
  void  rightMess(float f);
  void  contoursMess(int f);
};
#endif
