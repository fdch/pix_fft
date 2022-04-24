////////////////////////////////////////////////////////
//
// pix_border
//
// Draws contours on thresholded pixes
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
#include "pix_border.h"
// #include "Utils/Functions.h"//for CLAMP
// #include <cmath>
// #include <complex>
// #define PLANFLAG FFTW_ESTIMATE
// #define DEF 64

CPPEXTERN_NEW_WITH_ONE_ARG(pix_border, t_floatarg, A_DEFFLOAT);
/////////////////////////////////////////////////////////
//
// pix_border
//
/////////////////////////////////////////////////////////
// Constructor
//
/////////////////////////////////////////////////////////
pix_border :: pix_border(t_floatarg n):
  m_top(0), 
  m_bottom(240), 
  m_left(0), 
  m_right(240),
  m_thresh(0.f), m_bwidth(3.0),
  m_contours(0)
{
  if (n) m_thresh = n;
  m_dataOut = outlet_new(this->x_obj, &s_list);
  m_size = std::abs(m_bottom-m_top)*std::abs(m_left-m_right);
  m_data = new unsigned char [m_size];
}
/////////////////////////////////////////////////////////
// Destructor
//
/////////////////////////////////////////////////////////
pix_border :: ~pix_border()
{
  outlet_free(m_dataOut);
  delete [] m_data;
}
/////////////////////////////////////////////////////////
// Process image (grey space only)
//
/////////////////////////////////////////////////////////
void pix_border :: processGrayImage(imageStruct &image)
{
  // shamelessly taken from pix_mano
  int xsize = image.xsize;
  int ysize = image.ysize;
  unsigned char *base = image.data;


  int xcount, ycount, xcoord, ycoord, ycoordm1, ycoordm2;
  int pix12, pix21, pix22, pix23, pix32, edge;

  // Limits for search
  m_right=m_right>=xsize?xsize:m_right;  
  m_bottom=m_bottom>=ysize?ysize:m_bottom;  

  int area = std::abs(m_bottom-m_top)*std::abs(m_left-m_right);

  if (m_size != area) {
    delete [] m_data;
    m_data = new unsigned char [area];
    m_size = area;
    post("Resized to current area of size %d pixels.", area);
  }

  // IMAGE TREATMENT, border, thresh and edge detection

  // ****************** make BLACK BORDER....
  for (ycount = m_top; ycount <= m_bottom; ycount++) {
    ycoord = image.csize * xsize * ycount;
    if ( (ycount <= (m_top + m_bwidth)) || 
         (ycount >= (m_bottom - m_bwidth)) ) {
      for (xcount = m_left; xcount <= m_right; xcount++) {
        xcoord = image.csize * xcount + ycoord;
        base[xcoord] =    0;
      }
    }
    else {
      for (xcount = m_left; xcount <= m_left + m_bwidth; xcount++) {
        xcoord = image.csize * xcount + ycoord;
        base[xcoord] =    0;
      }
      for (xcount = m_right - m_bwidth; xcount <= m_right; xcount++) {
        xcoord = image.csize * xcount + ycoord;
        base[xcoord] =    0;
      }
    }
  }
  if(m_contours) {
    // ****************** RASTER
    for (ycount = m_top; ycount <= m_bottom; ycount++) {
      // only specifying coords for further calculations
      ycoord    = image.csize * xsize * ycount;
      ycoordm1  = image.csize * xsize * (ycount - 1);
      ycoordm2  = image.csize * xsize * (ycount - 2);  
      for (xcount = m_left; xcount <= m_right; xcount++) {
        xcoord    = image.csize * xcount + ycoord;
        // THRESH the IMAGE
        base[xcoord] = base[xcoord] < m_thresh ? 0 : 80;
        // prevent calculations outside area
        if (  (xcount > (m_left    + 2)) &&
              (xcount < (m_right  - 1))  &&
              (ycount > (m_top    + 2))  &&
              (ycount < (m_bottom - 2))    )  {
          // matrix, EDGE.
          pix12  =  base[image.csize * (xcount - 2) + ycoordm2];
          pix21  =  base[image.csize * (xcount - 3) + ycoordm1];
          pix22  =  base[image.csize * (xcount - 2) + ycoordm1];
          pix23  =  base[image.csize * (xcount - 1) + ycoordm1];
          pix32  =  base[image.csize * (xcount - 2) + ycoord  ];

          edge = ((pix12 * -1) + (pix21 * -1) + (pix22 * 4) + (pix23 * -1) + (pix32 * -1));

          edge = edge >= 80 ? 255: 0;

          base[image.csize * (xcount - 2) + ycoordm2] =    edge;
        }
      }
    }
  }
}

void pix_border :: contoursMess(int f)
{
  m_contours=f;
}
void pix_border :: threshMess(float f)
{
  m_thresh=f;
}
void pix_border :: headMess(float f)
{
  m_top=f<0?0:f;
}
void pix_border :: leftMess(float f)
{
  m_left=f<0?0:f;
}
void pix_border :: bottomMess(float f)
{
  m_bottom=f<0?0:f;
}
void pix_border :: rightMess(float f)
{
  m_right=f<0?0:f;
}
void pix_border :: borderWidthMess(float f)
{
  m_bwidth=f<0?0:f;
}
/////////////////////////////////////////////////////////
// static member function
//
/////////////////////////////////////////////////////////
void pix_border :: obj_setupCallback(t_class *classPtr) 
{
  CPPEXTERN_MSG1(classPtr, "thresh", threshMess, float);
  CPPEXTERN_MSG1(classPtr, "top", headMess, float);
  CPPEXTERN_MSG1(classPtr, "left", leftMess, float);
  CPPEXTERN_MSG1(classPtr, "right", rightMess, float);
  CPPEXTERN_MSG1(classPtr, "bottom", bottomMess, float);
  CPPEXTERN_MSG1(classPtr, "border_width", borderWidthMess, float);
  CPPEXTERN_MSG1(classPtr, "contours", contoursMess, int);
}