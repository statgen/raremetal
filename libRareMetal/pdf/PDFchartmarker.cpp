////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartmarker.cpp 
// (c) 2000-2007 Goncalo Abecasis (c) 2002-2007 Jan Wigginton
// 
// This file is distributed as part of the PEDSTATS source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile PEDSTATS.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Tuesday December 18, 2007
// 
 
/* Written by Jan Wigginton */

#include "PDFchartmarker.h"
#include "Constant.h"
#include "Error.h"

PDFDataMarker::PDFDataMarker()
   {
   radius = _NAN_;
   shape = lmDiamond;
   red = blue = green = gray = 0.0;
   }


PDFDataMarker::~PDFDataMarker()
   {
   }

void PDFDataMarker::DrawAt(PDF & pdf, double x_0, double y_0)
   {
   double r = radius;
   double r_half = radius * 0.5;

   pdf.page.SetFillColor(red,green,blue);
   switch(shape)
      {
      case lmCircle:
         pdf.page.FillCircle(x_0, y_0, r);
         break;
      case lmSquare:
         pdf.page.FillRectangle(x_0 - r_half, y_0 - r_half, x_0 + r_half,  \
         y_0 + r_half);
         break;
      case lmCross:
         pdf.page.DrawLine(x_0 - 0.75 * r, y_0 - 0.75 * r, x_0 + 0.75 * r,   \
         y_0 + 0.75 *r);
         pdf.page.DrawLine(x_0  - 0.75 * r, y_0 + 0.75 * r, x_0 + 0.75 * r,  \
         y_0 - 0.75 * r);
         break;
      case lmDiamond:
         {
         double x_coord[4] = { x_0 - r_half, x_0, x_0 + r_half, x_0 };
         double y_coord[4] = { y_0, y_0 + r_half, y_0, y_0 - r_half };
         pdf.page.DrawPolygon(x_coord, y_coord, 4);
         }
         break;
      default:
         error("PDFDataMarker: Illegal marker type specified");
         break;
      }
   }


void PDFDataMarker::operator=(const PDFDataMarker & rhs)
   {
   radius = rhs.radius;
   shape = rhs.shape;
   red = rhs.red;
   blue = rhs.blue;
   green = rhs.green;
   }

 
 
 
 
 
 
 
 
 
 
