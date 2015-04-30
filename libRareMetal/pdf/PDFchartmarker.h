////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartmarker.h 
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

#ifndef DATAMARKER_H
#define DATAMARKER_H

#include "PDF.h"

enum PDFMarkerShapes {lmCircle = 0, lmSquare = 1, lmCross = 2, lmDiamond = 3};

class PDFDataMarker
   {
   friend class PDFChartLine;
   friend class PDFLineChart;
   friend class PDFLineChartwithPolygon;
  friend class PDFmanhattan;

   private:

      PDFMarkerShapes shape;
      double  radius;
      double  red, blue, green;
      double  gray;

      PDFDataMarker();
      ~PDFDataMarker();

      void DrawAt(PDF & pdf, double x_0, double y_0);
      void operator=(const PDFDataMarker & rhs);
   };

#endif
 
 
 
 
 
 
 
 
 
 
