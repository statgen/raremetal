////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartobject.h 
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
 
#ifndef __CHARTOBJECT_H__
#define __CHARTOBJECT_H__

#include "PDFpage.h"

class PDFChartBar;
class PDFChartLine;

class PDFChartObject
   {
   friend class PDFLineChart;
  friend class PDFmanhattan;
  friend class PDFLineChartwithPolygon;
   friend class PDFHistogram;
   friend class PDFChartBasics;
   friend class PDFChartLegend;

   public:

      PDFChartObject();
      virtual ~PDFChartObject();

      // settings

      void SetGray(double line_gray);
      void SetColor(double object_red, double object_green, double object_blue);

      void SelectColorBrush(PDF & pdf);
      void SelectGrayBrush(PDF & pdf);
      virtual void SelectColorPen(PDF & pdf);
      virtual void SelectGrayPen(PDF & pdf);

   protected:

      double  red, green, blue;
      double  gray;
      bool    isInitialized;

      void InitializeGray(int i);
      void InitializeColor(int i);
   };


#endif

 
 
 
 
 
 
 
 
 
 
