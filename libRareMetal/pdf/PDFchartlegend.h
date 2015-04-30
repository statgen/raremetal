////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartlegend.h 
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
 
#ifndef CHARTLEGEND_H
#define CHARTLEGEND_H

#include "StringArray.h"
#include "PDFchartobject.h"
#include "PDF.h"

#define LEGEND_MIN_FONT 4.0

class PDFChartLegend
   {
   friend class PDFLineChart;
   friend class PDFLineChartwithPolygon;
  friend class PDFmanhattan;
   friend class PDFHistogram;
   friend class PDFChartBasics;

   public:

      PDFFonts font;

      PDFChartLegend() ;
     ~PDFChartLegend();

   private:

      StringArray labels;
      PDFChartLine *lines;
      int numSeries;
      
      double width, height;
      double widthMult;
      double x0, y0;
      double rowHeight, labelWidth;
      double fontSize;

      void InitializeLabels(int num_series);
      void Initialize(int num_series);
      void InitializeLines(int num_series);
      void ResizeLines(int num_series);

      void CalculateWidth(PDF & pdf, double page_width);
      void CalculateHeights(PDF & pdf, double page_height);
      void CalculateOrigin(PDF & pdf, double chart_y0, double page_x0, double page_width);

      void Reset();
      void Draw(PDF & pdf, double space, bool use_color = true);
   };

#endif
 
 
 
 
 
 
 
 
 
 
