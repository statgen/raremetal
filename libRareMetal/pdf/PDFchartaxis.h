////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartaxis.h 
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
 
// Written by Jan Wigginton
 
#ifndef __CHARTAXIS_H__
#define __CHARTAXIS_H__

#include "MathConstant.h"
#include "StringArray.h"

class PDF;

#define PDF_CHARTAXIS_EPS 1e-6

class PDFChartAxis
   {
   friend class PDFLineChart;
   friend class PDFLineChartwithPolygon;
  friend class PDFmanhattan;
   friend class PDFHistogram;
   friend class PDFChartBasics;
   friend class PDFGrid;

   public:
       
      double   minMax, maxMax;
      double   minMin, maxMin;
      bool     centerData, frameData;
      bool     useTicks;
      bool     useAlternatingLabels;
      bool     useDiscreteValues;
      bool     labelsOffset;
      bool     keepFontSize;
      bool     usePow10Labels;

      double   tickSize;
      double   tickLabelSize;

      int      minDigits, maxDigits;
      double   labelFontSize, tickLabelFontSize;
      PDFFonts font;

      String   label;
      StringArray stringTickLabels;
      IntArray  tickLabelHighlights;

      PDFChartAxis(bool v = true);
      ~PDFChartAxis();

      void   SetMin(double new_min);
      void   SetMax(double new_max);
      void   SetStep(double new_step);
      void   SetStepStart(double step_start);
      void   SetStepMin(double step_min);
      void   SetTickLabel(int i, const char * tick_label);
      void   SetStringTickLabels(bool use_string);
      void   SetNumericTickLabels(bool use_numeric);
      void   SetTickLabelFontSize(double font_size);
      void   SetTickLabelHighlight(int i, bool highlight);
      void   SetMinDigits(int digits);
      void   SetMaxDigits(int digits);
      
   protected:

      double   min, max;
      double   step, stepStart, stepMin;
      bool     useNumericLabels, useStringLabels;
      bool     useScientific;

      // variables to ensure user-defined settings for step, min, and max
      // are preserved across calls to Draw()

      bool     keepStep, keepMin, keepMax;

      double   height, width;
      bool     vert, includesZero;
      int      digits;

      void GetDataRange(const Matrix & values);

      void Reset();
      void Init();

      void Center(const Matrix & values);
      void Frame(const Matrix & values);
      
      void   CalculateStep();
      void   CalculateDigits();
      void   CalculateTick();
    
      double Eps() {return PDF_CHARTAXIS_EPS * (max - min);}
      double FloorEPS(double x){ return floor(x * (1.0 + PDF_CHARTAXIS_EPS));  }
      double CeilEPS(double x) { return  ceil(x * (1.0 - PDF_CHARTAXIS_EPS)); }

    };

#endif
 
 
 
 
 
 
 
