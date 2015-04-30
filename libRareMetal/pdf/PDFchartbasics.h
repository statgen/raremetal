////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartbasics.h 
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

#ifndef __PDFCHARTBASICS_H__
#define __PDFCHARTBASICS_H__

#include "IntArray.h"
#include "MathMatrix.h"
#include "PDF.h"
#include "PDFchartaxis.h"
#include "PDFchartlegend.h"
// PDFChartBasics.h
// Written by Jan Wigginton

#include "PDFchartmarker.h"
#include "PDFpageobject.h"
#include "StringBasics.h"
#include "Constant.h"

#define PDFCHART_EDGE 0.1

class PDFChartBasics : public PDFPageObject
   {
   public:

      PDFChartAxis   xAxis, yAxis;
      PDFChartLegend legend;

      bool        useLegend, useColor, useHeader;
      bool        drawHGrid, drawXY, drawQuadrants;
      bool        drawVConnector;
      bool        symmetricAxes;

      //Alias for legend.lines pointer
      PDFChartLine * & lines;

      PDFChartBasics();
      virtual ~PDFChartBasics();

      void SetSeriesColor(int series, double red, double green, double blue);
      void SetSeriesGray(int series, double gray);
      void SetSeriesLabel(int series, const char * label);

      virtual void DrawInBox(PDF & pdf, double x_0, double y_0, double x_1, double y_1, bool close = true);

   protected:

      double      xScale, yScale;
      bool        noData;
      bool        keepLegend;

      // Drawing subroutines...
      virtual void DrawBody(PDF & pdf)= 0;

      void DrawAxis(PDF & pdf,  PDFChartAxis & axis);
      virtual void  DrawTicks(PDF & pdf,  PDFChartAxis & axis);
      virtual void DrawTickLabels(PDF & pdf, PDFChartAxis & axis);
      void DrawAxisLabels(PDF & pdf, PDFChartAxis & axis);
      void DrawTitle(PDF & pdf);
      virtual void DrawLegend(PDF & pdf) = 0;

      // Initialization subroutines
      void InitializePage(PDF & pdf);
      void InitializeBlankGraph(PDF & pdf);

      virtual void Close(PDF & pdf);
      virtual bool Open(PDF & pdf) = 0 ;
      void Reset();
      void Init();

      int  ReadData(const char *filename, Matrix &temp_values);

      // Scaling, sizing subroutines
      void SetDimensions(PDF & pdf);
      virtual void SetComponentDimensions(PDF & pdf);

      void SetTickSize(PDFChartAxis & axis);
      void SetAxisDimensions(PDF & pdf, PDFChartAxis & axis);
      void SetTitleDimensions(PDF & pdf);
      void SetLegendDimensions(PDF & pdf);
      void SetAxisScaling();

      // Utility functions
      double MapX(double pos);
      double MapY(double pos);
   };

#endif






