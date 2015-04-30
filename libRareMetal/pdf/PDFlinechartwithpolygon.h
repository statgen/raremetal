////////////////////////////////////////////////////////////////////// 
// janpdf/PDFlinechart.h 
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

#ifndef __PDFLINECHARTWITHPOLYGON_H__
#define __PDFLINECHARTWITHPOLYGON_H__

#include "IntArray.h"
#include "MathMatrix.h"
#include "PDF.h"
#include "PDFchartaxis.h"
#include "PDFchartlegend.h"
#include "PDFchartmarker.h"
#include "StringBasics.h"
#include "Constant.h"
#include "PDFchartline.h"
#include "PDFchartbasics.h"

class PDFLineChartwithPolygon : public PDFChartBasics
   {
   public:

      PDFLineChartwithPolygon();
      virtual ~PDFLineChartwithPolygon();

      // data input/chart initialization

      // input format for Read() is one series per column, first column contains
      // x values (i.e. x1 y11 y21 y31 y41
      //                x2 y12 y22 y32 y42...

      double * poly_x;
      double * poly_y;
      int poly_n;
String demo;
void DrawDemo(PDF & pdf);

      void ReadData(const char * infile);

      // input format for SetDataValues() is one series per row,
      // for SetDataValues(Matrix & ), first row of matrix argument contains x-values
      void SetDataValues(Matrix & input_values);
      void SetDataValues(Vector & x_vals, Matrix & y_vals);
      void SetDataValue(int series, int datapt, double value);

      void Dimension(int series, int columns);
      void Reset();

      // chart preferences/settings
      void SetLineStyle(int series, PDFLineStyle style);
      void SetLineWeight(int series, double line_weight);
      void SetMarkerShape(int series, PDFMarkerShapes shape);
      void SetMarkerBlack(int series, bool is_black);
      void SetMarkerColor(int series, double red, double green, double blue);
      void SetMarkerRadius(int series, double radius);

      void ShowMarker(int series, bool show_marker);
      void ShowLine(int series, bool show_line);
      void ShowMarkers(bool show_markers);
      void ShowLines(bool show_lines);

   private:

      // (x_i, y_ji) = (xValues[0][i], yValues[j][i]) for  j = 0,... yValues.rows
      Matrix xValues, yValues;

      // drawing subroutines...
      void DrawQuadrants(PDF & pdf);
      virtual void DrawBody(PDF & pdf);
      virtual void DrawLegend(PDF & pdf);

      bool CheckForNoData();
      void BalanceAxes();

      virtual void SetComponentDimensions(PDF & pdf);
      virtual void Close(PDF & pdf);
      virtual bool Open(PDF & pdf);
   };

#endif











