////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartline.h 
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
 
#ifndef CHARTLINE_H
#define CHARTLINE_H

#include "PDFpage.h"
#include "PDFchartmarker.h"
#include "PDFchartobject.h"

class PDFChartLine : public PDFChartObject
  {
  friend class PDFLineChart;
  friend class PDFLineChartwithPolygon;
  friend class PDFmanhattan;
  friend class PDFChartBasics;
  friend class PDFHistogram;
  friend class PDFChartLegend;

  public:

     PDFChartLine();
     virtual ~PDFChartLine();

  private:

     PDFDataMarker   marker;
     double       line_weight;
     PDFLineStyle line_style;

      bool    hasMarkers;
      bool    showLine;
      bool    skipSeries;

     void operator=(const PDFChartLine & rhs);

     void SetMarkerColor(double marker_red, double marker_green, double marker_blue);
     void SetMarkerShape(PDFMarkerShapes marker_shape);
     void SetMarkerRadius(double marker_radius);

     void SetShowLine(bool show_line);
     void SetShowMarker(bool show_marker);

     void SetStyle(PDFLineStyle style);
     void SetWeight(double weight);

     void InitializeMarker(int i);
     void InitializeStyle(int i);

     virtual void SelectColorPen(PDF & pdf);
     virtual void SelectGrayPen(PDF & pdf);

     virtual void DrawMarkerAt(PDF & pdf, double x, double y) { marker.DrawAt(pdf, x, y); }

  };

#endif
 
 
 
 
 
 
 
 
 
 
