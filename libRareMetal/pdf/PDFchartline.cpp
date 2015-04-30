////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartline.cpp 
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

#include "PDFchartline.h"
#include "PDF.h"
#include "Constant.h"

PDFChartLine::PDFChartLine()
   {
   showLine = true;
   hasMarkers = true;
   skipSeries = false;
   line_style = lsSolid;
   line_weight = 1.0;
   }

PDFChartLine::~PDFChartLine()
   {
   }

void PDFChartLine::operator=(const PDFChartLine & rhs)
   {
   red = rhs.red;
   green = rhs.green;
   blue = rhs.blue;
   gray = rhs.gray;
   isInitialized = rhs.isInitialized;

   line_style = rhs.line_style;
   line_weight = rhs.line_weight;
   hasMarkers = rhs.hasMarkers;
   showLine = rhs.showLine;
   skipSeries = rhs.skipSeries;
   marker = rhs.marker;
   }

void PDFChartLine::InitializeMarker(int i)
   {
   marker.red = (((i % 7) & 4) == 0);
   marker.green =(((i % 7) & 2) == 0);
   marker.blue = (((i % 7) & 1) == 0);
   }

void PDFChartLine::SetMarkerColor(double marker_red, double marker_green, double marker_blue)
   {
   marker.red = marker_red;
   marker.green = marker_green;
   marker.blue = marker_blue;
   }

void PDFChartLine::SetMarkerShape(PDFMarkerShapes marker_shape)
   {
   marker.shape = marker_shape;
   }

void PDFChartLine::SetMarkerRadius(double marker_radius)
   {
   marker.radius = marker_radius;
   }

void PDFChartLine::SetShowLine(bool show_line)
   {
   showLine = show_line;
   }

void PDFChartLine::SetShowMarker(bool show_marker)
   {
   hasMarkers = show_marker;
   }

void PDFChartLine::InitializeStyle(int i)
   {
   PDFLineStyle styles[3] = { lsSolid, lsDashed, lsDotted };
   line_style = styles[(i) % 3];
   }

void PDFChartLine::SetStyle(PDFLineStyle style)
   {
   line_style = style;
   }

void PDFChartLine::SetWeight(double weight)
   {
   line_weight = weight;
   }

void PDFChartLine::SelectColorPen(PDF & pdf)
   {
   pdf.page.SetLineColor(red, green, blue);
   pdf.page.SetLineStyle(line_style);
   if (line_weight != _NAN_) pdf.page.SetLineWidth(line_weight);
   }

void PDFChartLine::SelectGrayPen(PDF & pdf)
   {
   pdf.page.SetLineGray(gray);
   pdf.page.SetLineStyle(line_style);
   if (line_weight != _NAN_) pdf.page.SetLineWidth(line_weight);
   }




 
 
 
 
 
 
