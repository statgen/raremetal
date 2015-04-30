////////////////////////////////////////////////////////////////////// 
// janpdf/PDFlinechartwithpolygon.cpp 
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

#include "MathMatrix.h"
#include "Constant.h"
#include "MathConstant.h"
#include "Parameters.h"
#include "StringArray.h"
#include "PDF.h"
#include "PDFlinechartwithpolygon.h"
#include "PDFchartline.h"
#include "QuickIndex.h"

#include <math.h>
#include <string.h>


PDFLineChartwithPolygon::PDFLineChartwithPolygon() : PDFChartBasics()
   {
   yValues.Dimension(0, 0);
   xValues.Dimension(0, 0);
   }

PDFLineChartwithPolygon::~PDFLineChartwithPolygon()
   {
   }

/// Data input/ chart initialization subroutines////////////////////////////////

void PDFLineChartwithPolygon::ReadData(const char * filename)
   {
   Matrix temp_values;

   PDFChartBasics::ReadData(filename, temp_values);

   if (!noData)
      {
      Dimension(temp_values.rows - 1, temp_values.cols);
      for (int i = 0; i < temp_values.rows; i++)
         if (temp_values[i][0] == _NAN_)
            {
            warning("PDFLineChartwithPolygon::Non-numeric x-input value read at line %d %lf,"
             "Some data points won't be drawn", i, temp_values[i][0]);

            for (int j = 0; j < temp_values.cols; j++)
               {
               for (int k = 0; k < temp_values.rows; k++)
                  printf("%lg  ", temp_values[k][j]);
               printf("\n");
               }
            break;
            }

      SetDataValues(temp_values);
      }
   else
      Dimension(0, 0);
   }


void PDFLineChartwithPolygon::SetDataValues(Vector & x_values, Matrix & y_values)
   {
   // check for (correct) prior call to Dimension
   if (y_values.rows != yValues.rows || y_values.cols != yValues.cols)
      error("PDFLineChartwithPolygon: Cannot set data values. Dimension of y matrix argument"
            "differs from chart data matrix" );

   if (x_values.dim != xValues.cols)
      error("PDFLineChartwithPolygon: Cannot set data values. Dimension of x vector argument"
            "differs from chart data vector" );

   //if we have data, check for missing data in x cooordinates
   if (!noData)
      for (int i = 0; i < x_values.dim; i++)
         if (x_values[i] == _NAN_)
            {
            warning("PDFLineChartwithPolygon::Missing x-value argument supplied, some data points "
              "won't be drawn");
            break;
            }
    // sort all data based on x_values
    QuickIndex index(x_values);

    for (int j = 0; j < x_values.dim; j++)
       {
       int out = index[j];
       xValues[0][j] = x_values[out];
       for (int i = 0; i < y_values.rows; i++)
          yValues[i][j] = y_values[i][out];
       }

    // mark data series with all data missing
    for (int i = 0; i < yValues.rows; i++)
      {
      lines[i].skipSeries = true;
      for (int j = 0 ; j < xValues.cols; j++)
         if (yValues[i][j] != _NAN_)
            {
            lines[i].skipSeries = false;
            break;
            }
      }
   }

void PDFLineChartwithPolygon::SetDataValues(Matrix & input_values)
   {
   Vector x_values;

   if (input_values.rows <= 0)
      {
      noData = true;
      Dimension(0, 0);
      return;
      }

   int cols = input_values[0].Length();
   x_values.Dimension(cols);
   for (int i = 0; i < cols; i++)
      x_values[i] = input_values[0][i];

   input_values.DeleteRow(0);

   SetDataValues(x_values, input_values);
   }

void PDFLineChartwithPolygon::SetDataValue(int series, int datapt, double value)
   {
   if (series < 0 || series >= yValues.rows)
      error("PDFLineChartwithPolygon: Cannot set data value. Series argument is out of range");

   if (datapt >= yValues.cols)
      error("PDFLineChartwithPolygon::Cannot set data value. Data point is out of range");

   yValues[series][datapt] = value;
   SetDataValues(xValues[0], yValues);
   }

void PDFLineChartwithPolygon::Reset()
   {
   PDFChartBasics::Reset();
   legend.Reset();
   legend.Initialize(yValues.rows);
   }

void PDFLineChartwithPolygon::Dimension(int series, int columns)
   {
   int old_series = yValues.rows;
   int old_cols = yValues.cols;

   xValues.Dimension(1, columns);
   yValues.Dimension(series, columns);

   for (int i = 0; i < yValues.rows; i++)
      for (int j = old_series; j < yValues.cols; j++)
            yValues[i][j] = _NAN_;
   for (int i = old_series; i < yValues.rows; i++)
      for (int j = 0; j < old_cols; j++)
         if (j < yValues.cols) yValues[i][j] = _NAN_;

   legend.Initialize(series);
   for (int i = 0; i < legend.numSeries; i++)
      {
      lines[i].SetWeight(1.0);
      lines[i].hasMarkers = true;
      lines[i].showLine = true;
      }

   noData = (series == 0 || columns == 0);
   }

/// Private subroutines called after Draw() ////////////////////////////////////

bool PDFLineChartwithPolygon::Open(PDF & pdf)
   {
   InitializePage(pdf);
   keepLegend = true;

   // if we think we've got some data, make sure this is the case
   if (!noData)
      noData = CheckForNoData();

   if (noData)
      {
      keepLegend = false;
      InitializeBlankGraph(pdf);
      return false;
      }

   xAxis.GetDataRange(xValues);
   yAxis.GetDataRange(yValues);

   if (symmetricAxes)
     BalanceAxes();

   // Styles specific to line charts
   for (int i = 0; i < legend.numSeries; i++)
      useColor ?  lines[i].InitializeStyle(0) : lines[i].InitializeStyle(i);

   return true;
   }

void PDFLineChartwithPolygon::BalanceAxes()
   {
   double set_max = ::max(xAxis.max, yAxis.max);
   double set_min = ::min(xAxis.min, yAxis.min);
   xAxis.SetMax(set_max);
   yAxis.SetMax(set_max);
   xAxis.SetMin(set_min);
   yAxis.SetMin(set_min);

   xAxis.GetDataRange(xValues);
   yAxis.GetDataRange(yValues);
   }

bool PDFLineChartwithPolygon::CheckForNoData()
   {
   for (int i = 0; i < yValues.rows; i++)
      for (int j = 0; j < yValues.cols; j++)
         if (yValues[i][j] != _NAN_)
            return false;

   return true;
   }

void PDFLineChartwithPolygon::SetComponentDimensions(PDF & pdf)
   {
   PDFChartBasics::SetComponentDimensions(pdf);

   for (int i = 0; i < yValues.rows; i++)
     {
     if (lines[i].marker.radius == _NAN_)
        SetMarkerRadius(i, 0.5 * Space());
     SetLineWeight(i, 0.5 * Space());
     }
   }
void PDFLineChartwithPolygon::DrawDemo(PDF & pdf)
{
if(demo=="") return;
pdf.page.SetFontSize((y1-y0)*0.04);
pdf.page.WriteText(MapX(xAxis.min+(xAxis.max-xAxis.min)*0.20),MapY(yAxis.max-(yAxis.max-yAxis.min)*0.05),(const char *)demo);
}

void PDFLineChartwithPolygon::DrawBody(PDF & pdf)
   {
   DrawQuadrants(pdf);
   //Add polygon here
   pdf.page.SetFillColor(0.66,0.66,0.66);
   for(int i=0;i<poly_n;i++)
   {
      poly_x[i] = MapX(poly_x[i]);
      poly_y[i] = MapY(poly_y[i]);
   }
DrawDemo(pdf);
   pdf.page.FillPolygon(poly_x,poly_y,poly_n);

   pdf.page.SetLineWidth(0.5 * Space());

   for (int i = 0; i < yValues.rows; i++)
      {
      useColor ? lines[i].SelectColorPen(pdf) : lines[i].SelectGrayPen(pdf);

      // if we use this series, draw markers at data points
      if (lines[i].hasMarkers && !lines[i].skipSeries)
         for (int j = 0; j < xValues.cols; j++)
            if (yValues[i][j] != _NAN_ && xValues[0][j] != _NAN_)
              lines[i].marker.DrawAt(pdf, MapX(xValues[0][j]), MapY(yValues[i][j]));

      // if we're connecting data points with line
      if (lines[i].showLine)
         {
         // move to first non-missing data point
         int j = 0;
         for ( ; j < xValues.cols; j++)
            if (yValues[i][j] != _NAN_ && xValues[0][j] != _NAN_)
               {
               pdf.page.PathMoveTo(MapX(xValues[0][j]), MapY(yValues[i][j]));
               break;
               }

         // then draw the line
         for ( ; j < xValues.cols; j++)
            if (yValues[i][j] != _NAN_ && xValues[0][j] != _NAN_)
               pdf.page.PathLineTo(MapX(xValues[0][j]),MapY(yValues[i][j]));

         pdf.page.PathStroke();
         }
      }

   if (drawVConnector)
      {
      double connect_min = _NAN_;
      double connect_max = _NAN_;

      pdf.page.SetLineColor(0.0, 0.0, 0.0);

      for (int j = 0; j < xValues.cols; j++)
         {
         if (xValues[0][j] != _NAN_)
            {
            int i = 0;
            for ( ; i < yValues.cols; i++)
               if (yValues[i][j] != _NAN_ && !lines[i].skipSeries)
                  {
                  connect_min = connect_max = yValues[i][j];
                  break;
                  }

            for (; i < yValues.rows; i++)
               if (yValues[i][j] != _NAN_ && !lines[i].skipSeries)
                  {
                  connect_min = min(yValues[i][j], connect_min);
                  connect_max = max(yValues[i][j], connect_max);
                  }

            if (connect_min != _NAN_ || connect_max != _NAN_)
               pdf.page.DrawLine(MapX(xValues[0][j]), MapY(connect_min), MapX(xValues[0][j]), MapY(connect_max));
            connect_min = connect_max = _NAN_;
            }
         }
      }
   }

void PDFLineChartwithPolygon::DrawLegend(PDF & pdf)
   {
   legend.Draw(pdf, Space(), useColor);
   }

void PDFLineChartwithPolygon::DrawQuadrants(PDF & pdf)
   {
   pdf.page.SetLineColor(0.0, 0.0, 0.0);
   pdf.page.SetLineWidth(0.25 * Space());

   if (yAxis.includesZero && drawQuadrants)
      pdf.page.DrawLine(x0, MapY(0.0), x1, MapY(0.0));

   if (xAxis.includesZero && drawQuadrants)
      pdf.page.DrawLine(MapX(0.0), y0, MapX(0.0), y1);

   if (drawXY)
      pdf.page.DrawLine(x0, y0, x1, y1);
   }


void PDFLineChartwithPolygon::Close(PDF & pdf)
   {
   PDFChartBasics::Close(pdf);

   for (int i = 0; i < yValues.rows; i++)
      lines[i].marker.radius = _NAN_;
   }

///////////////// Chart preferences/settings///////////////////////////////////

void PDFLineChartwithPolygon::SetLineStyle(int series, PDFLineStyle style)
   {
   if (!lines)
      error("PDFLineChartwithPolygon: No data loaded, cannot specify line colors");

   if (series >= 0 && series < legend.numSeries)
      lines[series].SetStyle(style);
   }

void PDFLineChartwithPolygon::SetLineWeight(int series, double line_weight)
   {
   if (!lines)
      error("PDFLineChartwithPolygon: No data loaded, cannot specify line weight");

   if (series >= 0 && series < legend.numSeries)
      lines[series].SetWeight(line_weight);
   }

void PDFLineChartwithPolygon::SetMarkerShape(int series, PDFMarkerShapes shape)
   {
   if (!lines)
      error("PDFLineChartwithPolygon: No data loaded, cannot specify line markers");

   if (series >= 0 && series < legend.numSeries)
      lines[series].SetMarkerShape(shape);
   }

void PDFLineChartwithPolygon::SetMarkerBlack(int series, bool is_black)
   {
   if (!lines)
      error("PDFLineChartwithPolygon: No data loaded, cannot specify line markers");

   if (series >= 0 && series < legend.numSeries)
      {
      if (is_black)
         lines[series].SetMarkerColor(0.0, 0.0, 0.0);
      else
         lines[series].SetMarkerColor(lines[series].red,
           lines[series].green, lines[series].blue);
      }
   }

void PDFLineChartwithPolygon::SetMarkerColor(int series, double red, double green, double blue)
   {
   if (!lines)
      error("PDFLineChartwithPolygon: No data loaded, cannot specify marker color");

   if (series >= 0 && series < legend.numSeries)
      lines[series].SetMarkerColor(red, green, blue);
   }

void PDFLineChartwithPolygon::SetMarkerRadius(int series, double radius)
   {
   if (!lines)
      error("PDFLineChartwithPolygon: No data loaded, cannot specify line colors");

   if (series >= 0 && series < legend.numSeries)
      lines[series].SetMarkerRadius(radius);
   }

void PDFLineChartwithPolygon::ShowMarker(int series, bool show_marker)
   {
   if (!lines)
      error("PDFLineChartwithPolygon: No data loaded, cannot specify marker use");

   if (series >= 0 && series < legend.numSeries)
      lines[series].SetShowMarker(show_marker);
   }

void PDFLineChartwithPolygon::ShowLine(int series, bool show_line)
   {
   if (!lines)
      error("PDFLineChartwithPolygon: No data loaded, cannot specify line properties");

   if (series >= 0 && series < legend.numSeries)
      lines[series].SetShowLine(show_line);
   }

void PDFLineChartwithPolygon::ShowMarkers(bool show_markers)
   {
   for (int i = 0; i < legend.numSeries; i++)
      ShowMarker(i,show_markers);
   }

void PDFLineChartwithPolygon::ShowLines(bool show_lines)
   {
   for (int i = 0; i < legend.numSeries; i++)
      ShowLine(i, show_lines);
   }


 
 
 
 
 
 
 
 
