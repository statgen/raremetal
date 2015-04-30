////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartbasics.cpp 
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

#include "PDFchartbasics.h"
#include "MathMatrix.h"
#include "Constant.h"
#include "MathConstant.h"
#include "Parameters.h"
#include "StringArray.h"
#include "PDF.h"
#include "PDFchartaxis.h"
#include "PDFchartline.h"

#include <math.h>
#include <string.h>

#define MAPX(x) (( (x) - xAxis.min) * xScale + x0)
#define MAPY(x) (( (x) - yAxis.min) * yScale + y0)


PDFChartBasics::PDFChartBasics() : xAxis(false), yAxis(true), lines(legend.lines)
   {
   noData = false;
   Init();
   }

void PDFChartBasics::Init()
   {
   useLegend = keepLegend = true;
   useColor = useHeader = true;
   symmetricAxes = false;
   drawHGrid = drawXY = drawQuadrants = false;
   drawVConnector = false;
   }

PDFChartBasics::~PDFChartBasics()
   {
   }

void PDFChartBasics::Reset()
  {
  xAxis.Reset();
  yAxis.Reset();
  legend.Reset();
  Init();
  }

int PDFChartBasics::ReadData(const char *filename, Matrix &temp_values)
   {
   StringArray input, tokens;
   char * flag;
   double value;
   int pts_read = 0;

   noData = true;

   input.Read(filename);

   for (int j = 0; j < input.Length(); j++)
      {
      tokens.Clear();
      tokens.AddTokens(input[j]);

       if (legend.labels.Length() == 0 && useLegend)
         {
         if (useHeader)
            {
            legend.labels = tokens;
            continue;
            }
         }

      if (tokens.Length() == 0)
         continue;

      if (temp_values.rows == 0)
         temp_values.Dimension(tokens.Length(), 0);

      pts_read++;

      if (tokens.Length() != temp_values.rows)
        error("PDFchartbasics::Unlabeled missing data points at line %d of %s\n", j, filename);

      for (int i = 0; i < temp_values.rows; i++)
         {
         value = strtod((const char *) tokens[i], &flag);
         if (*flag)
            value = _NAN_;
         else
            noData = false;
         temp_values[i].Push(value);
         }
      }

   temp_values.Dimension(temp_values.rows, pts_read);
   if (legend.labels.Length() < temp_values.rows - 1)
      legend.InitializeLabels(temp_values.rows - 1);

   return pts_read;
   }


//////Chart scaling, sizing subroutines/////////////////////////////////////////


void PDFChartBasics::SetComponentDimensions(PDF & pdf)
   {
   bool done = false;
   bool x0_done = false;
   bool y0_done = false;
   bool y1_done = false;

   while (!done)
      {
      if (useLegend) SetLegendDimensions(pdf);
      x1 = x0Page + width - legend.width - PDFCHART_EDGE * width;

      SetAxisDimensions(pdf, xAxis);
      SetAxisDimensions(pdf, yAxis);

      if (yAxis.width > 0.5 * (x0 - x0Page)  && !x0_done)
         {
         x0 = x0Page + yAxis.width + 0.5 * PDFCHART_EDGE * width;
         x0_done = true;
         continue;
         }

      if (xAxis.height > 0.5 * (y0 - y0Page) && !y0_done)
         {
         y0 = y0Page + xAxis.height + 0.5 * PDFCHART_EDGE * height;
         y0_done = true;
         continue;
         }

      SetTitleDimensions(pdf);

      if (titleHeight > 0.5 * (y0Page + height - y1) && !y1_done)
         {
         y1_done = true;
         y1 = y0Page + height;
         y1 -= (titleHeight + PDFCHART_EDGE * height);
         continue;
         }

      done = true;
      }
   SetAxisScaling();
   }

void PDFChartBasics::SetTitleDimensions(PDF & pdf)
   {
   titleFontSize = (y1 - y0) * 0.08;
   pdf.page.SetFont(titleFont);
   pdf.page.SetFontSize(titleFontSize);

   int count = 0;
   double width_title;

   while (true)
      {
      width_title = pdf.page.TextWidth(title);

      if (width_title < (x1 - x0))
         break;

      titleFontSize *= 0.9;
      pdf.page.SetFontSize(titleFontSize);
      count++;

      if (count > 10)
         break;
      }
   titleHeight = pdf.page.TextHeight(title) * 2.0;
   }

void PDFChartBasics::SetAxisDimensions(PDF & pdf, PDFChartAxis & axis)
   {
   char label[1000];

   SetTickSize(axis);
   if (!axis.keepFontSize) axis.tickLabelFontSize = axis.tickSize * 2.5;
   axis.labelFontSize = axis.tickSize * 2.5;

   pdf.page.SetFont(axis.font);
   pdf.page.SetFontSize(axis.tickLabelFontSize);

   if (axis.vert)
      {
      double start = (axis.stepStart != _NAN_ ? axis.stepStart : axis.min);

      double max_pos = axis.FloorEPS(axis.max / axis.step) * axis.step + axis.Eps();
      double pos = axis.CeilEPS(start / axis.step ) * axis.step;

      double width_curr, label_width_max = 0.0;

      for (; pos <= max_pos; pos += axis.step)
         {
         if (fabs(pos) < 1e-6 && axis.step > 1e-5)
            sprintf(label, "%.*f", axis.digits, fabs(pos));
         else
            sprintf(label, axis.useScientific ? "%.*e" : "%.*f" , axis.digits, pos);

         width_curr = pdf.page.TextWidth(label);
         label_width_max = max(width_curr, label_width_max);
         }

      double min_label_width = pdf.page.TextWidth("99.99");
      axis.tickLabelSize = max( min_label_width, label_width_max );

      pdf.page.SetFontSize(axis.tickLabelFontSize);
      pdf.page.SetFont(axis.font, true);

      axis.width = axis.tickLabelSize + axis.tickSize + Space();
      if (axis.label.Length() > 0)
         axis.width += (pdf.page.TextHeight(axis.label) + 4.0 * Space());

      axis.width = max(axis.width, PDFCHART_EDGE * width);
      axis.width = min(axis.width, 0.5 * width);
      axis.height = y1 - y0;
      }
   else
      {
      sprintf(label, "%.*f", axis.digits, axis.min);
      axis.tickLabelSize = pdf.page.TextHeight(label);

      pdf.page.SetFontSize(axis.labelFontSize);
      pdf.page.SetFont(axis.font, true);

      axis.height = axis.tickSize + axis.tickLabelSize + Space();
      if (axis.label.Length() > 0)
         axis.height  += (pdf.page.TextHeight(axis.label) + 4.0 * Space());

      axis.height = max( axis.height, PDFCHART_EDGE * height);
      axis.height = min(axis.height,  0.5 * height);
      axis.width = x1 - x0;
      }
   }

void PDFChartBasics::SetTickSize(PDFChartAxis & axis)
   {
   double wide = x1 - x0;
   double high = y1 - y0;
   axis.tickSize = min(wide * 0.02, high * 0.02);
   }

void PDFChartBasics::SetLegendDimensions(PDF & pdf)
   {
   legend.fontSize = 2.0 * Space();
   pdf.page.SetFont(legend.font);
   pdf.page.SetFontSize(legend.fontSize);

   legend.CalculateWidth(pdf, width);
   legend.CalculateHeights(pdf,height);
   legend.CalculateOrigin(pdf, y0, x0Page, width);
   }

void PDFChartBasics::SetAxisScaling(void)
   {
   if (yAxis.max == yAxis.min)
      error("PDFChartBasics: Y-axis maximum and minimum are same");

   yScale =  (y1 - y0) / (yAxis.max - yAxis.min);
   xScale = xAxis.max != xAxis.min ? (x1 - x0) / (xAxis.max - xAxis.min) : 1.0;
   }

///////Chart drawing subroutines////////////////////////////////////////////////

void PDFChartBasics::DrawInBox(PDF & pdf, double x_0, double y_0, double x_1, double y_1, bool close)
   {
   PDFPageObject::DrawInBox(pdf, x_0, y_0, x_1, y_1, false);

   DrawOutline(pdf);
   pdf.page.SetClipRectangle(x0Page, y0Page, x1 + 0.5 * ( x0Page + width - x1), y0Page + height);
   DrawAxis(pdf, xAxis);
   pdf.page.SetClipRectangle(x0Page, y0Page, x0Page + width, y1 + 0.5 * ( y0Page + height - y1));
   DrawAxis(pdf, yAxis);
   if (useLegend && keepLegend) DrawLegend(pdf);
   if (close)
      Close(pdf);
   }

void PDFChartBasics::DrawAxis(PDF & pdf, PDFChartAxis & axis)
   {
   pdf.page.hTextAlignment = axis.vert ? taRight : taCenter;
   pdf.page.vTextAlignment = axis.vert ? taMiddle : taBelow;

   pdf.page.SetFontSize(axis.tickLabelFontSize);
   pdf.page.SetFont(axis.font);

   if (axis.useTicks)
      DrawTicks(pdf, axis);

   if (axis.useNumericLabels || axis.useStringLabels )
      DrawTickLabels(pdf,axis);

   pdf.page.SetFont(axis.font, true);
   DrawAxisLabels(pdf, axis);
   }

void PDFChartBasics::DrawTicks(PDF & pdf, PDFChartAxis & axis)
   {
   if (noData) return;

   double start = (axis.stepStart != _NAN_ ? axis.stepStart: axis.min);

   double max_pos = axis.FloorEPS(axis.max / axis.step)* axis.step + axis.Eps();
   double pos = axis.CeilEPS(start / axis.step) * axis.step;

   for (; pos <= max_pos; pos += axis.step)
      if (axis.vert)
         {
         pdf.page.DrawLine(x0, MapY(pos), x0 - axis.tickSize, MapY(pos));
         if (drawHGrid) pdf.page.DrawLine(x0, MapY(pos), x1, MapY(pos));
         }
      else
         if (pos >=  axis.min ) pdf.page.DrawLine(MapX(pos), y0, MapX(pos), y0 - axis.tickSize);
   }

void PDFChartBasics::DrawTickLabels(PDF & pdf, PDFChartAxis & axis)
   {
   if (noData) return;

   char label[100];
   String string_label;

   double start = (axis.stepStart != _NAN_ ? axis.stepStart : axis.min);

   double max_pos = axis.FloorEPS(axis.max / axis.step) * axis.step + axis.Eps();
   double pos = axis.CeilEPS(start / axis.step ) * axis.step;

   if (axis.labelsOffset) pos += axis.step * 0.5;

   int i = 0;

   for (; pos <= max_pos; pos += axis.step)
      {
      if (axis.useNumericLabels)
         {
         if (fabs(pos) < 1e-6 && axis.step > 1e-5)
            sprintf(label, "%.*f", axis.digits, fabs(pos));
         else
            sprintf(label, axis.useScientific ? "%.*e" : "%.*f" , axis.digits, pos);

         string_label = label;
         }
      else
         string_label = axis.useStringLabels && i < axis.stringTickLabels.Length() ?
           (const char *) axis.stringTickLabels[i] : "";

      double red_fill = axis.tickLabelHighlights.Length() > i && axis.tickLabelHighlights[i] == 1 ? 1.0 : 0.0;
      pdf.page.SetFillColor(red_fill, 0.0, 0.0);
      i++;

      pdf.page.hTextAlignment = axis.vert ? taRight : taCenter;
      pdf.page.vTextAlignment = axis.vert ? taMiddle : taBelow;

      if (axis.vert)
         {
         pdf.page.WriteText(x0 - axis.tickSize - Space(), MapY(pos), (const char *) string_label);
         continue;
         }

      if (pos >= axis.min)
         {
         double mult = i%2 == 1 || !axis.useStringLabels || !axis.useAlternatingLabels ? 1.0 : 4.0;
         pdf.page.WriteText(MapX(pos), y0 - axis.tickSize - mult * Space(), (const char *) string_label);
         }
      }
   }

void PDFChartBasics::DrawAxisLabels(PDF & pdf, PDFChartAxis & axis)
   {
   pdf.page.hTextAlignment = axis.vert ? taRight : taCenter;
   pdf.page.vTextAlignment = axis.vert ? taMiddle : taBelow;

   if (axis.vert)
      {
      pdf.page.SetFontOrientation(90.0);
      pdf.page.WriteText(x0 - axis.tickSize - axis.tickLabelSize - 4.0 * Space(),
        0.5 * (y0 + y1) , axis.label);
      pdf.page.SetFontOrientation(0.0);
      }
   else
      pdf.page.WriteText((x0 + x1) * 0.5,
        y0 - axis.tickSize - axis.tickLabelSize - 4.0 * Space(), axis.label);
   }

void PDFChartBasics::Close(PDF & pdf)
   {
   pdf.page.ClearClipRectangle();
   }

void PDFChartBasics::InitializeBlankGraph(PDF & pdf)
   {
   if (noData)
      {
      if (!xAxis.keepMin)
         xAxis.min = 0.0;
      if (!xAxis.keepMax)
         xAxis.max = (xAxis.keepMin ? xAxis.min + 1.0 : 1.0);

      if (!yAxis.keepMin)
         yAxis.min = 0.0;
      if (!yAxis.keepMax)
         yAxis.max = (yAxis.keepMax ? xAxis.min + 1.0 : 1.0);

      xAxis.step = yAxis.step = 1.0;
      }
   }

void PDFChartBasics::InitializePage(PDF & pdf)
   {
   PDFPageObject::InitializePage(pdf);

   legend.width = legend.height = _NAN_;
   legend.x0 = legend.y0 = _NAN_;
   }

////Chart preferences, settings ///////////////////////////////////////////////

void PDFChartBasics::SetSeriesColor(int series, double red, double green, double blue)
   {
   if (!lines)
      error("PDFLineChart: No data loaded, cannot specify line colors");

   if (series >= 0 && series < legend.numSeries)
      {
      lines[series].SetColor(red, green, blue);
      lines[series].SetMarkerColor(red, green, blue);
      }
   }

void PDFChartBasics::SetSeriesGray(int series, double gray)
   {
   if (!lines)
      error("PDFChart: No data loaded, cannot specify line colors");

   if (series >= 0 && series < legend.numSeries)
      lines[series].SetGray(gray);
   }

void PDFChartBasics::SetSeriesLabel(int series, const char * label)
   {
   if (!lines)
      error("PDFChartBasics: No legend label array");

   if (series >= 0 && series < legend.numSeries)
      legend.labels[series] = label;
   }


///////Utility functions////////////////////////////////////////////////////////

double PDFChartBasics::MapX(double pos)
   {
   return ((pos - xAxis.min) * xScale + x0);
   }

double PDFChartBasics::MapY(double pos)
   {
   return ((pos - yAxis.min) * yScale + y0);
   }













 
 
 
 
 
