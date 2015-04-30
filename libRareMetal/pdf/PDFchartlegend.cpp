////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartlegend.cpp 
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

#include "PDFchartlegend.h"
#include "Constant.h"
#include "Error.h"
#include "PDFchartline.h"
#include "PDFchartbar.h"


PDFChartLegend::PDFChartLegend()
   {
   width = height = 0.0;
   rowHeight = labelWidth = 0.0;
   x0 = y0 = _NAN_;
   fontSize = 0.0;
   font = fHelvetica;
   widthMult = 2.6;

   lines = NULL;
   numSeries = 0;
   }

PDFChartLegend::~PDFChartLegend()
   {
   if (lines)
      delete [] lines;
   }

void PDFChartLegend::CalculateWidth(PDF & pdf, double page_width)
   {
   double curr_width;

   pdf.page.SetFontSize(fontSize);
   pdf.page.SetFont(font);

   double max_width = pdf.page.TextWidth(labels[0]);

   for (int i = 0; i < labels.Length(); i++)
      {
      curr_width = pdf.page.TextWidth(labels[i]);
      if (curr_width > max_width) max_width = curr_width;
      }
   labelWidth = max_width;

   width =  widthMult * max_width;
   }

void PDFChartLegend::CalculateHeights(PDF & pdf, double page_height)
   {
   pdf.page.SetFontSize(fontSize);
   pdf.page.SetFont(font);

   double label_height = pdf.page.TextHeight(labels[0]);

   rowHeight = label_height * 2.0;
   if ((height = rowHeight * (labels.Length() + 1)) < page_height) return;

   while (true)
      {
      // first try reducing row height
      rowHeight = label_height * 1.5;
      if ( (height = rowHeight * (labels.Length() + 1)) < page_height) break;

      // then try reducing fontSize
      if (fontSize > LEGEND_MIN_FONT)
         {
         fontSize *= 0.8;
         label_height = pdf.page.TextHeight(labels[0]);
         }

      // just cut off values - won't fit
      else
         {
         height = page_height;
         break;
         }
      }
   }

void PDFChartLegend::CalculateOrigin(PDF & pdf, double chart_y0, double page_x0, double page_width)
   {
   y0 = chart_y0;
   x0 = page_x0 + page_width - width - 0.05 * page_width;
   }

void PDFChartLegend::Draw(PDF & pdf, double space, bool use_color)
   {
   // if there aren't any data series don't draw legend
   if (numSeries == 0) return;

   pdf.page.SetFont(font,true);
   pdf.page.SetFontSize(fontSize);
   pdf.page.hTextAlignment = taRight;
   pdf.page.vTextAlignment = taMiddle;

   pdf.page.SetLineWidth(0.25 * space);
   pdf.page.SetLineColor(0.0,0.0, 0.0);

   double y1 = y0 + height;
   double x1 = x0 + width;
   pdf.page.DrawRectangle(x0, y0, x1, y1);

   pdf.page.SetClipRectangle(x0, y0, x1, y1);
   pdf.page.SetLineWidth(0.5 * space);

   int j = 0;
   for (int i = 0; i < numSeries ; i++)
      {
      if (lines[i].skipSeries) continue;

      double y_pos =  y1 - ( (j + 1) * rowHeight);

      pdf.page.SetFillColor(0.0, 0.0, 0.0);
      pdf.page.WriteText(x1 - labelWidth * 1.4, y_pos, labels[i]);

      if (use_color)
         lines[i].SelectColorPen(pdf);
      else
         lines[i].SelectGrayPen(pdf);

      // if line isn't hidden, draw it
      if (lines[i].showLine)
         pdf.page.DrawLine(x1 - labelWidth * 1.2, y_pos, x1 - labelWidth * 0.2, y_pos);

      // otherwise, draw extra markers to symbolize series
      else if (lines[i].hasMarkers)
         {
         lines[i].DrawMarkerAt(pdf, x1 - labelWidth * 1.1, y_pos);
         lines[i].DrawMarkerAt(pdf, x1 - labelWidth * 0.3, y_pos);
         }

      if (lines[i].hasMarkers)
         lines[i].DrawMarkerAt(pdf, x1 - labelWidth * 0.7, y_pos);

      j++;
      }

   pdf.page.ClearClipRectangle();
   }

void PDFChartLegend::Reset()
  {
  for (int i = 0; i < numSeries; i++)
     lines[i].isInitialized = false;

  labels.Clear();
  Initialize(numSeries);
  }

void PDFChartLegend::Initialize(int num_series)
  {
  InitializeLabels(num_series);
  InitializeLines(num_series);
  }

void PDFChartLegend::InitializeLabels(int num_series)
   {
   int num_labels = labels.Length();
   if (numSeries == 0)
      num_labels = 0;

   if (num_series > 0)
      labels.Dimension(num_series);
   else
      {
      labels.Dimension(1);
      labels[0] = " ";
      }

   if (num_labels < num_series)
      for (int j = num_labels; j < num_series; j++)
         labels[j].printf("Series_%d", j);
   }

void PDFChartLegend::InitializeLines(int num_series)
   {
   if (!lines && num_series > 0)
     lines = new PDFChartLine[num_series];
   else if (numSeries != num_series)
      ResizeLines(num_series);

   for (int i = 0; i < num_series; i++)
      if (!lines[i].isInitialized)
         {
         lines[i].InitializeColor(i);
         lines[i].InitializeMarker(i);
         lines[i].isInitialized = true;
         }
   numSeries = num_series;
   }

void PDFChartLegend::ResizeLines(int num_series)
   {
   PDFChartLine * new_lines;
   new_lines = new PDFChartLine[num_series + (num_series == 0)];

   for (int i = 0; i < num_series; i++)
      if (i < numSeries)
         new_lines[i] = lines[i];

   if (lines)
      delete [] lines;

   lines = new_lines;
   }






 
 
 
 
 
 
