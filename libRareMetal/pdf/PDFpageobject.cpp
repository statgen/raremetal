////////////////////////////////////////////////////////////////////// 
// janpdf/PDFpageobject.cpp 
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
 
#include "PDFpageobject.h"

// Written by Jan Wigginton

PDFPageObject::PDFPageObject()
   {
   x0 = x1 = y0 = y1 = _NAN_;
   width = height = _NAN_;
   Init();
   }

void PDFPageObject::Init()
   {
   titleFont = fHelvetica;
   title = "";
   subTitle = "";
   titleFontSize = 0.0;
   titleHeight = 0.0;
   }

void PDFPageObject::SetDimensions(PDF & pdf)
   {
   if (x0 != _NAN_)
      {
      x0Page = x0;
      y0Page = y0;

      width = x1 - x0;
      height = y1 - y0;
      }
   else
      {
      x0Page = 0.0;
      y0Page = 0.0;

      height = pdf.page.GetHeight();
      width = pdf.page.GetWidth();
      }

   x0 = x0Page + width * PDFCHART_EDGE;
   y0 = y0Page + height * PDFCHART_EDGE;
   x1 = x0Page + width * (1.0 - PDFCHART_EDGE);
   y1 = y0Page + height * (1.0 - PDFCHART_EDGE);

   SetComponentDimensions(pdf);
   }

void PDFPageObject::Draw(PDF & pdf)
   {
   pdf.page.OpenPage();
   DrawInBox(pdf, _NAN_, _NAN_, _NAN_, _NAN_);
   pdf.page.ClosePage();
   }

void PDFPageObject::DrawInUpperLeft(PDF & pdf, double bump)
   {
   double page_width = pdf.page.GetWidth();
   double page_height = pdf.page.GetHeight();

   double x_0 = page_width * 0.0;
   double y_0 = page_height * (0.45 + bump);
   double x_1 = page_width * 0.5;
   double y_1 = page_height * (0.9 + bump);

   DrawInBox(pdf, x_0, y_0, x_1, y_1);
   }

void PDFPageObject::DrawInUpperRight(PDF & pdf, double bump)
   {
   double page_width = pdf.page.GetWidth();
   double page_height = pdf.page.GetHeight();

   double x_0 = page_width * 0.5;
   double y_0 = page_height * (0.45 + bump);
   double x_1 = page_width * 1.0;
   double y_1 = page_height * (0.9 + bump);

   DrawInBox(pdf, x_0, y_0, x_1, y_1);
   }

void PDFPageObject::DrawInLowerLeft(PDF & pdf, double bump)
   {
   double page_width = pdf.page.GetWidth();
   double page_height = pdf.page.GetHeight();

   double x_0 = page_width * 0.0;
   double y_0 = page_height * (0.0 + bump);
   double x_1 = page_width * 0.5;
   double y_1 = page_height * (0.45 + bump);

   DrawInBox(pdf, x_0, y_0, x_1, y_1);
   }

void PDFPageObject::DrawInLowerRight(PDF & pdf, double bump)
   {
   double page_width = pdf.page.GetWidth();
   double page_height = pdf.page.GetHeight();

   double x_0 = page_width * 0.5;
   double y_0 = page_height * (0.0 + bump);
   double x_1 = page_width * 1.0;
   double y_1 = page_height * (0.45 + bump);

   DrawInBox(pdf, x_0, y_0, x_1, y_1);
   }

void PDFPageObject::DrawInGrid(PDF & pdf, int row_to_use, int col_to_use, int rows, int cols,
double spacer)
   {
   double page_width = pdf.page.GetWidth();
   double page_height = pdf.page.GetHeight();
   //bool legend_setting = useLegend;

   double x_0 = 0.99 * page_width * (double(col_to_use)/cols);
   double y_0 = 0.9 * page_height * (double(row_to_use)/rows + spacer);
   double x_1 = 0.99 * page_width * ((col_to_use + 1.0)/ cols );
   double y_1 = 0.9 * page_height * ((row_to_use + 1.0)/rows) - spacer;

   //useLegend = false;
   DrawInBox(pdf, x_0, y_0, x_1, y_1);
   //useLegend = legend_setting;
   }
/*
void PDFPageObject::DrawInBox(PDF & pdf, double x_0, double y_0, double x_1, double y_1)
   {
   bool draw_body = OpenChart(pdf);

   x0 = x_0;
   x1 = x_1;
   y0 = y_0;
   y1 = y_1;

   SetChartDimensions(pdf);

   pdf.page.SetClipRectangle(x0Page, y0Page,x0Page + width, y0Page + height);
   if (!title.IsEmpty()) DrawTitle(pdf);
   DrawOutline(pdf);
   pdf.page.SetClipRectangle(x0Page, y0Page, x1 + 0.5 * ( x0Page + width - x1), y0Page + height);
   DrawAxis(pdf, xAxis);
   pdf.page.SetClipRectangle(x0Page, y0Page, x0Page + width, y1 + 0.5 * ( y0Page + height - y1));
   DrawAxis(pdf, yAxis);
   pdf.page.SetClipRectangle(x0, y0, x1, y1);
   if (draw_body) DrawBody(pdf);
   pdf.page.ClearClipRectangle();
   if (useLegend && keepLegend) DrawLegend(pdf);
   CloseChart(pdf);
   }
  */

void PDFPageObject::DrawInBox(PDF & pdf, double x_0, double y_0, double x_1, double y_1, bool close)
   {
   bool draw_body = Open(pdf);

   x0 = x_0;
   x1 = x_1;
   y0 = y_0;
   y1 = y_1;

   SetDimensions(pdf);

   pdf.page.SetClipRectangle(x0Page, y0Page,x0Page + width, y0Page + height);
   if (!title.IsEmpty()) DrawTitle(pdf);

   pdf.page.SetClipRectangle(x0, y0, x1, y1);
   if (draw_body) DrawBody(pdf);
   pdf.page.ClearClipRectangle();

   if (close)
      Close(pdf);
   }

void PDFPageObject::DrawTitle(PDF & pdf)
   {
   pdf.page.SetFont(titleFont, true);
   pdf.page.SetFontSize(titleFontSize);
   pdf.page.hTextAlignment = taCenter;

   double title_height = pdf.page.TextHeight(title);

   pdf.page.vTextAlignment = taMiddle;
   pdf.page.WriteText(0.5 * (x0 + x1), 0.5 * (y0Page + height + y1), (const char *) title);

   if (subTitle.Length() > 0)
      {
      pdf.page.SetFontSize(titleFontSize * 0.6);
      pdf.page.vTextAlignment = taBelow;
      pdf.page.WriteText(0.5 * (x0 + x1),
                         0.5 *( y0Page + height + y1) - 1.0 * title_height,
                         (const char *) subTitle);
      }
   }

void PDFPageObject::DrawOutline(PDF & pdf)
   {
   pdf.page.SetLineColor(0.0, 0.0, 0.0);
   pdf.page.SetLineWidth(0.25 * Space());
   pdf.page.DrawRectangle(x0, y0, x1, y1);
   }

void PDFPageObject::InitializePage(PDF & pdf)
   {
   pdf.page.SetFontOrientation(0.0);

   x0 = x1 = y0 = y1 = _NAN_;
   x0Page = y0Page = _NAN_;
   }

void PDFPageObject::Close(PDF & pdf)
   {
   pdf.page.ClearClipRectangle();
   }

double PDFPageObject::Space(void)
   {
   return 0.005 * (x1 - x0 + y1 - y0);
   }


 
 
 
 
 
