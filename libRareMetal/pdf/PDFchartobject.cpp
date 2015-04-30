////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartobject.cpp 
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
 
/* This file written by Jan Wigginton */
 
#include "PDFchartobject.h"
#include "PDF.h"
#include "Constant.h"

PDFChartObject::PDFChartObject()
  {
  red = green = blue = 0.0;
  gray = 1.0;
  isInitialized = false;
  }

PDFChartObject::~PDFChartObject()
  {
  }

void PDFChartObject::InitializeGray(int i)
  {
  gray = (0.5 * (((i) % 6) / 3));
  }

void PDFChartObject::InitializeColor(int i)
  {
  red = (i % 7) & 4;
  green = (i % 7) & 2;
  blue = (i % 7) & 1;
  }

void PDFChartObject::SelectColorBrush(PDF & pdf)
   {
   pdf.page.SetFillColor(red, green, blue);
   }

void PDFChartObject::SelectGrayBrush(PDF & pdf)
   {
   pdf.page.SetFillGray(gray);
   }

void PDFChartObject::SelectColorPen(PDF & pdf)
  {
  pdf.page.SetLineColor(red, green, blue);
  pdf.page.SetLineStyle(lsSolid);
  pdf.page.SetLineWidth(1.0);
  }

void PDFChartObject::SelectGrayPen(PDF & pdf)
  {
  pdf.page.SetLineGray(gray);
  pdf.page.SetLineStyle(lsSolid);
  pdf.page.SetLineWidth(1.0);
  }

void PDFChartObject::SetColor(double object_red, double object_green, double object_blue)
  {
  red = object_red;
  green = object_green;
  blue = object_blue;
  }

void PDFChartObject::SetGray(double object_gray)
  {
  gray = object_gray;
  }

 
 
 
 
 
 
 
 
 
 
