////////////////////////////////////////////////////////////////////// 
// janpdf/PDFpageobject.h 
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
 
#ifndef __PDFPAGEOBJECT_H__
#define __PDFPAGEOBJECT_H__

#include "IntArray.h"
#include "MathMatrix.h"
#include "PDF.h"
#include "StringBasics.h"
#include "Constant.h"
// PDFPageobject.h
// Written by Jan Wigginton

#define PDFCHART_EDGE 0.1

class PDFPageObject
   {
   public:

      String      title, subTitle;
      PDFFonts    titleFont;
      double      titleFontSize;

      PDFPageObject();
      virtual ~PDFPageObject() {};

      void Draw(PDF & pdf);

      virtual void DrawInBox(PDF & pdf, double x_0, double y_0, double x_1,
        double y_1, bool close = true);
      virtual void DrawInGrid(PDF & pdf, int row_to_use, int col_to_use,
        int rows, int cols, double spacer = 0.05);

      void DrawInUpperLeft(PDF & pdf, double bump = 0.0);
      void DrawInUpperRight(PDF & pdf, double bump = 0.0);
      void DrawInLowerLeft(PDF & pdf, double bump = 0.0);
      void DrawInLowerRight(PDF & pdf, double bump = 0.0);

   protected:

      double      height, width;
      double      titleHeight;
      double      x0, y0, x1, y1;
      double      x0Page, y0Page;

      void DrawTitle(PDF & pdf);
      void DrawOutline(PDF & pdf);
      virtual void DrawBody(PDF & pdf)= 0;

      virtual void Close(PDF & pdf);
      virtual bool Open(PDF & pdf) = 0;
      void Init();

      virtual void InitializePage(PDF & pdf);
      void SetDimensions(PDF & pdf);
      virtual void SetComponentDimensions(PDF & pdf) = 0;

      double Space();
   };

#endif
 
 
 
 
 
