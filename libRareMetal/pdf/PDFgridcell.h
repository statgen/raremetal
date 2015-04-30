////////////////////////////////////////////////////////////////////// 
// janpdf/PDFgridcell.h 
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

#ifndef __PDFGRIDCELL_H__
// PDFgridcell.h
// Written by Jan Wigginton

#define __PDFGRIDCELL_H__

#include "PDF.h"
#include "PDFchartobject.h"

class PDFGridCell : public PDFChartObject
   {
   friend class PDFGrid;

   public:

      PDFGridCell();
      ~PDFGridCell();

      void operator=(const PDFGridCell & rhs);
      void SelectBorderBrush(PDF & pdf);
      void SelectTextBrush(PDF & pdf);

   private:

      double xPos, yPos;
      double height, width;

      String text_value;
      double numeric_value;

      double textRed, textGreen, textBlue;
      double borderRed, borderGreen, borderBlue;

      void SetTextColor(double text_red, double text_green, double text_blue);
      void SetBorderColor(double border_red, double border_green, double border_blue);

     // double numText;
     //  bool useText;
   };


#endif


 
 
 
 
 
 
 
 
