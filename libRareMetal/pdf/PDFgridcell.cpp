////////////////////////////////////////////////////////////////////// 
// janpdf/PDFgridcell.cpp 
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

#include "PDFgridcell.h"
#include "Constant.h"


PDFGridCell::PDFGridCell()
   {
   xPos = yPos = _NAN_;
   height = width = _NAN_;

   borderRed = borderGreen = borderBlue = 0.0;
   textRed = textGreen = textBlue = 0.0;
   red = green = blue = 1.0;

   numeric_value = _NAN_;
   }

PDFGridCell::~PDFGridCell()
   {
   }

void PDFGridCell::operator=(const PDFGridCell & rhs)
   {
   red = rhs.red;
   green = rhs.green;
   blue = rhs.blue;
   gray = rhs.gray;
   isInitialized = rhs.isInitialized;

   xPos = rhs.xPos;
   yPos = rhs.yPos;
   height = rhs.height;
   width = rhs.width;

   numeric_value = rhs.numeric_value;
   text_value = rhs.text_value;

   textRed = rhs.textRed;
   textGreen = rhs.textGreen;
   textBlue = rhs.textBlue;

   borderRed = rhs.borderRed;
   borderGreen = rhs.borderGreen;
   borderBlue = rhs.borderBlue;
   }

void PDFGridCell::SetTextColor(double text_red, double text_green, double text_blue)
   {
   textRed = text_red;
   textGreen = text_green;
   textBlue = text_blue;
   }

void PDFGridCell::SetBorderColor(double border_red, double border_green, double border_blue)
   {
   borderRed = border_red;
   borderGreen = border_green;
   borderBlue = border_blue;
   }

void PDFGridCell::SelectBorderBrush(PDF & pdf)
   {
   pdf.page.SetFillColor(borderRed, borderGreen, borderBlue);
   }

void PDFGridCell::SelectTextBrush(PDF & pdf)
   {
   pdf.page.SetFillColor(textRed, textGreen, textBlue);
   }




 
 
 
 
 
 
 
 
