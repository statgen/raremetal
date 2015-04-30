////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartbar.cpp 
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

#include "PDFchartbar.h"
#include "PDF.h"
#include "Constant.h"

PDFChartBar::PDFChartBar()
  {
  count = 0;
  hasData = false;
  lowerBound = upperBound = _NAN_;
  barTag = "";
  isTagged = false;
  }

PDFChartBar::~PDFChartBar()
  {
  }

void PDFChartBar::operator=(const PDFChartBar & rhs)
   {
   red = rhs.red;
   green = rhs.green;
   blue = rhs.blue;
   gray = rhs.gray;
   isInitialized = rhs.isInitialized;

   hasData = rhs.hasData;
   lowerBound = rhs.lowerBound;
   upperBound = rhs.upperBound;
   count = rhs.count;
   
   isTagged = rhs.isTagged;
   barTag = rhs.barTag;
   }

void PDFChartBar::SetTag(const char * tag)
   {
   barTag = tag;
   isTagged = true;
   } 
 
 
 
 
 
 
 
