////////////////////////////////////////////////////////////////////// 
// WritePDF.h 
// (c) 2012-2013 Shuang Feng, Dajiang Liu, Goncalo Abecasis
// 
// This file is distributed as part of the RareMetalWorker source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile RareMetalWorker.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Wednesday November 28, 2012
// 
 
#ifndef __WRITEPDF_H__
#define __WRITEPDF_H__

#include "PDF.h"
#include "PDFlinechart.h"
#include "PDFlinechartwithpolygon.h"
#include "PDFmanhattan.h"

class WritePDF
{
   public:
      WritePDF();
      static int resolution;
      static bool thinPoints;
      //Draw QQ plots of three series
      void Draw(PDF & pdf,StringArray & geneLabels,Vector & pvalueAll, Vector & pvalue1,Vector & pvalue5, StringArray & chr, Vector & pos,String title, String demo1,String demo2,String demo3,bool useDataLabels);
      void Draw(PDF & pdf,StringArray & geneLabels, Vector & pvalueAll, StringArray & chr, Vector & pos,String title,String extraTitle,String demo,bool useDataLabels);
      void DrawOverlayQQplot(PDF & pdf,Vector & pvalueAll,Vector & pvalue1,Vector & pvalue5,String title);
      //Draw QQ plots of single series
      void DrawQQplot(PDF & pdf, Vector & pvalueAll,String title, String plot_title,String demo);
      void DrawManhattanPlot(PDF & pdf, StringArray & geneLabels,Vector & pvalueAll, StringArray & chr, Vector & pos,String title,bool useDataLabels);
};

#endif


