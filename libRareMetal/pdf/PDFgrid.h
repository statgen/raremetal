////////////////////////////////////////////////////////////////////// 
// janpdf/PDFgrid.h 
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

#ifndef __PDFGRID_H__
#define __PDFGRID_H__

#include "MathMatrix.h"
#include "PDF.h"
#include "PDFchartbasics.h"
#include "PDFgridcell.h"

class PDFGrid : public PDFChartBasics
   {
   public:

      PDFGridCell ** cells;

      Vector rowHeightPcts;
      Vector colWidthPcts;

      bool printLower, printUpper;

      PDFGrid();
      virtual ~PDFGrid();

      void Dimension(int num_rows, int num_cols);

      void SetRowHeightPct(int row, double height);
      void SetColWidthPct(int row, double width);

      void SetCellText(int row, int col, const char * text);
      void SetCellText(int row, int col, double val, int precision = -1);

      void SetCellColor(int row, int col, double red, double green, double blue);
      void SetRowColor(int row, double red, double green, double blue);
      void SetColumnColor(int col, double red, double green, double blue);
      void SetCellBorderColor(int row, int col, double red, double green, double blue);

   private:

      int  nRows, nCols;
      bool keepRowHeights, keepColWidths;

      void CheckRowColSettings();

      virtual bool Open(PDF & pdf);
      virtual void DrawBody(PDF & pdf);
      virtual void DrawLegend(PDF & pdf);

      void AllocateCellMatrix(int num_rows, int num_cols);
      void ResizeCellMatrix(int num_rows, int num_cols);
      void DimensionRowSettings(int num_rows);
      void DimensionColSettings(int num_cols);

      void InitializeAxes(int num_rows, int num_cols);
      void InitializeCellPositions();

      void SetComponentDimensions(PDF & pdf);
   };

#endif





 
 
 
 
 
