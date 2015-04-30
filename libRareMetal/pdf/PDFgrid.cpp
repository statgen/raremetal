////////////////////////////////////////////////////////////////////// 
// janpdf/PDFgrid.cpp 
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

#include "PDFgrid.h"
#include "PDFgridcell.h"

PDFGrid::PDFGrid()
   {
   cells = NULL;
   nRows = 0;
   nCols = 0;
   keepRowHeights = keepColWidths = false;
   useLegend = false;
   printLower = printUpper = true;
   }

PDFGrid::~PDFGrid()
   {
   if (cells)
      {
      for (int i = 0; i < nRows; i++)
         if (cells[i])
            delete [] cells[i];
      delete [] cells;
      cells = NULL;
      }
   }

void PDFGrid::SetRowHeightPct(int row, double height_pct)
   {
   if (height_pct < 0.0 || height_pct > 100.0)
      error("PDFGrid::Height percentage given must be between 0.0 and 100.0");

   rowHeightPcts[row] = height_pct;
   }

void PDFGrid::SetColWidthPct(int col, double width_pct)
   {
   if (width_pct < 0.0 || width_pct > 100.0)
      error("PDFGrid::Width percentage given must be between 0.0 and 100.0");

   colWidthPcts[col] = width_pct;
   }

void PDFGrid::Dimension(int num_rows, int num_cols)
   {
   if (num_rows <= 0)
      error("PDFGrid::Number of rows must be a positive integer");
   if (num_cols <= 0)
      error("PDFGrid::Number of columns must be a positive integer");

   DimensionRowSettings(num_rows);
   DimensionColSettings(num_cols);
   AllocateCellMatrix(num_rows, num_cols);
   }

void PDFGrid::DimensionRowSettings(int num_rows)
   {
   if (num_rows == nRows)
      return;

   if (!keepRowHeights)
      {
      rowHeightPcts.Clear();
      rowHeightPcts.Dimension(num_rows);
      for (int i = 0; i < num_rows; i++)
         rowHeightPcts[i] =  100.0/ (double) num_rows;
      }
   }

void PDFGrid::DimensionColSettings(int num_cols)
   {
   if (num_cols == nCols)
      return;

   if (!keepColWidths)
      {
      colWidthPcts.Clear();
      colWidthPcts.Dimension(num_cols);
      for (int i = 0; i < num_cols; i++)
         colWidthPcts[i] =  100.0/(double) num_cols;
      }
   }

void PDFGrid::AllocateCellMatrix(int num_rows, int num_cols)
   {
   if (!cells && num_rows > 0 && num_cols > 0)
      {
      cells = new PDFGridCell *[num_rows];
      for (int i = 0; i < num_rows; i++)
           cells[i] = new PDFGridCell[num_cols];
      }
   else
      ResizeCellMatrix(num_rows, num_cols);

   nRows = num_rows;
   nCols = num_cols;
   }

void PDFGrid::ResizeCellMatrix(int new_rows, int new_cols)
   {
   PDFGridCell ** new_cells = NULL;

   if (new_rows > 0)
      new_cells = new PDFGridCell *[new_cols];

   for (int i = 0; i < new_rows; i++)
      if (new_cols > 0)
         new_cells[i] = new PDFGridCell[new_cols];
      else
         new_cells[i] = NULL;

   for (int i = 0; i < new_rows; i++)
      if (i < nRows)
         for (int j = 0; j < new_cols; j++)
            if ( j < nCols)
               new_cells[i][j] = cells[i][j];

   for (int i = 0; i < nRows; i++)
       if (cells[i]) delete [] cells[i];

   if (cells) delete [] cells;

   cells = new_cells;
   }

bool PDFGrid::Open(PDF & pdf)
   {
   CheckRowColSettings();
   InitializeAxes(nRows, nCols);
   InitializeCellPositions();

   return true;
   }

void PDFGrid::SetComponentDimensions(PDF & pdf)
   {
   PDFChartBasics::SetComponentDimensions(pdf);

   if ( nCols > 20 || nRows > 20 )
      {
      double x_font, y_font;

      x_font = MapX(1.5) - MapX(0.5);
      y_font = MapY(1.5) - MapY(0.5);

      double axis_font = ::min(x_font, y_font);

      xAxis.SetTickLabelFontSize(axis_font);
      yAxis.SetTickLabelFontSize(axis_font);
      }
   }

void PDFGrid::CheckRowColSettings()
   {
   double sum = 0.0;
   sum = rowHeightPcts.Sum();
   if (sum > 100.01 || sum < 99.99)
      {
      if (keepRowHeights)
         error("PDFGrid::Row height percentages do not add to l00%");
      else
         DimensionRowSettings(nRows);
      }

   sum = colWidthPcts.Sum();
   if (sum > 100.01 || sum < 99.99 )
      {
      if (keepColWidths)
         error("PDFGrid::Column width percentages do not add to l00%");
      else
         DimensionColSettings(nCols);
       }
   }

void PDFGrid::InitializeAxes(int num_rows, int num_cols)
   {
   xAxis.SetMin(0.0);
   xAxis.SetMax(nCols);
   yAxis.SetMin(0.0);
   yAxis.SetMax(nRows);
   yAxis.SetStepStart(-0.5);
   xAxis.SetStepStart(-0.5);

   xAxis.labelsOffset = yAxis.labelsOffset = true;
   xAxis.useDiscreteValues = yAxis.useDiscreteValues = true;
   xAxis.useTicks = yAxis.useTicks = false;

   xAxis.SetStringTickLabels(true);
   yAxis.SetStringTickLabels(true);

//   yAxis.SetTickLabelFontSize(8.0 * Space() * (5.0 / (num_cols + 5.0)));
//   xAxis.SetTickLabelFontSize(8.0 * Space() * (5.0 / (num_cols + 5.0)));

   String label;
   int n_xlabels = xAxis.stringTickLabels.Length();
   for (int i = num_cols - 1; i >= n_xlabels ; i--)
      {
      label = i;
      xAxis.SetTickLabel(i, (const char *) label);
      }

   int n_ylabels = yAxis.stringTickLabels.Length();
   for (int i = num_rows - 1; i >= n_ylabels ; i--)
      {
      label = i;
      yAxis.SetTickLabel(i, (const char *) label);
      }

   xAxis.SetStep(1.0);
   yAxis.SetStep(1.0);
   }

void PDFGrid::InitializeCellPositions()
   {
   for (int i = 0; i < nRows; i++)
     for (int j = 0; j < nCols; j++)
        {
        cells[i][j].xPos = 0.5 + i;
        cells[i][j].yPos = 0.5 + j;
        cells[i][j].width = colWidthPcts[j] * (xAxis.max - xAxis.min) * 0.01;
        cells[i][j].height = rowHeightPcts[i] * (yAxis.max - yAxis.min) * 0.01;
        }
   }

void PDFGrid::DrawLegend(PDF & pdf)
   {
   }

void PDFGrid::DrawBody(PDF & pdf)
   {
   pdf.page.SetLineWidth(0.5 * Space());
   pdf.page.vTextAlignment = taMiddle;
   pdf.page.hTextAlignment = taCenter;

   double rect_x0, rect_x1;
   double rect_y0, rect_y1;

   for (int i = 0; i < nRows; i++)
      for (int j = 0; j < nCols; j++)
         {
         if (!printUpper && j > i) continue;
         if (!printLower && j < i) continue;

         useColor ? cells[i][j].SelectColorBrush(pdf) : cells[i][j].SelectGrayBrush(pdf);

         rect_x0 = cells[i][j].xPos + 100 * xAxis.Eps() - cells[i][j].width * 0.5;
         rect_x1 = cells[i][j].xPos + 100 * xAxis.Eps() + cells[i][j].width * 0.5;
         rect_y0 = cells[i][j].yPos + 100 * yAxis.Eps() - cells[i][j].height * 0.5;
         rect_y1 = cells[i][j].yPos + 100 * yAxis.Eps() + cells[i][j].height * 0.5;

         pdf.page.FillRectangle(MapX(rect_x0), MapY(rect_y0), MapX(rect_x1), MapY(rect_y1));

         cells[i][j].SelectBorderBrush(pdf);
         pdf.page.DrawRectangle(MapX(rect_x0), MapY(rect_y0), MapX(rect_x1), MapY(rect_y1));

         if (cells[i][j].text_value != "")
            {
            cells[i][j].SelectTextBrush(pdf);
            pdf.page.SetFontSize(6.0 * Space() * ( 5.0 / (nCols + 5.0)));
            pdf.page.WriteText(MapX(cells[i][j].xPos), MapY(cells[i][j].yPos),
                              (const char * )cells[i][j].text_value);
            }
         }
   }

void PDFGrid::SetCellText(int row, int col, const char * text)
   {
   if (row < 0 || row > nRows)
      error("PDFGrid::Row argument is out of bounds in set cell data");
   if (col < 0 || col > nCols)
      error("PDFGrid::Column argument is out of bounds in set cell data");

   cells[row][col].text_value = text;
   cells[row][col].numeric_value = _NAN_;
   }

void PDFGrid::SetCellText(int row, int col, double val, int precision)
   {
   int digits = precision;

   if (row < 0 || row > nRows)
      error("PDFGrid::Row argument is out of bounds in set cell data");
   if (col < 0 || col > nCols)
      error("PDFGrid::Column argument is out of bounds in set cell data");

   if (digits == -1)
      digits = fabs(val) > 0.1 || val == 0.0 ? 1 : 4;

   cells[row][col].text_value.printf("%.*f", digits, val);
   cells[row][col].numeric_value = val;
   }

void PDFGrid::SetRowColor(int row, double red, double green, double blue)
   {
   if (row < 0 || row > nRows)
      error("PDFGrid::Attempt to set row color for row %d - this row number is out of range", row);

   for (int j = 0; j < nCols; j++)
      cells[row][j].SetColor(red, green, blue);
   }

void PDFGrid::SetColumnColor(int col, double red, double green, double blue)
   {
   if (col < 0 || col > nCols)
      error("PDFGrid::Attempt to set column color for column %d - this column number is out of range", col);

   for (int i = 0; i < nRows; i++)
      cells[i][col].SetColor(red, green, blue);
   }

void PDFGrid::SetCellColor(int row, int col, double red, double green, double blue)
   {
   if (col < 0 || col > nCols)
      error("PDFGrid::Attempt to set column color for column %d - this column number is out of range", col);

   if (row < 0 || row > nRows)
      error("PDFGrid::Attempt to set row color for row %d - this row number is out of range", row);

   cells[row][col].SetColor(red, green, blue);
   }

void PDFGrid::SetCellBorderColor(int row, int col, double red, double green, double blue)
   {
   if (col < 0 || col > nCols)
      error("PDFGrid::Attempt to set cell border color for row %d column %d - this column number is out of range", row, col);

   if (row < 0 || row > nRows)
      error("PDFGrid::Attempt to set cell border color for row %d col %d - this row number is out of range", row, col);

   cells[row][col].SetBorderColor(red, green, blue);
   }



 
 
 
 
 
