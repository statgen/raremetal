////////////////////////////////////////////////////////////////////// 
// janpdf/PDFhistogram.h 
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

#ifndef __PDFHISTOGRAM_H__
#define __PDFHISTOGRAM_H__

#include "IntArray.h"
#include "MathMatrix.h"
#include "PDF.h"
#include "PDFchartaxis.h"
#include "PDFchartlegend.h"
#include "PDFchartmarker.h"
#include "StringBasics.h"
#include "Constant.h"
#include "PDFchartbar.h"
#include "PDFchartbasics.h"

class PDFHistogram : public PDFChartBasics
   {
   public:

      bool  stacked, narrowBars;

      PDFHistogram();
      virtual ~PDFHistogram();

      // data input/chart initialization
      void ReadData(const char * infile);
      void ReadBinnedData(const char * infile);

      void SetDataCounts(const Matrix & data_counts);
      void SetCategoricalDataCounts(const Matrix & data_counts);

      // simple counts, categories ordered from 0 to data.Length() - 1
      void SetCategoricalDataCounts(const Vector & data);

      void SetDataValues(int series, const Vector & v);
      void SetDataValues(int series, const IntArray & v);
      void SetDataValue(int series, int data_pt, double value);
      void SetDataValues(const Matrix & values);

      void Dimension(int num_series, int data_pts_per_series);
      void Reset();

      // chart set tings/preferences
      void SetBarWidth(double bar_width);
      void SetBarStart(double bar_start);
      void SetBarEnd(double bar_end);
      void SetBarSpacing(double bar_spacing);

      void SetDataValueTag(int series, double value, const char * tag);

      void SetUseSeriesPercentages(bool use_series_pct);
      void SetUseTotalPercentages(bool use_total_pct);
      void SetCategoricalBars(bool categorical, double first_bar_category = _NAN_,
                                 double last_bar_category = _NAN_);

   private:

      // dataValues[i] holds the ith column of data. Matrix is private because
      // changes to data require recalculating the bars to update the chart.
      Matrix   dataValues;

      PDFChartBar ** bars;

      // used for binned input data to initialize bar ranges
      Matrix   xRanges;

      double   barStart, barEnd;
      double   barWidth, barSpacing, barSpacingPct;
      double   dataMin, dataMax;
      int      numBars;
      bool     keepBarWidth, keepBarSpacing;
      bool     useSeriesPercentages, useTotalPercentages;

      int      barSeries;

      bool tooManyBars;

      // flag for new data input
      bool dirty;

      //flag for data read in as counts instead of raw data
      bool useBinData;

      //values on histogram to tag, strings to use as tags
      Vector      taggedValues;
      StringArray valueTags;

      // drawing subroutines
      virtual void DrawLegend(PDF & pdf);
      virtual void DrawBody(PDF & pdf);

      void DrawBarsStacked(PDF & pdf);
      void DrawBarsAligned(PDF & pdf);
      void DrawNoDataLabel(PDF & pdf, double x_pos);
      void DrawErrorGraph(PDF & pdf);

      // chart initialization, reinitialization for new data
      virtual bool Open(PDF & pdf);
      bool CheckForNoData();

      void InitializeForNewData();
      void InitializeBars();
      void InitializeBarsFromBins();
      void InitializeAxes();
      void Init();

      void AllocateBarMatrix(int num_series, int num_bars_per_series);
      void ResizeBarMatrix(int old_num_bars = 0, int old_num_series = 0);
      void AssignBarTags();

      int    CalculateBarRanges(int times_tried);
      int    CalculateBarWidth(int times_tried);
      double CalculateBarCounts(const Matrix & sorted);
      double CalculateBarCountsFromBins(const Matrix & data_values);
      void   ScaleToSeriesPct(int data_series);
      void   ScaleToTotalPct();

   };

#endif





 
 
 
 
 
