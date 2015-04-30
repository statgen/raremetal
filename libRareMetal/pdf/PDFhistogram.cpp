////////////////////////////////////////////////////////////////////// 
// janpdf/PDFhistogram.cpp 
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

#include "MathMatrix.h"
#include "Constant.h"
#include "MathConstant.h"
#include "Parameters.h"
#include "StringArray.h"
#include "PDF.h"
#include "PDFhistogram.h"
#include "PDFchartaxis.h"
#include "PDFchartline.h"

#include <math.h>
#include <string.h>

#define PDF_MAX_NARROW_BAR_CT 20
#define PDFHIST_MAX_BARS _NAN_

PDFHistogram::PDFHistogram()
   {
   Init();
   bars = NULL;
   numBars = barSeries = 0;
   dataMin = dataMax = _NAN_;
   yAxis.SetMin(0.0);
   }

PDFHistogram::~PDFHistogram()
   {
   if (bars != NULL)
      {
      for (int i = 0; i < barSeries; i++)
         if (bars[i] != NULL)
            delete [] bars[i];

      delete [] bars;
      bars = NULL;
      }
   }

void PDFHistogram::Init()
   {
   useBinData = false;
   barWidth = barStart = barEnd = barSpacing = barSpacingPct = _NAN_;
   keepBarWidth = keepBarSpacing = false;
   useSeriesPercentages = useTotalPercentages  = false;
   stacked = true;
   narrowBars = false;
   tooManyBars = false;
   dirty = true;
   }

/// Data input/ chart initialization routines///////////////////////////////////

void PDFHistogram::ReadData(const char * filename)
   {
   Matrix tempValues;
   int pts_per_series = PDFChartBasics::ReadData(filename, tempValues);
   Dimension(tempValues.rows, pts_per_series);
   SetDataValues(tempValues);
   }

void PDFHistogram::ReadBinnedData(const char * filename)
   {
   Matrix tempValues;
   PDFChartBasics::ReadData(filename, tempValues);
   SetDataCounts(tempValues);
   }

void PDFHistogram::Dimension(int num_series, int pts_per_series)
   {
   int old_series = dataValues.rows;
   int old_pts_per_series = dataValues.cols;

   dataValues.Dimension(num_series, pts_per_series);

   for (int i = 0; i < dataValues.rows; i++)
      for (int j = old_series; j < dataValues.cols; j++)
            dataValues[i][j] = _NAN_;

   for (int i = old_series; i < dataValues.rows; i++)
      for (int j = 0; j < old_pts_per_series; j++)
         if (j < dataValues.cols) dataValues[i][j] = _NAN_;

   legend.Initialize(num_series);

   for (int i = 0; i < legend.numSeries; i++)
      {
      lines[i].SetWeight(1.0);
      lines[i].hasMarkers = false;
      lines[i].showLine = true;
      }

   taggedValues.Dimension(num_series);
   valueTags.Dimension(num_series);
   taggedValues.Set(_NAN_);

   noData = (num_series == 0 || pts_per_series == 0);
   }

void PDFHistogram::SetDataCounts(const Matrix & data_counts)
   {
   if (data_counts.rows < 2)
      error("PDFHistogram::No data ranges given in set data counts");

   if (data_counts.rows < 3)
      error("PDFHistogram::No data values given in set data counts");

   useBinData = true;

   xRanges.Dimension(2, data_counts.cols);
   for (int i = 0; i < 2; i++)
      xRanges[i].Copy(data_counts[i]);

   Dimension(data_counts.rows - 2, data_counts.cols);

   Matrix new_counts(data_counts.rows - 2, data_counts.cols);
   for (int i = 0;  i < data_counts.rows - 2; i++)
   for (int j = 0; j < data_counts.cols; j++)
      new_counts[i][j] = data_counts[i + 2][j];

   SetDataValues(new_counts);
   }

void PDFHistogram::SetCategoricalDataCounts(const Matrix & data_counts)
   {
   if (data_counts.rows < 1)
      error("PDFHistogram::No data bins given in set categorical data counts");

   useBinData = true;

   xRanges.Dimension(2, data_counts.cols);
   for (int i = 0; i < data_counts.cols; i++)
      {
      xRanges[0][i] = data_counts[0][i];
      xRanges[1][i] = data_counts[0][i];
      }

   SetCategoricalBars(true, data_counts[0].SafeMin(), data_counts[0].SafeMax());
   Dimension(data_counts.rows - 1, data_counts.cols);

   Matrix new_counts;
   new_counts = data_counts;
   new_counts.DeleteRow(0);
   SetDataValues(new_counts);
   }

void PDFHistogram::SetCategoricalDataCounts(const Vector & counts)
   {
   Matrix categorical_data;
   categorical_data.Dimension(2, counts.Length());

   for (int i = 0; i < counts.Length(); i++)
      {
      categorical_data[0][i] = 0;
      categorical_data[1][i] = counts[i];
      }
   SetCategoricalDataCounts(categorical_data);
   }

void PDFHistogram::SetDataValues(int series, const Vector & values)
   {
   if (series < 0 || series >= dataValues.rows)
      error("PDFHistogram: Cannot set data values. Series argument is out of range");

   if (values.dim > dataValues[series].dim)
      error("PDFHistogram: Cannot set data values. Dimension of data vector greater than chart data matrix");

   for (int i = 0; i < values.dim; i++)
      dataValues[series][i] = values[i];

   for (int i = values.dim; i < dataValues[series].dim; i++)
      dataValues[series][i] = _NAN_;

   dirty = true;
   }

void PDFHistogram::SetDataValues(int series, const IntArray & values)
   {
   if (series < 0 || series >= dataValues.rows)
      error("PDFHistogram: Cannot set data values. Series argument is out of range");

   int len = values.Length();

   if (len > dataValues[series].dim)
      error("PDFHistogram: Cannot set data values. Dimension of data vector greater than chart data matrix");

   for (int i = 0; i < len; i++)
      dataValues[series][i] = values[i];

   for (int i = len; i < dataValues[series].dim; i++)
      dataValues[series][i] = _NAN_;

   dirty = true;
   }

void PDFHistogram::SetDataValue(int series, int datapt, double value)
   {
   if (series < 0 || series >= dataValues.rows)
      error("PDFHistogram: Cannot set data value. Series argument is out of range");

   if (datapt > dataValues.cols)
       error("PDFHistogram:: Cannot set data value. Datapt argument is out of range");

   dataValues[series][datapt] = value;
   dirty = true;
   }

void PDFHistogram::SetDataValues(const Matrix & values)
   {
   if (values.rows != dataValues.rows || values.cols != dataValues.cols)
      error("PDFHistogram: Cannot set data values. Dimension of matrix argument"
      "differs from chart data matrix" );

   dataValues.Copy(values);
   dirty = true;
   }

void PDFHistogram::Reset()
   {
   PDFChartBasics::Reset();
   Init();

   for (int i = 0; i < legend.numSeries; i++)
      lines[i].isInitialized = false;

   legend.labels.Clear();
   legend.InitializeLabels(legend.numSeries);
   }

/// Private subroutines called to initialize the chart////////////////////

bool PDFHistogram::Open(PDF & pdf)
   {
   if (useBinData)
      InitializeBarsFromBins();

   if (dirty)
      InitializeForNewData();

   if (noData)
      {
      InitializePage(pdf);
      AllocateBarMatrix(0, 0);
      InitializeBlankGraph(pdf);
      return false;
      }

   if (!useBinData)
      InitializeBars();

   InitializeAxes();
   InitializePage(pdf);

   return true;
   }

void PDFHistogram::InitializeForNewData()
   {
   dataMin = _NAN_;
   dataMax = _NAN_;

   if (!xAxis.keepMin) xAxis.min = _NAN_;
   if (!xAxis.keepMax) xAxis.max = _NAN_;
   if (!yAxis.keepMax) yAxis.max = _NAN_;
   if (!yAxis.keepMin) yAxis.min = _NAN_;
   if (!keepBarWidth) barWidth = _NAN_;
   if (!xAxis.keepStep) xAxis.step = _NAN_;

   if (!noData)
     noData = CheckForNoData();

   dirty = false;
   }

bool PDFHistogram::CheckForNoData()
   {
   if (!useBinData)
      {
      for (int i = 0; i < dataValues.rows; i++)
         for (int j = 0; j < dataValues.cols; j++)
            if (dataValues[i][j] != _NAN_)
               return false;
      }
   else
      {
      for (int i = 0; i < dataValues.rows; i++)
         for (int j = 0; j < dataValues.cols; j++)
            if (dataValues[i][j] != 0.0)
               return false;
      }

   return true;
   }

void PDFHistogram::InitializeAxes()
   {
   xAxis.GetDataRange(useBinData ? xRanges : dataValues);

   if (!xAxis.keepMin)
      {
      xAxis.min = min( xAxis.min, dataMin - barWidth * 0.5);
      if (xAxis.minMin != _NAN_)
         xAxis.min = max(xAxis.min, xAxis.minMin);
      }

   if (!xAxis.keepMax)
      {
      xAxis.max = max(xAxis.max, dataMax + barWidth * 0.5);
      if (xAxis.maxMax != _NAN_)
         xAxis.max = min(xAxis.max, xAxis.maxMax);
      }

   // For y-Axis, determine the range of the counts, which depends on the type
   // of histogram (stacked vs unstacked), the x-axis range settings and the data.

   Matrix count_values(stacked ? 1 : dataValues.rows, numBars);
   count_values.Zero();

   int j = 0;

   for (int k = 0; k < dataValues.rows; k++)
      {
      for (int i = 0; i < numBars; i++)
         if (bars[k][i].lowerBound + barWidth >= xAxis.min
          && bars[k][i].upperBound - barWidth <= xAxis.max )
            count_values[j][i] += bars[k][i].count;
         if (!stacked) j++;
      }

   yAxis.GetDataRange(count_values);
   if (!yAxis.keepMax)
      {
      double max = count_values.SafeMax();
      if (yAxis.max <= max + 0.1 * yAxis.step)
         yAxis.max += 0.5 * yAxis.step;
      if (yAxis.maxMax != _NAN_)
         yAxis.max = min(yAxis.max, yAxis.maxMax);
      }
   }

void PDFHistogram::InitializeBars()
   {
   int    pass = 0;
   double max_bar_count;

   Matrix sortedData;
   sortedData.Copy(dataValues);

   for (int i = 0; i < sortedData.rows; i++)
      sortedData[i].Sort();

   do
      {
      int new_num_bars = CalculateBarRanges(pass);
      AllocateBarMatrix(dataValues.rows, new_num_bars);
      max_bar_count = CalculateBarCounts(sortedData);
      pass++;
      }
   while (pass <=3 && narrowBars && max_bar_count >= PDF_MAX_NARROW_BAR_CT);

   AssignBarTags();
   }

void PDFHistogram::InitializeBarsFromBins()
   {
   // first calculate initialize bar ranges from the bin ranges given
   bool fixedMin = dataMin != _NAN_;
   bool fixedMax = dataMax != _NAN_;

   if (!fixedMin)
      dataMin = ::min(xRanges[0].SafeMin(), xRanges[1].SafeMin());

   if (!fixedMax)
      dataMax = ::max(xRanges[0].SafeMax(), xRanges[1].SafeMax());

   if (dataMax < dataMin)  { double tmp = dataMax; dataMax = dataMin; dataMin = tmp; }

   if (!fixedMin) dataMin = xAxis.FloorEPS(dataMin / barWidth) * barWidth;
   if (!fixedMax) dataMax = xAxis.CeilEPS(dataMax / barWidth) * barWidth;

   AllocateBarMatrix(dataValues.rows, dataValues.cols);
   CalculateBarCountsFromBins(dataValues);
   AssignBarTags();
   }

double PDFHistogram::CalculateBarCountsFromBins(const Matrix & data_values)
   {
   //   double start = barStart != _NAN_ ? barStart : dataMin ;
   double max_count = 0.0;

   for (int k = 0; k < data_values.rows; k++) // each row past two is a new bar series
      {
      for (int i = 0; i < data_values.cols; i++)   // each column is a bar in a bar series.
         {
         bars[k][i].hasData = false;

         if (!xAxis.useDiscreteValues)
            {
            bars[k][i].lowerBound = xRanges[0][i];
            bars[k][i].upperBound = xRanges[1][i];
            }
         else    // if using categorical, lower and upper is set to the same...
            {
            bars[k][i].lowerBound = xRanges[0][i] - 0.5;
            bars[k][i].upperBound = xRanges[1][i] + 0.5;
            }
         bars[k][i].count = data_values[k][i];
         bars[k][i].isInitialized = true;
         bars[k][i].hasData = true;
         max_count = max(bars[k][i].count, max_count);
         }
      if (useSeriesPercentages ) ScaleToSeriesPct(k);
      }
   if (useTotalPercentages) ScaleToTotalPct();

   if (useSeriesPercentages || useTotalPercentages) max_count = 1.0;

   return max_count;
   }

void PDFHistogram::AssignBarTags()
   {
   for (int series = 0; series < dataValues.rows; series++)
      {
      if (taggedValues[series] == _NAN_) continue;
      if (taggedValues[series] < bars[series][0].lowerBound) continue;
      if (taggedValues[series] > bars[series][numBars - 1].upperBound) continue;

      int bar = 0;
      while (taggedValues[series] > bars[series][bar].upperBound)
         bar++;

      if (bars[series][bar].count == 0) continue;

      bars[series][bar].SetTag((const char *) valueTags[series]);
      }
   }

int PDFHistogram::CalculateBarRanges(int times_tried)
   {
   int new_num_bars;
   bool fixedMin = dataMin != _NAN_;
   bool fixedMax = dataMax != _NAN_;

   if (!fixedMin)
      dataMin = dataValues.SafeMin();

   if (!fixedMax)
      dataMax = dataValues.SafeMax();

   if (dataMax < dataMin)  { double tmp = dataMax; dataMax = dataMin; dataMin = tmp; }

   new_num_bars = CalculateBarWidth(times_tried);

   if (!fixedMin) dataMin = xAxis.FloorEPS(dataMin / barWidth) * barWidth;
   if (!fixedMax) dataMax = xAxis.CeilEPS(dataMax / barWidth) * barWidth;

   new_num_bars = CalculateBarWidth(times_tried);

   return new_num_bars;
   }

int PDFHistogram::CalculateBarWidth(int times_tried)
   {
   double start = barStart != _NAN_ ?  barStart : dataMin;
   double end = barEnd != _NAN_ ? barEnd : dataMax;

   double range = fabs(end - start);
   if (range == 0.0)
      range = 1.0;

   double num_bars;

   if (keepBarWidth)
      {
      num_bars = range / barWidth;
      tooManyBars = (num_bars > PDFHIST_MAX_BARS && PDFHIST_MAX_BARS != _NAN_);  // turn off feature for now
      }
   else
      {
      double mag = floor(log10(range));
      double factor_up[]   = {2.0, 2.0, 1.25, 2.0};
      double factor_down[] = {0.5, 0.8, 0.5, 0.5};
      int lower_bound = 5 + 5 * times_tried;
      int upper_bound = 12 + 5 * times_tried;
      int cycle = 0;

      barWidth = exp(mag * log(10.0));
      num_bars = range / barWidth;

      bool tooLarge = num_bars < lower_bound;
      bool tooSmall = num_bars > upper_bound;

      if (tooLarge || tooSmall)
         {
         while (true)
            {
            barWidth *= tooLarge ? factor_down[cycle % 4] : factor_up[cycle % 4];
            num_bars = range / barWidth;
            if (tooLarge &&  num_bars  >= lower_bound) break;
            if (tooSmall && num_bars <= upper_bound) break;
            cycle++;
            }
         }
      }

   return int(ceil(num_bars));
   }

void PDFHistogram::AllocateBarMatrix(int num_series, int num_bars_per_series)
   {
   if (!bars && num_series > 0 && num_bars_per_series > 0)
      {
      bars = new PDFChartBar *[num_series];
      for (int i = 0; i < num_series; i++)
         bars[i] = new PDFChartBar[num_bars_per_series];
      }
   else
      ResizeBarMatrix(num_series, num_bars_per_series);

   numBars = num_bars_per_series;
   barSeries = num_series;  
   }

void PDFHistogram::ResizeBarMatrix(int new_series, int new_bars_per_series)
   {
   PDFChartBar ** new_bars = NULL;

   if (new_series > 0)
      new_bars = new PDFChartBar *[new_series];

   for (int i = 0; i < new_series; i++)
      new_bars[i] = new_bars_per_series > 0 ? new PDFChartBar[new_bars_per_series] : NULL;

   for (int i = 0; i < new_series; i++)
      if (i < barSeries)
         for (int j = 0; j < new_bars_per_series; j++)
            if (j < numBars)
               new_bars[i][j] = bars[i][j];

   for (int i = 0; i < barSeries; i++)
      if (bars[i])
         delete [] bars[i];

   if (bars)
      delete [] bars;

   bars = new_bars;
   }

double PDFHistogram::CalculateBarCounts(const Matrix & sorted)
   {
   double start = barStart != _NAN_ ? barStart : dataMin ;
   double max_count = 0.0;

   for (int k = 0; k < sorted.rows; k++)
      {
      int    datapt = 0;
      double curr_value = sorted[k][0];

      for (int i = 0; i < numBars; i++)
         {
         bars[k][i].count = 0;
         bars[k][i].hasData = false;

         bars[k][i].lowerBound = start + i * barWidth;
         bars[k][i].upperBound = start + (i + 1) * barWidth;

         double tol = min( 1.0,  1e-6 * (dataMax - start));
         double bar_upper_tol = bars[k][i].upperBound + tol;

         while (curr_value < bar_upper_tol && datapt < sorted.cols)
            {
            if (curr_value != _NAN_)
               {
               bars[k][i].count += 1.0;
               bars[k][i].hasData = true;
               }
            datapt++;
            curr_value =  datapt < sorted.cols ? sorted[k][datapt] : _NAN_;
            }

         bars[k][i].isInitialized = true;
         max_count = max(bars[k][i].count, max_count);
         }
      if (useSeriesPercentages ) ScaleToSeriesPct(k);
      }

   if (useTotalPercentages) ScaleToTotalPct();
   if (useSeriesPercentages || useTotalPercentages) max_count = 1.0;

   return max_count;
   }

void PDFHistogram::ScaleToSeriesPct(int data_series)
   {
   double row_count = 0.0;
   for (int i = 0; i < numBars; i++)
      row_count += bars[data_series][i].count;

   row_count = row_count > 0.0 ? 1.0 / (double) row_count : 0.0;

   for (int i = 0; i < numBars; i++)
      bars[data_series][i].count *= row_count;
   }

void PDFHistogram::ScaleToTotalPct()
   {
   double overall_total = 0.0;

   for (int k = 0; k < barSeries; k++)
      for (int i = 0; i < numBars; i++)
          overall_total += bars[k][i].count;

   overall_total = overall_total > 0.0 ? 1.0 / overall_total : 0.0;

   for (int k = 0; k < barSeries; k++)
      for (int i = 0; i < numBars; i++)
         bars[k][i].count *= overall_total;
   }

/// Private drawing subroutines called directly or indirectly from base class Draw()////

void PDFHistogram::DrawBody(PDF & pdf)
   {
   pdf.page.SetLineWidth(0);
   pdf.page.SetLineColor(0, 0, 0);
   pdf.page.SetLineStyle(lsSolid);
   pdf.page.SetFillColor(0, 0, 0);

   if (!tooManyBars)
      stacked ?  DrawBarsStacked(pdf) : DrawBarsAligned(pdf);
   else
      DrawErrorGraph(pdf);
   }

void PDFHistogram::DrawErrorGraph(PDF & pdf)
   {
   pdf.page.hTextAlignment = taCenter;
   pdf.page.vTextAlignment = taMiddle;
   pdf.page.WriteText(0.5 * (x0 + x1), 0.5 * (y0 + y1), "Data doesn't fit");
   }

void PDFHistogram::DrawBarsStacked(PDF & pdf)
   {
   double * bar_y0 = new double[numBars];
   double   rect_x0, rect_y0, rect_x1, rect_y1;

   barSpacing = keepBarSpacing ? barSpacingPct * barWidth : 0.0;

   for (int k = 0; k < dataValues.rows; k++)
      {
      for (int i = 0; i < numBars; i++)
         {
         useColor ? lines[k].SelectColorBrush(pdf) : lines[k].SelectGrayBrush(pdf);

         if (k == 0)
            bar_y0[i] = 0.0;
         else if (bars[k - 1][i].count > 0.0)
            bar_y0[i] += bars[k-1][i].count;

         rect_x0 = bars[k][i].lowerBound + 100 * xAxis.Eps() + barSpacing * 0.5;
         rect_x1 = bars[k][i].upperBound - 100 * xAxis.Eps() - barSpacing * 0.5;
         rect_y0 = bar_y0[i];
         rect_y1 = bars[k][i].count + bar_y0[i];

         if (bars[k][i].hasData)
            {
            pdf.page.FillRectangle(MapX(rect_x0), MapY(rect_y0), MapX(rect_x1), MapY(rect_y1));
            pdf.page.DrawRectangle(MapX(rect_x0), MapY(rect_y0), MapX(rect_x1), MapY(rect_y1));
            if (bars[k][i].isTagged)
               {
               double last_font_size = pdf.page.GetFontSize();
               pdf.page.SetFillColor(0.0, 0.0, 0.0);
               pdf.page.SetFontSize(last_font_size * 2.0);
               pdf.page.hTextAlignment = taCenter;
               pdf.page.vTextAlignment = taMiddle;
               double label_height =  pdf.page.TextHeight((const char *) bars[k][i].barTag);
               pdf.page.WriteText(MapX((rect_x0 + rect_x1) * 0.5), MapY(rect_y1) + 0.5 * label_height,
                                  (const char *) bars[k][i].barTag);
               pdf.page.SetFontSize(last_font_size);
               }
            }
         }
      }

   if (bar_y0)
      delete [] bar_y0;
   }

void PDFHistogram::DrawBarsAligned(PDF & pdf)
   {
   double bar_x0, bar_x1;

   // set spacing between bars; use less spacing as the number of bar series increases
   barSpacing =  (keepBarSpacing ? barSpacingPct : 0.2) * barWidth * (3.0 / (barSeries + 3.0));

   double subBarWidth = (barWidth - barSpacing)/barSeries;
   double barEps = 100 * xAxis.Eps();

   for (int k = 0; k < dataValues.rows; k++)
      {
      bool row_has_data = false;
      for (int i = 0; i < numBars; i++)
         row_has_data |= (bars[k][i].hasData);

      for (int i = 0; i < numBars; i++)
         {
         useColor ? lines[k].SelectColorBrush(pdf) : lines[k].SelectGrayBrush(pdf);

         bar_x0 = bars[k][i].lowerBound + barSpacing * 0.5 + k * subBarWidth;
         bar_x1 = bar_x0 + subBarWidth;

         if (bars[k][i].hasData)
            {
            pdf.page.FillRectangle(MapX(bar_x0 + barEps), MapY(0.0),
               MapX(bar_x1 - barEps), MapY(bars[k][i].count));
            pdf.page.DrawRectangle(MapX(bar_x0 + barEps), MapY(0.0),
               MapX(bar_x1 - barEps), MapY(bars[k][i].count));
            if (bars[k][i].isTagged)
               {
               double last_font_size = pdf.page.GetFontSize();
               pdf.page.SetFillColor(0.0, 0.0, 0.0);
               pdf.page.SetFontSize(last_font_size * 2.0);
               pdf.page.WriteText(MapX(bar_x1 - barEps), MapY(bars[k][i].count), (const char *) bars[k][i].barTag);
               pdf.page.SetFontSize(last_font_size);
               }
            }
         else if (!row_has_data  && bars[k][i].lowerBound <= dataMax)
            DrawNoDataLabel(pdf, MapX( 0.5 * (bar_x0 + bar_x1)));
         }
      }
   }

void PDFHistogram::DrawNoDataLabel(PDF & pdf, double x_pos)
   {
   // Set font size to at most 0.5 titleFontSize, half-maximal when barWidth is
   // 1/20 of xAxis range

   pdf.page.SetFontSize(0.5 * titleFontSize * barWidth / (barWidth + 0.05 * (xAxis.max  - xAxis.min)));

   double label_width = pdf.page.TextWidth("NO DATA");

   pdf.page.SetFontOrientation(90.0);
   pdf.page.vTextAlignment = taMiddle;
   pdf.page.hTextAlignment = taCenter;

   pdf.page.WriteText(x_pos ,MapY(0.0) + label_width * 1.2, "NO DATA");
   pdf.page.SetFontOrientation(0.0);
   }

void PDFHistogram::DrawLegend(PDF & pdf)
   {
   legend.Draw(pdf, Space(), useColor);
   }

///// Chart settings/preferences////////////////////////////////////////////////

void PDFHistogram::SetBarWidth(double bar_width)
   {
   barWidth = bar_width;
   if (!yAxis.keepMax)  yAxis.max = _NAN_;
   if (!yAxis.keepStep) yAxis.step = _NAN_;
   keepBarWidth = (barWidth != _NAN_);
   }

void PDFHistogram::SetBarStart(double bar_start)
   {
   barStart = bar_start;
   if (!yAxis.keepMax) yAxis.max = _NAN_;
   if (!yAxis.keepStep) yAxis.step = _NAN_;
   }

void PDFHistogram::SetBarEnd(double bar_end)
   {
   barEnd = bar_end;
   if (!yAxis.keepMax) yAxis.max = _NAN_;
   if (!yAxis.keepStep) yAxis.step = _NAN_;
   }

void PDFHistogram::SetBarSpacing(double bar_spacing)
   {
   barSpacingPct = bar_spacing;
   keepBarSpacing = barSpacingPct != _NAN_;
   }

void PDFHistogram::SetUseSeriesPercentages(bool use_series_pct)
   {
   useSeriesPercentages = use_series_pct;
   if (useSeriesPercentages)
      {
      useTotalPercentages = false;
      if (!yAxis.keepMax) yAxis.max = _NAN_;
      if (!yAxis.keepStep) yAxis.step = _NAN_;
      }
   }

void PDFHistogram::SetUseTotalPercentages(bool use_total_pct)
   {
   useTotalPercentages = use_total_pct;
   if (useTotalPercentages)
      {
      useSeriesPercentages = false;
      if (!yAxis.keepMax) yAxis.max = _NAN_;
      if (!yAxis.keepStep) yAxis.step = _NAN_;
      }
   }

void PDFHistogram::SetDataValueTag(int series, double value, const char * tag)
   {
   if (series < 0 || series > dataValues.rows)
      error("PDFhistogram::Can't tag value %lf in data series %d - "
            "data series doesn't exist", value, series);

   taggedValues[series] = value;
   valueTags[series] = tag;
   }

void PDFHistogram::SetCategoricalBars(bool categorical, double first_bar_category, double last_bar_category)
   {
   xAxis.useDiscreteValues = categorical;

   if (categorical)
      {
      SetBarWidth(1.0);
      SetBarSpacing(0.4);
      SetBarEnd(last_bar_category);

      xAxis.minMax = last_bar_category != _NAN_ ? last_bar_category + 0.5 : _NAN_;
      xAxis.maxMin = first_bar_category != _NAN_ ? first_bar_category - 0.5 : _NAN_;

      SetBarStart(xAxis.maxMin);
      xAxis.SetStepStart(xAxis.maxMin);
      }
   else
      {
      SetBarStart(_NAN_);
      SetBarWidth(_NAN_);
      SetBarSpacing(_NAN_);
      SetBarEnd(_NAN_);
      xAxis.SetStepStart(_NAN_);
      xAxis.minMin = _NAN_;
      xAxis.maxMin = _NAN_;
      }
   }










 
 
 
 
 
