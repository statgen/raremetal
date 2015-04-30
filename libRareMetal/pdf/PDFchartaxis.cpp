////////////////////////////////////////////////////////////////////// 
// janpdf/PDFchartaxis.cpp 
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
#include "PDF.h"
#include "PDFchartaxis.h"
#include "MathConstant.h"
#include "Constant.h"
#include <math.h>


PDFChartAxis::PDFChartAxis(bool v) :
   vert(v)
   {
   Init();

   digits = 1;
   tickSize = _NAN_;
   tickLabelSize = 0.0;
   width =  height = 0.0;
   includesZero = false;
   }

PDFChartAxis::~PDFChartAxis()
   {
   }

void PDFChartAxis::Init()
   {
   min = max = _NAN_;
   minMin = maxMin = minMax = maxMax = _NAN_;
   step = stepStart = stepMin = _NAN_;
   keepStep = keepMin = keepMax = false;
   centerData = false;
   frameData = false;
   useDiscreteValues = false;

   font = fHelvetica;
   tickLabelFontSize = 0.0;
   labelFontSize = 0.0;
   keepFontSize = false;
   useScientific = false;
   minDigits = maxDigits = -1;
   usePow10Labels = false;
 
   useTicks = useNumericLabels = true;
   useStringLabels = useAlternatingLabels = false;
   labelsOffset = false;
   }

void PDFChartAxis::SetTickLabelFontSize(double font_size)
   {
   tickLabelFontSize = font_size;
   keepFontSize = true;
   }

void PDFChartAxis::SetNumericTickLabels(bool use_numeric)
   {
   useNumericLabels = use_numeric;
   if (use_numeric && useStringLabels)
      useStringLabels = false;
   }

void PDFChartAxis::SetStringTickLabels(bool use_string)
   {
   useStringLabels = use_string;
   if  (use_string && useNumericLabels)
      useNumericLabels = false;
   }

void PDFChartAxis::Reset()
   {
   Init();
   label = "";
   stringTickLabels.Clear();
   }

void PDFChartAxis::GetDataRange(const Matrix & values)
   {
   bool fixedMin = min != _NAN_;
   bool fixedMax = max != _NAN_;

   if (!fixedMin)
      min = values.SafeMin();

   if (!fixedMax)
      max = values.SafeMax();

   if (max < min)  { double tmp = max; max = min; min = tmp; }
   if (max == min)
      {
      max += 0.5;
      min -= 0.5;
      }

   if (centerData) Center(values);
   if (frameData) Frame(values);

   if (minMin != _NAN_ && min < minMin)  min = minMin;
   if (maxMin != _NAN_ && min > maxMin)  min = maxMin;
   if (minMax != _NAN_ && max < minMax) max = minMax;
   if (maxMax != _NAN_ && max > maxMax) max = maxMax;

   if (!keepStep) CalculateStep();

   if (!fixedMin) min = FloorEPS(min /step) * step;
   if (!fixedMax) max = CeilEPS(max /step) * step;

   CalculateDigits();
   if (min < 0.0 && max > 0.0) includesZero = true;
   }

void PDFChartAxis::Center(const Matrix & values)
   {
   double mean = values.SafeMean();
   double diff_min = mean - min;
   double diff_max = max - mean;

   double delta = ::max(diff_min, diff_max);
   double min_delta = ::min((max - min) * 0.1, 1.0);
   delta = ::max(delta, min_delta);
   min = mean - delta;
   max = mean + delta;
   }

void PDFChartAxis::Frame(const Matrix & values)
   {
   double data_min = values.SafeMin();
   double data_max = values.SafeMax();

   if (max == data_max || min == data_min)
      {
      double diff_max = fabs(max - data_max);
      double diff_min = fabs(data_min - min);

      double delta = ::max(diff_min, diff_max);
      double min_delta = ::min((max - min)* 0.1, 0.1);

      delta = ::max(delta, min_delta);

      if (min == data_min) min -= delta;
      if (max == data_max) max += delta;
      }
   }

void PDFChartAxis::SetMin(double new_min)
   {
   minMin = new_min;
   maxMin = new_min;
   min = new_min;
   keepMin = (min != _NAN_);
   }

void PDFChartAxis::SetMax(double new_max)
   {
   minMax = new_max;
   maxMax = new_max;
   max = new_max;
   keepMax = (max != _NAN_);
   }

void PDFChartAxis::SetStep(double new_step)
   {
   step = new_step;
   keepStep = (step != _NAN_);
   }

void PDFChartAxis::SetStepStart(double step_start)
   {
   stepStart = step_start;
   }

void PDFChartAxis::SetStepMin(double step_min)
   {
   stepMin = step_min;
   }

void PDFChartAxis::CalculateStep()
   {
   double factor_up[]   = {2.0, 2.0, 1.25, 2.0};
   double factor_down[] = {0.5, 0.8, 0.5, 0.5};

   double start = (stepStart != _NAN_ ? stepStart: min);

   if (keepStep) return;

   double range = fabs(max - start);
   double mag = floor(log10(range));
   step = exp(mag * log(10.0));
   double steps = range / step;

   bool tooLarge = steps <= 3;
   bool tooSmall = steps > 6;

   int cycle = 0;
   if (tooLarge || tooSmall)
      {
      while (true)
         {
         step *= tooLarge ? factor_down[cycle % 4] : factor_up[cycle % 4];
         steps = range / step;
         if (tooLarge && steps >= 3) break;
         if (tooSmall && steps <= 6) break;
         cycle++;
         }
      }
   if ( stepMin != _NAN_ && step < stepMin )
      step = stepMin;
   if (useDiscreteValues)
      step = ::max(1.0, FloorEPS(step));

   }

void PDFChartAxis::CalculateDigits()
   {
   int mag_step = int(usePow10Labels ? step : log10(step));
   double axis_min = usePow10Labels ? pow(10, min) : min;
   double axis_max = usePow10Labels ? pow(10, max) : max;

   digits = 1;

   if (mag_step > 0 && (fabs(axis_min) > 1e5 || fabs(axis_max) > 1e5))
      {
      double abs_max = fabs(axis_max);
      int mag_max = (int) (abs_max != 0.0 ? FloorEPS(log10(abs_max)) : 0);

      double abs_min = fabs(axis_min);
      int mag_min = (int) (abs_min != 0.0 ? FloorEPS(log10(abs_min)): 0);
      int magnitude = ::max(mag_max, mag_min);

      useScientific = true;
      digits = (magnitude > mag_step ? 0 : 1);
      }
   else
      {
      if (mag_step < -5)
         {
         digits = 1;
         useScientific = true;
         }
      else if (mag_step < 0)
         digits = int(-mag_step) + 2;
      else if (mag_step > 1)
         digits = 0;
      }

   if (maxDigits != -1) 
     digits = ::min(digits, maxDigits);

   if (minDigits != -1)
     digits = ::max(digits, minDigits); 
   }

void PDFChartAxis::SetMinDigits(int digits)
   {
   if (digits < 0) 
      error("PDFchartaxis::Can't set axis label digits to a negative number");

   minDigits = digits;
   }

void PDFChartAxis::SetMaxDigits(int digits)
   {
   if (digits < 0)
      error("PDFchartaxis::Can't set axis label digits to a negative number");

   maxDigits = digits;
   }

void PDFChartAxis::SetTickLabel(int i, const char * tick_label)
  {
  if (i >= stringTickLabels.Length())
     stringTickLabels.Dimension(i+1);

  stringTickLabels[i] = tick_label;
  }

void PDFChartAxis::SetTickLabelHighlight(int i, bool highlight)
  {
  if (i >= tickLabelHighlights.Length())
     tickLabelHighlights.Dimension(i+1);

  tickLabelHighlights[i] = highlight ? 1 : 0;
  }










 
 
 
 
 
 
 
 
