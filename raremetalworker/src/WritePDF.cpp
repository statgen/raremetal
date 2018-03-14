#include <math.h>
#include "WritePDF.h"

#define MATHLIB_STANDALONE

#include <Rmath.h>
#include "QuickIndex.h"
#include "PDFchartline.h"

WritePDF::WritePDF()
{}

int WritePDF::resolution = 1000;
bool WritePDF::thinPoints = false;

void
WritePDF::Draw(PDF &pdf, StringArray &geneLabels, Vector &pvalueAll, Vector &pvalue1, Vector &pvalue5, StringArray &chr,
               Vector &pos, String title, String demo, bool useDataLabels)
{
    // Setup Q-Q plot
    if (pvalueAll.Length() > 1000000000 || !thinPoints)
    {
        thinPoints = true;
    }
    String plot_title = "";
    if (pvalueAll.Length() > 0)
    {
        DrawQQplot(pdf, pvalueAll, title, plot_title, demo);
    }
    plot_title = "(maf<0.01)";
    String dummy = "";
    if (pvalue1.Length() > 0)
    {
        DrawQQplot(pdf, pvalue1, title, plot_title, dummy);
    }
    plot_title = "(maf<0.05)";
    if (pvalue5.Length() > 0)
    {
        DrawQQplot(pdf, pvalue5, title, plot_title, dummy);
    }
    if (pvalueAll.Length() > 0)
    {
        DrawManhattanPlot(pdf, geneLabels, pvalueAll, chr, pos, title, useDataLabels);
    }
}

void WritePDF::Draw(PDF &pdf, StringArray &geneLabels, Vector &pvalueAll, StringArray &chr, Vector &pos, String title,
                    String extraTitle, String demo, bool useDataLabels)
{
    // Setup Q-Q plot
    thinPoints = true;
    if (pvalueAll.Length() > 0)
    {
        DrawQQplot(pdf, pvalueAll, title, extraTitle, demo);
        DrawManhattanPlot(pdf, geneLabels, pvalueAll, chr, pos, title, useDataLabels);
    }
}

void WritePDF::DrawOverlayQQplot(PDF &pdf, Vector &pvalueAll_, Vector &pvalue1_, Vector &pvalue5_, String title)
{
    pdf.page.SetSize(800, 800);
    //Setup pvalue
    Vector pvalueAll, pvalue1, pvalue5;
    pvalueAll.Copy(pvalueAll_);
    pvalue1.Copy(pvalue1_);
    pvalue5.Copy(pvalue5_);
    int N = pvalueAll.Length();
    int N1 = pvalue1.Length();
    int N5 = pvalue5.Length();
    int col = N + N1 + N5;
    Matrix y;
    y.Dimension(5, col, _NAN_);

    pvalueAll.Sort();
    pvalue1.Sort();
    pvalue5.Sort();
    int step = N / 100;

    for (int i = 0; i < N; i++)
    {
        y[0][i] = -log10((i * 1.0 + 0.5) / N);
        if (i % step == 0)
        {
            y[1][i] = y[0][i];
        }
        y[2][i] = -log10(pvalueAll[i]);
    }
    for (int i = N; i < N + N1; i++)
    {
        y[0][i] = -log10((i * 1.0 - N + 0.5) / N1);
        y[3][i] = -log10(pvalue1[i - N]);
    }
    for (int i = N + N1; i < N + N1 + N5; i++)
    {
        y[0][i] = -log10((i * 1.0 - N - N1 + 0.5) / N5);
        y[4][i] = -log10(pvalue5[i - N - N1]);
    }
    // Setup Q-Q plot
    PDFLineChart chart;
    chart.Reset();

    chart.title.printf(title);
    chart.xAxis.label.printf("expected -log10(pvalue)");
    chart.yAxis.label.printf("observed -log10(pvalue)");

    chart.xAxis.SetMin(0.0);
    chart.xAxis.SetMax(-log10(0.5 / N));
    chart.yAxis.SetMin(0.0);
    chart.yAxis.SetMax(-log10(pvalueAll.Min()));

    chart.Dimension(4, N);
    chart.ShowMarker(0, false);
    chart.ShowMarker(1, true);
    chart.ShowMarker(2, true);
    chart.ShowMarker(3, true);
    chart.ShowLine(0, true);
    chart.ShowLine(1, false);
    chart.ShowLine(2, false);
    chart.ShowLine(3, false);
    //PDFLineStyle style[3]={ lsSolid, lsDashed, lsDotted };
    //chart.SetLineStyle(0,style[0]);
    // Draw Q-Q plot
    chart.SetSeriesLabel(0, "reference");
    chart.SetSeriesLabel(1, "all vars");
    chart.SetSeriesLabel(2, "maf<0.01");
    chart.SetSeriesLabel(3, "maf<0.05");

    chart.SetLineWeight(0, 0.01);
    chart.SetMarkerShape(1, lmCircle);
    chart.SetMarkerShape(2, lmSquare);
    chart.SetMarkerShape(3, lmCross);
    chart.SetMarkerColor(1, 1.0, 0.0, 0.0);
    chart.SetMarkerColor(2, 0.0, 1.0, 0.0);
    chart.SetMarkerColor(3, 0.0, 0.0, 0.0);
    chart.SetMarkerRadius(1, 1.2);
    chart.SetMarkerRadius(2, 1.5);
    chart.SetMarkerRadius(3, 1.2);
    pdf.page.SetFontSize(14);

    chart.SetDataValues(y);
    chart.Draw(pdf);
}

void WritePDF::DrawQQplot(PDF &pdf, Vector &pvalueAll_, String title, String extraTitle, String demo)
{
    pdf.page.SetSize(800, 800);
    QuickIndex idx(pvalueAll_);

    PDFLineChartwithPolygon chart;

    chart.demo = demo;
    int N = pvalueAll_.Length();
    Matrix y;
    y.Dimension(3, N, _NAN_);

    for (int i = 0; i < N; i++)
    {
        //y[0][i] = floor((-log10((i*1.0+0.5)/N))*100+0.5)/100;
        //y[2][i] =  floor((-log10(pvalueAll_[idx[i]]))*100+0.5)/100;
        double p = -log10(pvalueAll_[idx[i]]);
        if (thinPoints)
        {
            if (p <= 2)
            {
                y[2][i] = fround(p, 2);
                y[0][i] = fround(-log10((i * 1.0 + 0.5) / N), 2);
            } else
            {
                y[2][i] = p;
                y[0][i] = -log10((i * 1.0 + 0.5) / N);
            }
        } else
        {
            y[2][i] = p;
            y[0][i] = -log10((i * 1.0 + 0.5) / N);
        }
    }

    //when ties plot the first point
    int bm = 0;
    for (int m = 0; m < N; m++)
    {
        if (y[2][m] <= 3)
        {
            bm = m;
            break;
        }
    }
    int tieCounter = 0;
    double tie = y[2][0];
    for (int k = bm; k < N; k++)
        //for(int k=0;k<N;k++)
    {
        if (y[2][k] == tie)
        {
            tieCounter++;
            if (tieCounter % 10 != 0)
            {
                y[2][k] = _NAN_;
            }
        } else
        {
            tie = y[2][k];
            tieCounter = 0;
        }
    }

    //setup the expected line
    y[1][0] = y[0][0];
    y[1][N - 1] = y[0][N - 1];

    chart.Reset();

    chart.title.printf(title + " " + extraTitle);
    chart.xAxis.label.printf("expected -log10(pvalue)");
    chart.yAxis.label.printf("observed -log10(pvalue)");

    chart.xAxis.SetMin(0.0);
    chart.xAxis.SetMax(-log10(0.5 / N) + 0.01 * (-log10(0.5 / N)));
    chart.yAxis.SetMin(0.0);
    chart.yAxis.SetMax(-log10(pvalueAll_.Min()) + 0.01 * (-log10(pvalueAll_.Min())));

    chart.Dimension(2, N);
    chart.ShowMarker(0, false);
    chart.ShowMarker(1, true);
    chart.ShowLine(0, true);
    chart.ShowLine(1, false);
    chart.useLegend = false;

    chart.SetLineWeight(0, 0.01);
    chart.SetMarkerShape(1, lmCircle);
    chart.SetMarkerColor(1, 1.0, 0.0, 0.0);
    chart.SetMarkerRadius(1, 2.0);

    chart.SetDataValues(y);

    //fill between the CI
    int n = N < 1000 ? N : 1000;
    double *poly_x = new double[2 * n];
    double *poly_y = new double[2 * n];
    //set values to poly_x and poly_y
    int step = N / n;
    double conf = 0.05;
    for (int i = 0; i < N; i++)
    {
        int mod = i % step;
        if (mod == 0)
        {
            int j = i / step;
            poly_x[j] = poly_x[2 * n - j - 1] = -log10((i * 1.0 + 0.5) / N);
            poly_y[j] = -log10(qbeta((1.0 - conf / 2), i + 1, N - i, 1, 0));
            poly_y[2 * n - 1 - j] = -log10(qbeta(conf / 2, i + 1, N - i, 1, 0));
        }
    }

    chart.poly_x = poly_x;
    chart.poly_y = poly_y;
    chart.poly_n = 2 * n;
    chart.Draw(pdf);

    if (poly_x)
    { delete[] poly_x; }
    if (poly_y)
    { delete[] poly_y; }
}

void WritePDF::DrawManhattanPlot(PDF &pdf, StringArray &geneLabels, Vector &pvalue, StringArray &chr, Vector &pos,
                                 String title, bool useDataLabels)
{
    pdf.page.SetSize(1200, 400);
    //Sanity check
    int N = pvalue.Length();
    if (pos.Length() != N || chr.Length() != N)
    {
        error("Vector of pvalue does not have the same length as chr and pos.\n");
    }

    PDFmanhattan chart;
    chart.Reset();
    if (useDataLabels)
    {
        useDataLabels = chart.useDataLabels;
    }
    //Setup data value for chart
    StringArray chr_cp;
    chr_cp = chr;
    chr_cp.Sort();
    StringArray tickLabels;
    tickLabels.Push(chr_cp[0]);
    String last_label = chr_cp[0];

    //this is for plots only

    for (int i = 1; i < chr.Length(); i++)
    {
        if (chr[i].AsDouble() > 22 || chr[i] == "Y" || chr[i] == "MT")
        {
            continue;
        }
        if (chr_cp[i] != last_label)
        {
            last_label = chr_cp[i];
            tickLabels.Push(chr_cp[i]);
        }
    }
    tickLabels.Sort();
    StringArray tmpLabels, tmp_del;
    tmpLabels.Dimension(tickLabels.Length());
    int tmp_i = 0;
    for (int i = 0; i < tickLabels.Length(); i++)
    {
        if (tickLabels[i].AsInteger() <= 22 && tickLabels[i].AsInteger() >= 1)
        {
            tmpLabels[tmp_i] = tickLabels[i];
            tmp_i++;
        } else
        {
            tmp_del.Push(tickLabels[i]);
        }
    }
    for (int i = 0; i < tmp_del.Length(); i++)
    {
        tmpLabels[tmp_i + i] = tmp_del[i];
    }
    for (int i = 0; i < tickLabels.Length(); i++)
    {
        tickLabels[i] = tmpLabels[i];
    }

    //first row is x, second row is threshold, then chr1, chr2...; last row has hits
    Matrix y;
    y.Dimension(tickLabels.Length() + 3, N, _NAN_);

    QuickIndex index(chr);

    int y_idx = 0;
    double previous_chr_end = 0.0;
    int previous_idx = 0;
    double threshold = -log10(0.05 / N);
    //double threshold = -log10(5e-08);
    for (int l = 0; l < tickLabels.Length(); l++)
    {
        chart.xAxis.stringTickLabels.Push(tickLabels[l]);
    }
    double total_range = 0.0;
    for (int i = 0; i < tickLabels.Length(); i++)
    {
        Vector pos_i, pvalue_i;
        StringArray geneLabels_i;
        for (int j = previous_idx; j < chr.Length(); j++)
        {
            int out = index[j];
            if (chr[out] != tickLabels[i])
            {
                continue;
            }
            double p = -log10(pvalue[out]);
            if (useDataLabels)
            {
                geneLabels_i.Push(geneLabels[out]);
            }
            if (thinPoints && p <= 2)
            {
                pvalue_i.Push(fround(p, 1));
            } else
            {
                pvalue_i.Push(p);
            }
            pos_i.Push(fround(pos[out] / resolution, 0));
        }
        //      previous_idx += pos_i.dim;

        //thin points pos_i pvalue_i pair
        QuickIndex pos_i_idx(pos_i);
        double pos_rep_start = pos_i[pos_i_idx[0]];
        Vector pvalue_same_pos;
        pvalue_same_pos.Push(pvalue_i[pos_i_idx[0]]);
        IntArray keep;
        keep.Push(pos_i_idx[0]);
        for (int k = 0; k < pos_i.Length(); k++)
        {
            if (pos_i[pos_i_idx[k]] != pos_rep_start)
            {
                pos_rep_start = pos_i[pos_i_idx[k]];
                pvalue_same_pos.Clear();
                pvalue_same_pos.Push(pvalue_i[pos_i_idx[k]]);
                keep.Push(pos_i_idx[k]);
                continue;
            }
            if (pvalue_same_pos.FastFind(pvalue_i[pos_i_idx[k]]) == -1)
            {
                keep.Push(pos_i_idx[k]);
            }
        }

        Vector pos_chr_i, pvalue_chr_i;
        StringArray geneLabels_chr_i;
        pos_chr_i.Dimension(keep.Length());
        pvalue_chr_i.Dimension(keep.Length());
        if (useDataLabels)
        {
            geneLabels_chr_i.Dimension(keep.Length());
        }
        for (int idx = 0; idx < keep.Length(); idx++)
        {
            pos_chr_i[idx] = pos_i[keep[idx]];
            pvalue_chr_i[idx] = pvalue_i[keep[idx]];
            if (useDataLabels)
            {
                geneLabels_chr_i[idx] = geneLabels_i[keep[idx]];
            }
        }

        QuickIndex pos_index(pos_chr_i);
        int last = pos_i_idx[pos_i.Length() - 1];
        int first = pos_i_idx[0];
        int max_pos = pos_i[last];
        int min_pos = pos_i[first];
        int range = max_pos - min_pos;
        total_range += range;
        //double chr_length = 100.0*range;
        //fill in the data matrix
        for (int k = 0; k < pos_chr_i.Length(); k++)
        {
            int pos_idx = pos_index[k];
            double rel_loc = 0.0;
            if (range == 0)
            {
                rel_loc = range;
            } else
            {
                rel_loc = pos_chr_i[pos_idx] - pos_chr_i[pos_idx - 1];
            }
            y[0][y_idx] = y[0][y_idx - 1] + rel_loc;
            y[i + 2][y_idx] = pvalue_chr_i[pos_idx];
            if (pvalue_chr_i[pos_idx] > threshold)
            {
                y[y.rows - 1][y_idx] = pvalue_chr_i[pos_idx];
                if (useDataLabels)
                {
                    chart.dataLabelsX.Push(y[0][y_idx]);
                    chart.dataLabelsY.Push(pvalue_chr_i[pos_idx]);
                    chart.dataLabels.Push(geneLabels_chr_i[pos_idx]);
                }
            }
            y_idx++;
        }
        previous_chr_end = y[0][y_idx - 1];
        //printf("pos for chr %s is %g\n",tickLabels[i].c_str(),previous_chr_end);
        chart.tickLabelsPos.Push(previous_chr_end - range / 2.0);
    }
    //delete _NAN_ columns in y
    for (int i = y.cols - 1; i >= 0; i--)
    {
        if (y[0][i] == _NAN_)
        {
            y.DeleteColumn(i);
        }
    }

    //set up the threshold data value
    N = y[0].dim;
    int step = 1;
    if (N > 10)
    {
        step = N / 10;
    }
    for (int i = 0; i < N; i++)
    {
        if (i % step == 0)
        {
            y[1][i] = threshold;
        }
    }
    y[1][N - 1] = threshold;
    int series = y.rows - 1;

    chart.Dimension(tickLabels.Length() + 2, y.cols);

    chart.title.printf(title);
    chart.xAxis.label.printf("chromosome");
    chart.yAxis.label.printf("-log10(pvalue)");

    chart.xAxis.SetMin(y[0][0] - (y[0][y.cols - 1] - y[0][0]) * 0.005);
    chart.xAxis.SetMax(y[0][y.cols - 1] + (y[0][y.cols - 1] - y[0][0]) * 0.005);
    chart.yAxis.SetMin(0.0);
    chart.yAxis.SetMax(-1.06 * log10(pvalue.Min()));

    //   chart.SetSeriesLabel(0,"expected");
    //   chart.SetSeriesLabel(1,legend.c_str());
    chart.useLegend = false;

    chart.xAxis.useTicks = true;
    chart.xAxis.SetNumericTickLabels(false);
    chart.xAxis.SetStringTickLabels(true);

    chart.SetDataValues(y);

    for (int i = 1; i < series - 1; i++)
    {
        chart.SetMarkerShape(i, lmCircle);
        /*
       if(i%6==0)
       chart.SetSeriesColor(i,0.0,0.0,0.0);
       else if(i%6==1)
       chart.SetSeriesColor(i,0.0,0.0,1.0);
       else if(i%6==2)
       chart.SetSeriesColor(i,0.0,1.0,0.0);
       else if(i%6==3)
       chart.SetSeriesColor(i,1.0,0.0,0.0);
       else if(i%6==4)
       chart.SetSeriesColor(i,0.0,1.0,1.0);
       else if(i%6==5)
       chart.SetSeriesColor(i,1.0,0.0,1.0);
         */
        if (i % 2 == 0)
        {
            chart.SetSeriesColor(i, 0.0, 0.0, 0.4);
        } else
        {
            chart.SetSeriesColor(i, 0.63, 0.67, 0.87);
        }
        chart.ShowMarker(i, true);
        chart.SetMarkerRadius(i, 2.0);
        chart.ShowLine(i, false);
    }
    chart.SetMarkerShape(series - 1, lmCircle);
    chart.SetMarkerRadius(series - 1, 3);
    chart.SetSeriesColor(series - 1, 0.0, 0.6, 0.0);
    chart.SetSeriesColor(0, 0.66, 0.66, 0.66);
    chart.ShowMarker(0, false);
    //show reference line
    chart.ShowLine(0, false);
    chart.ShowLine(y.rows - 1, false);
    chart.SetLineWeight(0, 0.001);
    chart.SetLineStyle(0, lsDashed);
    chart.xAxis.useTicks = false;

    chart.Draw(pdf);
}
