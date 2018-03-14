#include "MetaUtility.h"
#include <math.h>
#include <complex>
//#include <gsl/gsl_integration.h> // integration for Dajiang's method

bool SetIfilePosition(IFILE &sfile, Tabix &myTabix, String Chr, int pos)
{
    const char *chr = Chr.c_str();
    uint64_t fstart = 0;
    bool st = myTabix.getStartPos(chr, pos, fstart);
    if (!st)
    {
        printf("Warning:[SetIfilePosition] Unable to locate %s:%d\n", chr, pos);
        return false;
    }
// seek position now with fstart
//      if ( fstart != (uint64_t)iftell(sfile) ) {
    if (ifseek(sfile, fstart, SEEK_SET) != true)
    {
        printf("Unable to seek position %s:%d\n", chr, pos);
        exit(1);
    }
//      }


    if (ifeof(sfile))
    {
        printf("Warning:[SetIfilePosition] Reached end!\n");
        return false;
    }
// then move on to that position
    String buffer;
    StringArray tmp;
    bool found = false;
    while (!ifeof(sfile))
    {
        buffer.ReadLine(sfile);
        tmp.Clear();
        tmp.AddTokens(buffer, "\t");
        if (tmp[1].AsInteger() >= pos)
        {
            found = true;
            break;
        }
    }
    return found;
}

// RMW: adjust = 0;
// RVtest: adjust = 1;
void tellRvOrRmw(String &buffer, bool &adjust, int marker_col, int cov_col)
{
    if (buffer.Find("RareMetalWorker") == -1)
    { // rvt
        adjust = 1;
    } else
    {  //rmw
        adjust = 0;
    }
    setFromRvOrRmwAdjust(adjust, marker_col, cov_col);
}

void setFromRvOrRmwAdjust(bool adjust, int marker_col, int cov_col)
{
    if (adjust)
    { // rv test
        marker_col = 4;
        cov_col = 5;
    } else
    { // rmw
        marker_col = 2;
        cov_col = 3;
    }
}


// open result file as prefix.meta.method.results
void openMetaResultFile(String &prefix, String &filename, IFILE &output, String &method)
{
    if (prefix == "")
    {
        filename = "meta." + method + ".results";
    } else if (prefix.Last() == '.' || prefix.Last() == '/')
    {
        filename = prefix + "meta." + method + ".results";
    } else
    {
        filename = prefix + ".meta." + method + ".results";
    }
    output = ifopen(filename, "w", InputFile::UNCOMPRESSED);
}


double GetGenomicControlFromPvalue(Vector &pvalue)
{
    if (pvalue.Length() < 1)
    {
        error("No p value available for calculating Genomic Control. Check if markers were filtered!");
    }
    Vector tmp;
    tmp.Copy(pvalue);
    tmp.Sort();
    double pvalue_median = tmp[0.5];
    double GC = qchisq(pvalue_median, 1, 0, 0);
    GC /= 0.456;
    return GC;
}


double GetBetaDensity(double a, double b, double x)
{
    double density;
    //density = exp(gammln(a+b)-gammln(a)-gammln(b)+(a-1.0)*log(x)+(b-1.0)*log(1.0-x));
    density = dbeta(x, a, b, 0);
    return density;
}


double CalculateCorrCoef(Vector &a, Vector &b)
{
    double coef;
    int n = a.Length();
    double sum_a = a.Sum();
    double sum_b = b.Sum();
    coef = n * a.InnerProduct(b) - sum_a * sum_b;
    coef /= sqrt((n * a.SumSquares() - sum_a * sum_a) * (n * b.SumSquares() - sum_b * sum_b));
    return coef;
}

void RevertAllele(String SNP, String &newSNP)
{
    StringArray tmp;
    tmp.AddTokens(SNP, ":");
    newSNP = tmp[0] + ":" + tmp[1] + ":" + tmp[3] + ":" + tmp[2];
}

void removeChrFromString(String &chr_token)
{
    if (chr_token.Find("chr") != -1)
    {
        chr_token = chr_token.SubStr(3);
    }
}

/*
a = log(y.binary.mean / (1-y.binary.mean))
f = exp(a+e) / (1+exp(a+e)^2 * dnorm(e))
m = integrate: f(-500,500)
*/
/*
// density function for normal distribution N(0,sigma)
double dnorm_mean0( double x )
{
    double p;
    double cst = 0.3989423;
    p = cst * exp(-x*x/2);
    return p;
}

// integration of c
double c_inte( double x, void * params )
{
    double a = *(double*) params;
    double f = exp(a+x) / ( 1+exp(a+x)*exp(a+x) * dnorm_mean0(x) );
    return f;
}

double getCorrectC( double ymean,double sigma)
{
    double alpha = std::log( ymean / (1-ymean));

    gsl_integration_workspace* w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
//    F.function = &f;
    F.params = &alpha;
    double result, error;
    gsl_integration_qags( &F,0,1,0,1e-7,1000,w,&result,&error );
    double c = result;
    gsl_integration_workspace_free(w);
    
    return c;
}
*/

