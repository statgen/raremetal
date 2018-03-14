// Title: Imhof (1961) algorithm
// Ref. (book or article): {J. P. Imhof, Computing the Distribution of Quadratic Forms in Normal Variables, Biometrika, Volume 48, Issue 3/4 (Dec., 1961), 419-426

// Description:
// Distribution function (survival function in fact) of quadratic forms in normal variables using Imhof's method.

//#include <R.h>
#include "Rmath.h"

extern "C" {

double theta(double *u, double *lambda, int *lambdalen, double *h, double *x, double *delta2)
{ // See Imhof (1961), p.423
    int i, m;
    double sum = 0.0;
    m = lambdalen[0];
    for (i = 1; i <= m; i = i + 1)
    {
        sum = sum + h[i - 1] * atan(lambda[i - 1] * u[0]) +
              delta2[i - 1] * lambda[i - 1] * u[0] /
              (1.0 + R_pow(lambda[i - 1] * u[0], 2.0));
    }
    sum = 0.5 * sum - 0.5 * x[0] * u[0];
    return (sum);
}

double rho(double *u, double *lambda, int *lambdalen, double *h, double *delta2)
{ // See Imhof (1961), p.423
    int i, m;
    double prod = 1.0;
    m = lambdalen[0];
    for (i = 1; i <= m; i = i + 1)
    {
        prod = prod * R_pow(1.0 + R_pow(lambda[i - 1] * u[0], 2.0), 0.25 * h[i - 1]) *
               exp(0.5 * delta2[i - 1] * R_pow(lambda[i - 1] * u[0], 2.0) /
                   (1.0 + R_pow(lambda[i - 1] * u[0], 2.0)));
    }
    return (prod);
}


double imhoffunc(double *u, double *lambda, int *lambdalen, double *h, double *x, double *delta2)
{ // This is the function under the integral sign in equation (3.2), Imhof (1961), p.422
    double theta(double *u, double *lambda, int *lambdalen, double *h, double *x, double *delta2);
    double rho(double *u, double *lambda, int *lambdalen, double *h, double *delta2);
    double res;
    res = (sin(theta(u, lambda, lambdalen, h, x, delta2))) / (u[0] * rho(u, lambda, lambdalen, h, delta2));
    return (res);
}


typedef void integr_fn(double *x, int n, void *ex);

void f(double *x, int n, void *ex)
{
    double imhoffunc(double *u, double *lambda, int *lambdalen, double *h, double *x, double *delta2);
    int i;
    double *xx;
    xx = new double[1];
    xx[0] = ((double *) ex)[0];
    int *lambdalen;
    lambdalen = new int[1];
    lambdalen[0] = (int) (((double *) ex)[1]);
    double *lambda;
    lambda = new double[lambdalen[0]];
    for (i = 1; i <= lambdalen[0]; i = i + 1)
    { lambda[i - 1] = ((double *) ex)[i + 1]; }
    double *h;
    h = new double[lambdalen[0]];
    for (i = 1; i <= lambdalen[0]; i = i + 1)
    { h[i - 1] = ((double *) ex)[lambdalen[0] + i + 1]; }
    double *delta2;
    delta2 = new double[lambdalen[0]];
    for (i = 1; i <= lambdalen[0]; i = i + 1)
    { delta2[i - 1] = ((double *) ex)[2 * lambdalen[0] + i + 1]; }
    double *u;
    u = new double[1];
    for (i = 1; i <= n; i = i + 1)
    {
        u[0] = x[i - 1];
        x[i - 1] = imhoffunc(u, lambda, lambdalen, h, xx, delta2);
    }

    delete[] xx;
    delete[] lambdalen;
    delete[] lambda;
    delete[] h;
    delete[] delta2;
    delete[] u;


}


void probQsupx(double *x, double *lambda, int *lambdalen, double *h, double *delta2, double *Qx, double *epsabs,
               double *epsrel, int *limit)
{
    int i;
    void f(double *x, int n, void *ex);
    void Rdqagi(integr_fn f, void *ex, double *bound, int *inf,
                double *epsabs, double *epsrel,
                double *result, double *abserr, int *neval, int *ier,
                int *limit, int *lenw, int *last,
                int *iwork, double *work);
    double *ex;
    ex = new double[2 + 3 * lambdalen[0]];
    ex[0] = x[0];
    ex[1] = (double) lambdalen[0];
    for (i = 1; i <= lambdalen[0]; i = i + 1)
    { ex[i + 1] = lambda[i - 1]; }
    for (i = 1; i <= lambdalen[0]; i = i + 1)
    { ex[lambdalen[0] + i + 1] = h[i - 1]; }
    for (i = 1; i <= lambdalen[0]; i = i + 1)
    { ex[2 * lambdalen[0] + i + 1] = delta2[i - 1]; }
    double *bound;
    bound = new double[1];
    bound[0] = 0.0;
    int *inf;
    inf = new int[1];
    inf[0] = 1;
    // double *epsabs;  
    // epsabs = new double[1];
    // *(epsabs+0) = 0.000001; //0.0001220703;
    // double *epsrel;  
    // epsrel = new double[1];
    // *(epsrel+0) = 0.000001; //0.0001220703;
    double *result;
    result = new double[1];
    double *abserr;
    abserr = new double[1];
    int *neval;
    neval = new int[1];
    int *ier;
    ier = new int[1];
    // int *limit;  
    //   limit = new int[1];
    // *(limit+0) = 10000;
    int *lenw;
    lenw = new int[1];
    *(lenw + 0) = 4 * *(limit + 0);
    int *last;
    last = new int[1];
    int *iwork;
    iwork = new int[*(limit + 0)];
    double *work;
    work = new double[*(lenw + 0)];
    Rdqagi(f, ex, bound, inf, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);
    Qx[0] = 0.5 + result[0] / M_PI;

    epsabs[0] = abserr[0];

    delete[] ex;
    delete[] bound;
    delete[] inf;
    //    delete[] epsabs;
    //    delete[] epsrel;
    delete[] result;
    delete[] abserr;
    delete[] neval;
    delete[] ier;
    //    delete[] limit;
    delete[] lenw;
    delete[] last;
    delete[] iwork;
    delete[] work;
    return;
}

}


