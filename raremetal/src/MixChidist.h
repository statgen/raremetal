#ifndef __MIXCHIDIST_H__
#define __MIXCHIDIST_H__

#include <cstring>
#include <iostream>
#include <stdlib.h>
#include "StringBasics.h"
#include "qfc.h"
#include <algorithm>
#include "math.h"

#define MATHLIB_STANDALONE

#include <Rmath.h>
//#include "imhof.h"

double MixChidist(double lambda[], int n, double q_obs, String method)
{
    double *nc1 = new double[n];
    int *n1 = new int[n];
    double *trace = new double[7];
    double sigma = 0.0;
    int lim1 = 10000;
    double acc = 0.0001;
    int ifault = 0;
    double res = 0.0;

    for (int i = 0; i < n; i++)
    {
        nc1[i] = 0.0; //noncentral parameters
        n1[i] = 1; //h in R function
    }
    for (int i = 0; i < 7; i++)
    {
        trace[i] = 0.0;
    }

    if (method == "Davies")
    {
        //probability will be saved in res.
        qfc(lambda, nc1, n1, &n, &sigma, &q_obs, &lim1, &acc, trace, &ifault, &res);
        res = 1.0 - res;
    }

    //From Lee's SKAT package
    if (method == "Liu")
    {
        double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0;
        for (int i = 0; i < n; i++)
        {
            c1 += lambda[i];
            c2 += lambda[i] * lambda[i];
            c3 += lambda[i] * lambda[i] * lambda[i];
            c4 += lambda[i] * lambda[i] * lambda[i] * lambda[i];
        }

//printf("%g,%g,%g,%g\n",c1,c2,c3,c4);
        double s1 = c3 / sqrt(c2 * c2 * c2);
        double s2 = c4 / (c2 * c2);
        double muQ = c1;
        double sigmaQ = sqrt(2.0 * c2);

        double tstar = (q_obs - muQ) / sigmaQ;
        double delta, l, a;

        if (s1 * s1 > s2)
        {
            a = 1.0 / (s1 - sqrt(s1 * s1 - s2));
            delta = s1 * a * a * a - a * a;
            l = a * a - 2.0 * delta;
        } else
        {
            a = 1.0 / s1;
            delta = 0.0;
            l = c2 * c2 * c2 / (c3 * c3);
        }

        double muX = l + delta;
        double sigmaX = sqrt(2.0) * a;
        double q_new = tstar * sigmaX + muX;

        //printf("q_new=%g, df=%g, ncp=%g\n",q_new,l,delta);
        double Qq;
        if (delta == 0)
        {
            Qq = pchisq(q_new, l, 0, 0);
        } else
        {
            Qq = pnchisq(q_new, l, delta, 0, 0);
        }

        res = Qq;
    }

    /*
       if(method=="Liu")
       {
       double c1 = 0.0;
       double c2 = 0.0;
       double c3 = 0.0;
       double c4 = 0.0;
       for(int i=0;i<n;i++)
       {
       c1 += n1[i] * lambda[i];
       c2 += n1[i] * lambda[i]*lambda[i];
       c3 += n1[i] * lambda[i]*lambda[i]*lambda[i];
       c4 += n1[i] * lambda[i]*lambda[i]*lambda[i]*lambda[i];
       }
       for(int i=0;i<n;i++)
       {
       c1 += nc1[i] * lambda[i];
       c2 += 2.0*nc1[i] * lambda[i]*lambda[i];
       c3 += 3.0*nc1[i] * lambda[i]*lambda[i]*lambda[i];
       c4 += 4.0*nc1[i] * lambda[i]*lambda[i]*lambda[i]*lambda[i];
       }
       double s1 = c3/(sqrt(c2*c2*c2));
       double s2 = c4/(c2*c2);
       double muQ = c1;
       double sigmaQ = sqrt(2.0*c2);
       double tstar = (q_obs-muQ)/sigmaQ;
       double delta,l,a,Qq;

       if(s1*s1>s2)
       {
    //      a = q_obs/(s1-sqrt(s1*s1-s2));
    a = 1.0/(s1-sqrt(s1*s1-s2));
    delta = s1*a*a*a - a*a;
    l = a*a-2.0*delta;
    }
    else
    {
    //a = 1.0/s1;
    l = 1.0/s2;
    a=sqrt(l);
    delta = 0.0;
    //l = c2*c2*c2/(c3*c3);
    }
    double muX = l + delta;
    double sigmaX = sqrt(2.0)*a;

    double stat = tstar*sigmaX+muX;
    printf("stat is:%g;df is:%g;delta is: %g\n",stat,l,delta);

    if(stat<0||l<0||delta<0)
    return -1.0;
    if(delta==0.0)
    Qq = chidist(stat,l);
    else
    Qq = chidist(stat,l,delta);

    return Qq;
    }
     */
    return res;
}

#endif
