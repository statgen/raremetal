#ifndef __MVTNORM_H__
#define __MVTNORM_H__

#include <cmath>
#include <iostream>
#include <stdio.h>
#include "MathMatrix.h"
#include "MathSVD.h"

#define INF 10000.0

#include <Eigen/Eigenvalues>

extern "C"
{
void mvtdst_(int *N, int *df, double lower[], double upper[], int infin[], double correl[], double delta[], int *maxpts,
             double *abseps, double *releps, double *error, double *value, int *inform);

}

void cov2cor(Matrix &sigma, Matrix &corr)
{
    for (int i = 0; i < sigma.rows; i++)
    {
        for (int j = i; j < sigma.cols; j++)
        {
            if (i == j)
            {
                corr[i][j] = 1.0;
            } else
            {
                corr[i][j] = corr[j][i]
                        = sigma[i][j] / (sqrt(sigma[i][i]) * sqrt(sigma[j][j]));
            }
        }
    }
}

void pmvnorm(double lower[], double upper[], double mean[], Matrix &sigma, bool print, Vector &result)
{
    result.Dimension(3);
    int n = sigma.rows;
    int *infin = new int[n];
    double *delta = new double[n];
    double *corrF = new double[n * (n - 1) / 2];

    Matrix corr;
    corr.Dimension(n, n);

    //reformat the var-cov matrix into correlation matrix
    cov2cor(sigma, corr);

    if (print)
    {
        printf("mean is:\n");
        for (int i = 0; i < n; i++)
        {
            printf("%g\t", mean[i]);
        }
        printf("\n");

        printf("Sigma matrix is:\n");
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                printf("%g\t", sigma[i][j]);
            }
            printf("\n");
        }
    }
/*
      printf("corr matrix is:\n");
      for(int i=0;i<n;i++)
      {
      for(int j=0;j<n;j++)
      printf("%g\t",corr[i][j]);
      printf("\n");
      }
  */
    //fill in the corrF array to pass to Fortran
    for (int j = 0; j < n; j++)
    {
        for (int i = j + 1; i < n; i++)
        {
            int k = j + 1 + ((i - 1) * i) / 2 - 1;
            corrF[k] = corr[i][j];
        }
    }
    /*
       printf("corrF vector is:\n");
       for(int i=0;i<n*(n-1)/2;i++)
       printf("%g\t",corrF[i]);
       printf("\n");
     */
    //recalculate lower and upper limit
    for (int i = 0; i < n; i++)
    {
        delta[i] = 0.0;

        if (lower[i] != INF && lower[i] != -INF)
        {
            lower[i] = (lower[i] - mean[i]) / sqrt(sigma[i][i]);
        }

        if (upper[i] != INF && upper[i] != -INF)
        {
            upper[i] = (upper[i] - mean[i]) / sqrt(sigma[i][i]);
        }

        if (lower[i] == -INF)
        {
            infin[i] = 0;
        }
        if (upper[i] == INF)
        {
            infin[i] = 1;
        }
        if (lower[i] == -INF && upper[i] == INF)
        {
            infin[i] = -1;
        }
        if (lower[i] != -INF && upper[i] != INF)
        {
            infin[i] = 2;
        }
        if (lower[i] == -INF)
        {
            lower[i] = 0;
        }
        if (upper[i] == INF)
        {
            upper[i] = 0;
        }
    }

    int inform = 0;
    double value = 0.0;
    double error = 0.0;

    int sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += infin[i];
    }
    if (sum == -n)
    {
        inform = 0;
        value = 1.0;
    } else
    {
        int df = 0;
        int maxpts = 50000;
        double abseps = 0.001;
        double releps = 0.0;
        mvtdst_(&n, &df, lower, upper, infin, corrF, delta, &maxpts, &abseps, &releps, &error, &value, &inform);
    }
/*
      if(inform >1)
      {
         //      printf("inform=%d. error=%g.value=%g\n",inform,error,value);
         value = -1.0;
      }
*/
    //printf("inform is:%d\n",inform);
    if (inform == 3)
    {
        //printf("reconstructing correlation matrix to make it positive definite using Brissette et al. (2007) method ... ");

        Eigen::MatrixXd newcorr(n, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                newcorr(i, j) = corr[i][j];
            }
        }

        int trial = 0;
        while (inform > 1 && trial < 100)
        {
            Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> es(newcorr);
            Eigen::VectorXd S = es.eigenvalues();
            for (int i = 0; i < n; i++)
            {
                if (S[i] < 0.0)
                {
                    S[i] = 0.0;
                }
            }
            Eigen::MatrixXd D = S.asDiagonal();
            Eigen::MatrixXd V = es.eigenvectors();
            newcorr = V * D * V.transpose();
            //normalize the new matrix
            Eigen::MatrixXd newcorr_diag(n, 1);
            for (int i = 0; i < n; i++)
            {
                newcorr_diag(i, 0) = newcorr(i, i);
            }
            Eigen::MatrixXd norm = newcorr_diag * newcorr_diag.transpose();

            for (int j = 0; j < n; j++)
            {
                for (int i = j + 1; i < n; i++)
                {
                    int k = j + 1 + ((i - 1) * i) / 2 - 1;
                    corrF[k] = newcorr(i, j) / sqrt(norm(i, j));
                }
            }

            int df = 0;
            int maxpts = 50000;
            double abseps = 0.001;
            double releps = 0.0;
            mvtdst_(&n, &df, lower, upper, infin, corrF, delta, &maxpts, &abseps, &releps, &error, &value, &inform);
            trial++;
        }

        if (inform > 1)
        {
            value = -1.0;
        }
    }

    result[0] = value;
    result[1] = error;
    result[2] = inform;

    if (infin)
    { delete[] infin; }
    if (corrF)
    { delete[] corrF; }
    if (delta)
    { delete[] delta; }

    //   return result;
}

#endif
