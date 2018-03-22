#include <iostream>
#include <math.h>
#include <stdlib.h> // exit
#include "logistic.h"
#include "MathSVD.h"

// test logistic regression
// use gradient & newton

void GetLogitCoeff(Vector &betaHat, Matrix &X, Vector &Y)
{
    // debug
    printf("\nX size = %d,%d\n", X.rows, X.cols);
    printf("X=");
    for (int i = 0; i < X.rows; i++)
    {
        printf(",%g", X[i][0]);
    }
    printf("\nY length = %d\n", Y.Length());
    printf("Y=");
    for (int i = 0; i < Y.Length(); i++)
    {
        printf(",%g", Y[i]);
    }
    printf("\n");

    int ncoeff = X.cols; // X.dimension = #persons, fixed_effects
    int persons = Y.Length(); // #persons

    double min_dist = 0.0001; // convergence condition
    double gamma = 0.00005; // gradient move step
    int max_iters = 2000;
    int iter = 0;
    bool newton = 0; // if true, use newton's method

    Matrix w_old;
    Matrix w_new;
    w_old.Dimension(ncoeff, 1, 0);
    w_new.Dimension(ncoeff, 1, 0);

    Matrix Hessian;
    if (newton)
    {
        Hessian.Dimension(ncoeff, ncoeff, 0);
    }

    double dist = min_dist + 1;
    Matrix grad1;
    grad1.Dimension(ncoeff, 1, 0);
    while (iter < max_iters)
    {
        for (int k = 0; k < ncoeff; k++)
        {
            for (int i = 0; i < persons; i++)
            {
                double wx = 0;
                for (int j = 0; j < ncoeff; j++)
                {
                    wx += w_old[j][0] * X[i][j];
                }
                double sm = sigmoid(-Y[i] * wx);
                grad1[k][0] += Y[i] * X[i][k] * sm;
                //debug
//				printf("g,sm = %g,%g\n", grad1[k][0],sm);
                if (newton)
                {
                    for (int j = 0; j < ncoeff; j++)
                    {
                        Hessian[k][j] -= X[i][j] * X[i][k] * sm * (1 - sm);
                    }
                }
            }
        }
        // debug
        printf("w_old = ");
        for (int x = 0; x < ncoeff; x++)
        {
            printf("%g ", w_old[x][0]);
        }
        printf("\n");
        if (newton)
        {
            double d = GetDet(Hessian.rows, Hessian);
            if (d == 0)
            {
                newton = 0;
            } // use gradiant instead
        }
        // update coefficients
        if (newton)
        {
            SVD svd;
            Matrix invHessian = Hessian;
            svd.InvertInPlace(invHessian);
            invHessian.Product(invHessian, grad1);
            if (invHessian.cols != 1)
            {
                error("invHessian col != 1!");
            }
            w_new = w_old;
            w_new.Add(invHessian);
        } else
        { // gradient method
            grad1.Multiply(gamma);
            w_new = w_old;
            w_new.Add(grad1);
        }
        // debug
        printf("iter=%d\n", iter);
        printf("w_new = ");
        for (int x = 0; x < ncoeff; x++)
        {
            printf("%g ", w_old[x][0]);
        }
        printf("\n");
        dist = NormVector(w_new, w_old);
        printf("iter=%d, dist=%g\n", iter, dist);
        iter++;
        if (dist <= min_dist)
        {
            break;
        } else
        {
            w_old = w_new;
        }
    }

    if (dist > min_dist)
    {
        printf("Warning: logistic regression does not converge!\n");
    }
    for (int i = 0; i < ncoeff; i++)
    {
        betaHat[i] = w_new[i][0];
    }

    // debug
    //exit(0);
}

double sigmoid(double x)
{
    double sm = 1 / (1 + exp(-x));
    return sm;
}

// convert Y to +1 and -1. Record the original value.
void convertBinaryTrait(Vector &Y, int &Yplus1, int &Yminus1)
{
    bool see2;
    double v1 = Y[0]; // +1
    double v2;
    if (Y.Length() < 2)
    {
        error("less than 2 samples have phenotype!\n");
    }
    for (int i = 1; i < Y.Length(); i++)
    {
        if (!see2)
        {
            if (Y[i] == v1)
            {
                continue;
            }
            if (Y[i] > v1)
            {
                v2 = v1;
                v1 = Y[i];
            } else
            {
                v2 = Y[i];
            }
            see2 = 1;
            continue;
        }
        if (Y[i] != v1 && Y[i] != v2)
        {
            std::cout << "\nERROR: 3rd phenotype " << Y[i] << " exists. Is that binary trait?" << std::endl;
            std::cerr << "\nERROR: 3rd phenotype " << Y[i] << " exists. Is that binary trait?" << std::endl;
            exit(1);
        }
    }
    // convert
    Yplus1 = v1;
    Yminus1 = v2;
    if (v1 == 1 && v2 == -1)
    {
        return;
    }
    for (int i = 0; i < Y.Length(); i++)
    {
        if (Y[i] == v1)
        {
            Y[i] = 1;
        } else
        {
            Y[i] = -1;
        }
    }
}

double NormVector(Matrix &v1, Matrix &v2)
{
    if (v1.rows != v2.rows || v1.cols * v2.cols != 1)
    {
        error("NormVector: v1 & v2 do not have same length or have >2 cols!");
    }

    double dist = 0;
    for (int i = 0; i < v1.rows; i++)
    {
        double x = v1[i][0] - v2[i][0];
        if (x < 0)
        {
            x *= -1;
        }
        dist += x;
    }
    return dist;
}

double GetDet(int n, Matrix &m)
{
    if (m.rows != m.cols)
    {
        error("Non-square matrix cannot be inverted!");
    }

    if (n == 1)
    {
        return m[0][0];
    }

    if (n == 2)
    {
        return (m[0][0] * m[1][1] - m[1][0] * m[0][1]);
    }

    double d = 0;
    Matrix submat;
    submat.Dimension(n, n);
    for (int c = 0; c < n; c++)
    {
        int subi = 0;
        for (int i = 1; i < n; i++)
        {
            int subj = 0;
            for (int j = 0; j < n; j++)
            {
                if (j == c)
                {
                    continue;
                }
                submat[subi][subj] = m[i][j];
                subj++;
            }
            subi++;
        }
        d = d + pow(-1, c) * m[0][c] * GetDet(n - 1, submat);
    }
    return d;
}

bool isBinaryY(Vector &Y)
{
    if (Y.Length() < 2)
    {
        return 1;
    }
    bool see2;
    double v1 = Y[0]; // +1
    double v2;
    for (int i = 0; i < Y.Length(); i++)
    {
        if (!see2)
        {
            if (Y[i] == v1)
            {
                continue;
            }
            v2 = v1;
            v1 = Y[i];
            see2 = 1;
            continue;
        }
        if (Y[i] != v1 && Y[i] != v2)
        {
            return 0;
        }
    }
    // only 2 values
    return 1;
}


