#ifndef __AUTOFIT2_H__
#define __AUTOFIT2_H__

#include "MathNormal.h"
#include "Pedigree.h"
#include "TransformResiduals.h"

class NormalEquations;

class AutoFit2
{
public:
    AutoFit2();

    void FitPolygenicModels(Pedigree &ped, KinshipEmp &kin_emp, FastTransform &trans, int traitNum, bool quiet = false);

    static bool useCovariates;

    Matrix covMatrix;

    //This function write covariance matrix to covMatrix after VCs have been estimated

    double mean, variance, sigma_e2, heritability, heritabilityX;
    Vector betaHat;
};

#endif


