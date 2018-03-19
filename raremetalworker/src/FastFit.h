///////////////////////////////////////////////////////////////////// 
// FastFit.h 
// (c) 2012-2013 Shuang Feng, Dajiang Liu, Goncalo Abecasis
// 
// This file is distributed as part of the RareMetalWorker source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile RareMetalWorker.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Wednesday November 28, 2012
// 

#ifndef __FASTFIT_H__
#define __FASTFIT_H__

#include "MathGold.h"
#include "TransformResiduals.h"
#include "MathCholesky.h"
#include "TraitTransformations.h"
#include "MathConstant.h"
#include "MathSVD.h"

//the same definition can be found in MathGold.cpp
#define GLIMIT             100      // Maximum magnification, per round
#define SHIFT(a, b, c, d)  (a)=(b); (b)=(c); (c)=(d)
#define SHFT3(a, b, c)     (a)=(b); (b)=(c)


//This class uses optimization method proposed in FaST paper by Lippert 
//from Nature Vol.8 No.10 833-835 2011
// Factored spectrally transformed linear mixed models (FaST-LMM)
class FastFit : public ScalarMinimizer
{
public:
    FastFit()
    {};

    ~FastFit()
    {};

    virtual double Brent(double tol, FastTransform &trans);

    static bool useCovariates;
    /**
     * If --inverseNormal is used, but not with --makeResiduals, then trait values are inverse normalized before
     *  fitting linear models.
     * If --inverseNormal and --makeResiduals are used together, then covariates are adjusted and inverse normalized
     *  residuals are used to fit linear models.
     */
    static bool inverseNormal;
    /**
     * Whether to adjust covariates before fitting linear models using residuals.
     */
    static bool makeResiduals;
    /**
     * A flag that indicates whether this program should ever attempt to account for relatedness in the sample. Determined indirectly based on other command line arguments and defaults to false.
     */
    static bool unrelated;
    /**
     *  Whether to fit a linear mixed model using only chromosome X kinship for analyses of chromosome X markers.
     *  To use this flag, the vcX option must also be specified
     */
    static bool separateX;
    static bool CleanKin;
    /**
     * Undocumented CLI parameter: whether to treat this trait as binary. Cannot be used with --makeResiduals flag
     */
    static bool binary;
    int traitNum;
    /**
     * The name of the trait to be analyzed. If not specified, then all traits included in PED/DAT files are analyzed.
     */
    static String traitName;

    int GetTraitID(Pedigree &ped, const char *name, FILE *log);

    //Estimated parameters are saved here
    double sigma_g2Hat, sigma_e2Hat, deltaHat, sigma_gXHat, sigma_s2Hat, logLikelihood, sigma2Hat;
    double heritability, Xheritability, heritabilityAuto;
    Vector beta;
    Vector betaHat;
    Vector residuals_unrelated;
    Matrix inv;

    //Estimation functions
    double Evaluate(double delta, FastTransform &trans, bool full);

    double EvaluateMulti(FastTransform &trans, bool full);

    void GetBeta(double delta, FastTransform &trans);

    double GetSigma_g(double delta, Vector &beta, FastTransform &trans);

    //This is where data is transformed and parameters are estimated
    void PreFit(IFILE SCOREoutput, IFILE SCOREoutput_rec, IFILE SCOREoutput_dom, Pedigree &ped, double tol,
                FastTransform &trans, FILE *log);

    void FitModels(IFILE SCOREoutput, IFILE SCOREoutput_rec, IFILE SCOREoutput_dom, Pedigree &ped, double tol,
                   FastTransform &trans, KinshipEmp &kin_emp, FILE *log);

    void FitSimpleLinearModels(Pedigree &ped, double tol, FastTransform &trans, KinshipEmp &kin_emp, FILE *log);

    void FitSimpleLogisticModels(Pedigree &ped, FastTransform &trans, FILE *log);

    void
    AutoFitLinearMixModels(Pedigree &ped, double tol, FastTransform &trans, KinshipEmp &kin_emp, FILE *log, bool forX);

    void
    AutoFitLinearMixModels2(Pedigree &ped, double tol, FastTransform &trans, KinshipEmp &kin_emp, FILE *log, bool forX);

    void
    FastFitPolyGenicModels(Pedigree &ped, double tol, FastTransform &trans, KinshipEmp &kin_emp, FILE *log, bool forX);

    void PrintPreface(IFILE SCOREoutput, Pedigree &ped, FastTransform &trans, FILE *log);

    //optimization
    void OptimizeDelta(double tol, FastTransform &trans);
};

#endif


