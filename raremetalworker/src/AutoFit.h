////////////////////////////////////////////////////////////////////// 
// PreMeta.h 
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

#ifndef __AUTOFIT_H__
#define __AUTOFIT_H__

#include "MathNormal.h"
#include "Pedigree.h"
#include "MathSVD.h"
#include "MathMatrix.h"

class NormalEquations;

class AutoFit
{
public:
    AutoFit();

    NormalSet mvn;

    void FitPolygenicModels(Pedigree &ped, int traitNum, bool quiet = false);

    static bool fitSharedEnvironment;
    static bool fitX;
    static bool useCovariates;
    static bool useProbands;

    double mean, sigma_gX, sigma_e2, variance, heritability, heritabilityX, sharedEnv;
    Vector betaHat;
};

#endif


