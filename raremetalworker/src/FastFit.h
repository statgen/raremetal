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
class FastFit : public ScalarMinimizer
{
	public:
		FastFit() {};
		~FastFit() {};
		virtual double Brent(double tol, FastTransform & trans);
		static bool useCovariates;
		static bool inverseNormal;
		static bool makeResiduals;
		static bool unrelated;
		static bool separateX;
		static bool CleanKin;
		static bool binary;
		int traitNum;
		static String traitName;

		int GetTraitID(Pedigree & ped, const char * name,FILE * log);

		//Estimated parameters are saved here
		double sigma_g2Hat,sigma_e2Hat,deltaHat,sigma_gXHat,sigma_s2Hat,logLikelihood,sigma2Hat;
		double heritability,Xheritability,heritabilityAuto;
		Vector beta;
		Vector betaHat;
		Vector residuals_unrelated;
		Matrix inv;

		//Estimation functions
		double Evaluate(double delta,FastTransform & trans,bool full);
		double EvaluateMulti(FastTransform & trans,bool full);
		void GetBeta(double delta, FastTransform & trans);
		double GetSigma_g(double delta,Vector & beta, FastTransform & trans);

		//This is where data is transformed and parameters are estimated 
		void PreFit(IFILE SCOREoutput,IFILE SCOREoutput_rec,IFILE SCOREoutput_dom,Pedigree & ped,double tol,FastTransform & trans,FILE * log);
		void FitModels(IFILE SCOREoutput,IFILE SCOREoutput_rec,IFILE SCOREoutput_dom,Pedigree & ped, double tol, FastTransform & trans,KinshipEmp & kin_emp, FILE * log);
		void FitSimpleLinearModels(Pedigree & ped, double tol, FastTransform & trans,KinshipEmp & kin_emp,FILE * log);
		void FitSimpleLogisticModels(Pedigree & ped, FastTransform & trans,FILE * log);

		void AutoFitLinearMixModels(Pedigree & ped, double tol, FastTransform & trans,KinshipEmp & kin_emp,FILE * log,bool forX);
		void AutoFitLinearMixModels2(Pedigree & ped, double tol, FastTransform & trans,KinshipEmp & kin_emp,FILE * log,bool forX);
		void FastFitPolyGenicModels(Pedigree & ped, double tol, FastTransform & trans,KinshipEmp & kin_emp,FILE * log, bool forX);

		void PrintPreface(IFILE SCOREoutput,Pedigree & ped,FastTransform & trans,FILE * log);

		//optimization
		void OptimizeDelta(double tol,FastTransform & trans);
};

#endif


