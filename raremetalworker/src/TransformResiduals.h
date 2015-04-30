////////////////////////////////////////////////////////////////////// 
// TransformResiduals.h 
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

#ifndef __TRANSFORMRESIDULAS_H__
#define __TRANSFORMRESIDULAS_H__

#include "Kinship.h"
#include "KinshipX.h"
#include "KinshipEmp.h"
#include "MathNormal.h"
#include "MathVector.h"
#include "StringBasics.h"
#include "MathMatrix.h"
#include "MathCholesky.h"
#include "AutoFit.h"
#include <Eigen/Eigenvalues>

class FastTransform
{
	public:
		FastTransform();
		~FastTransform();

		static bool empKin;
		static bool pedKin;
		static bool mergedVCFID;
		static String readInEmp;
		static String readInEmpX;

		//# of useful families and individuals
		int families,persons,fixedEffects,analyzedFounders,numFounder,totalN,numFamily;
		double traitVar,traitMean;

		//Hash of sample ids from vcf file, if reading from vcf is needed
		StringIntHash sampleVCFIDHash,samplePEDIDHash,foundersHash;

		IntArray * pPheno;
		IntArray genotypedSamplePED, genotypedSampleVCF;
		StringArray founders;
		//Saved Matrices for likelihood calculation
		Eigen::MatrixXf transU; 
		Matrix transU_del;
		Matrix X,UX;
		Vector D,UY,Y; 
		//Matrix UDX,UD,inv; //here D is D^1/2
		//Vector UDY; //D is D^1/2

		//Reformat LL* into UDU*U
		void CleanUp()
		{
		   //X.Dimension(0,0);
		   UX.Dimension(0,0);
		};

		int GetCovarianceMatrix(AutoFit & engine, Pedigree &ped, int f, int count,Matrix & omega);
		void ReformCholesky(Matrix & L, Vector & D);
		int GetTraitID(Pedigree & ped, const char * name);
		//These functions are transformation work horses
		void GetFixedEffectsCount(Pedigree & ped, bool useCovariates);
		void Prepare(Pedigree & ped, int traitNum,bool useCovariates,FILE * log,bool shortVersion);
		void SubSet(Matrix & kin,Pedigree & ped,Matrix & tmp);
		void SelectSamplesPED(Pedigree & ped,bool useCovariates);
		void SelectSamplesVCF(Pedigree & ped, bool useCovariates);
		void ScreenSampleID(Pedigree & ped, bool useCovariates);

		void TransformPedkin(Pedigree & ped, int traitNum,bool useCovaraites,FILE * log);
		void TransformPedkinSepX(Pedigree & ped, int traitNum,bool useCovaraites,FILE * log);
		void TransformEmpkin(Pedigree & ped, int traitNum,bool useCovariates,KinshipEmp & kin_emp,Matrix & input,FILE * log);
		void TransformPedkinX(AutoFit & engine, Pedigree &ped);
		void TransformEmpkinX(Matrix & covMatrix);
		//void FinalizeTransform(double sigma_g2,double sigma_e2);
		int EigenDecompose(Matrix & matrix_in, Matrix & U, Vector & D);
		bool FinalizeProducts();
};

#endif
