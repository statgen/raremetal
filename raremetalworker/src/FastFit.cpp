////////////////////////////////////////////////////////////////////// 
// FastFit.cpp 
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

#include "FastFit.h"
#include "PreMeta.h"
#include "AutoFit2.h"
#include "logistic.h"

int traitNum = 0;
bool FastFit::useCovariates = false;
bool FastFit::inverseNormal = false;
bool FastFit::makeResiduals = false;
bool FastFit::unrelated = false;
bool FastFit::separateX = false;
bool FastFit::CleanKin = false;
bool FastFit::binary = false;
String FastFit::traitName = "";

void FastFit::OptimizeDelta(double tol,FastTransform & trans)
{
   double inc = 0.2;
   for(int i=0;i<100;i++)
   {
      //printf("Start optimization ...\n");
      a = exp(-10.0+i*inc); fa=-1.0*Evaluate(a,trans,false);
      c = exp(-10.0+(i+1)*inc); fc=-1.0*Evaluate(c,trans,false);
      b = (a+c)/2.0; fb=-1.0*Evaluate(b,trans,false);

      double ll = -Brent(tol,trans);
      if(i==0) 
      {
	 logLikelihood = ll;
	 deltaHat = min;
      }
      if(ll>logLikelihood) 
      {
	 logLikelihood = ll;
	 deltaHat = min;
      }
      //printf("logLikelihood is:%g,deltaHat is:%g\n",logLikelihood,deltaHat);
   }
}

void FastFit::GetBeta(double delta,FastTransform & trans)
{
   Matrix A;
   Vector B;
   A.Dimension(trans.fixedEffects,trans.fixedEffects,0.0);
   B.Dimension(trans.fixedEffects);
   B.Zero();

   //Now fill in matrix A
   Matrix subUXXU;
   subUXXU.Dimension(trans.fixedEffects,trans.fixedEffects);
   subUXXU.Zero();
   //fill in the subUXXU matrix
   for(int i=0;i<trans.persons;i++)
   {
      for(int r=0;r<trans.fixedEffects;r++) 
      {
	 for(int c=0;c<trans.fixedEffects;c++)
	 {
	    if(r==c) {
	       subUXXU[r][c] = trans.UX[i][c]*trans.UX[i][c];
	       subUXXU[r][c]= subUXXU[r][c]/(trans.D[i]+delta);
	    } 
	    else 
	    {
	       subUXXU[r][c]= trans.UX[i][r]*trans.UX[i][c];
	       subUXXU[r][c]= subUXXU[r][c]/(trans.D[i]+delta);
	       subUXXU[c][r]=subUXXU[r][c];
	    }
	 }
      }
      A.Add(subUXXU);
   }

   SVD svd; 
   svd.InvertInPlace(A);

   //Now fill in vector B
   Vector tmpAddToA;
   for(int i=0;i<trans.persons;i++)
   {
      tmpAddToA.Copy(trans.UX[i]);
      double factor = trans.UY[i]/(trans.D[i]+delta);
      tmpAddToA.Multiply(factor);
      B.Add(tmpAddToA);
   }
   beta.Product(A,B);
}

double FastFit::GetSigma_g(double delta,Vector & beta,FastTransform & trans)
{
   double sigma =0.0;
   for(int i=0;i<trans.persons;i++)
   {
      double factor = trans.UY[i] - trans.UX[i].InnerProduct(beta);
      sigma += factor*factor/(trans.D[i]+delta);
   }
   sigma /=trans.persons;
   return sigma;
}

double FastFit::Evaluate(double delta,FastTransform & trans,bool full)
{
   double logLikelihood;
   //Fill up beta vector
   GetBeta(delta,trans);	
   double sigma = GetSigma_g(delta, beta,trans);
   //calculate the constant
   double constant =0.0;
   for(int i=0;i<trans.persons;i++)
      constant += log(trans.D[i] + delta);
   if(full)
      logLikelihood = -0.5*(trans.persons*log(2*M_PI)+constant + trans.persons + trans.persons*log(sigma));
   else logLikelihood = -0.5*(constant + trans.persons + trans.persons*log(sigma));
   return logLikelihood;
}

double FastFit::EvaluateMulti(FastTransform & trans,bool full)
{
   double logLikelihood;
   double constant =0.0;
   for(int i=0;i<trans.persons;i++)
      constant += log(trans.D[i]);

   Vector residuals,fixedEff,weightedRes;

   fixedEff.Product(trans.UX,betaHat);
   residuals.Copy(trans.UY);
   residuals.Subtract(fixedEff);
   weightedRes.Copy(residuals);

   for(int i=0;i<trans.persons;i++)
      weightedRes[i] /= trans.D[i];

   if(full)
      logLikelihood = -0.5*(trans.persons*log(2*M_PI)+ constant + weightedRes.InnerProduct(residuals));
   else
   {
      logLikelihood = -0.5*(constant + weightedRes.InnerProduct(residuals));
   } 
   return logLikelihood;
}

void FastFit::PrintPreface(IFILE output,Pedigree & ped,FastTransform & trans,FILE * log)
{
   ifprintf(output,"##ProgramName=RareMetalWorker\n");
   ifprintf(output,"##Version=%s\n",VERSION);
   ifprintf(output,"##Samples=%d\n",trans.totalN);
   ifprintf(output,"##AnalyzedSamples=%d\n",trans.persons);
   ifprintf(output,"##Families=%d\n",trans.numFamily);
   ifprintf(output,"##AnalyzedFamilies=%d\n",trans.families);
   ifprintf(output,"##Founders=%d\n",trans.numFounder);
   if(makeResiduals)
      ifprintf(output, "##MakeResiduals=True\n");
   else
      ifprintf(output, "##MakeResiduals=False\n");
   ifprintf(output,"##AnalyzedFounders=%d\n",trans.analyzedFounders);

   if((useCovariates || makeResiduals) && ped.covariateNames.Length()>0)
   {
      ifprintf(output,"##Covariates=");
      ifprintf(output,"%s",ped.covariateNames[0].c_str());
      for(int i=1;i<ped.covariateCount;i++)
         ifprintf(output,",%s",ped.covariateNames[i].c_str());
      ifprintf(output,"\n");
      ifprintf(output,"##CovariateSummaries\tmin\t25th\tmedian\t75th\tmax\tmean\tvariance\n");
      for(int i=0;i<ped.covariateCount;i++) {
         Vector covariate;
         for (int f = 0; f < ped.familyCount; f++) {
            for(int j=0;j<trans.pPheno[f].Length();j++)
               covariate.Push(ped[trans.pPheno[f][j]].covariates[i]);
         }
        covariate.Sort();
         //covariate.RemoveDuplicates();
         int first,second,middle;
         first = (int) covariate.Length()/4;
         second = first*3;
         middle = (int) (first+second)/2;
         ifprintf(output,"##%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",ped.covariateNames[i].c_str(),covariate.Min(),covariate[first],covariate[middle],covariate[second],covariate.Max(),covariate.Average(),covariate.Var());
      }
   }
   else 
   {
      ifprintf(output,"##Covariates=\n");
   }
   ifprintf(output,"##InverseNormal=%s\n",inverseNormal?"ON":"OFF");
   int traitNum = GetTraitID(ped,traitName.c_str(),log);
   ifprintf(output,"##TraitSummaries\tmin\t25th\tmedian\t75th\tmax\tmean\tvariance\n");
   Vector y;
   for (int f = 0; f < ped.familyCount; f++) 
   {    
      for(int j=0;j<trans.pPheno[f].Length();j++)
      {
	 y.Push(ped[trans.pPheno[f][j]].traits[traitNum]);
      }         
   }
   y.Sort();
   int first,second,middle;
   first = (int) y.Length()/4;
   second = first*3;
   middle = (int) (first+second)/2;
   ifprintf(output,"##%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",traitName.c_str(),y.Min(),y[first],y[middle],y[second],y.Max(),y.Average(),y.Var());
}

void FastFit::PreFit(IFILE SCOREoutput,IFILE SCOREoutput_rec,IFILE SCOREoutput_dom,Pedigree & ped,double tol,FastTransform & trans,FILE * log)
{
   traitNum = GetTraitID(ped,traitName.c_str(),log);
   trans.GetFixedEffectsCount(ped,useCovariates);
   betaHat.Dimension(trans.fixedEffects);
   trans.Prepare(ped,traitNum,useCovariates,log,false);

   PrintPreface(SCOREoutput,ped,trans,log);
   if(PreMeta::recessive)
      PrintPreface(SCOREoutput_rec,ped,trans,log);
   if(PreMeta::dominant)
      PrintPreface(SCOREoutput_dom,ped,trans,log);

   if(makeResiduals) {
      FastFit::useCovariates=true;
      Matrix inv;
      Matrix transX;
      transX.Transpose(trans.X);

      //Get beta_hat and residuals
      Vector tmp;
      tmp.Dimension(trans.fixedEffects);
      for(int j=0;j<trans.fixedEffects;j++)
         tmp[j] =trans.Y.InnerProduct(transX[j]);

      inv.Product(transX,trans.X);
      SVD svd;
      svd.InvertInPlace(inv);

      for(int i=0;i<trans.fixedEffects;i++)
         betaHat[i] = inv[i].InnerProduct(tmp);

      residuals_unrelated.Dimension(trans.persons);
      for(int i=0;i<trans.persons;i++)
         residuals_unrelated[i] = trans.Y[i] - betaHat.InnerProduct(trans.X[i]);
      //replace trait values into residuals
      int index=0;
      for(int i=0;i<ped.count;i++)
         ped[i].traits[traitNum] = _NAN_;

      for (int f = 0; f < ped.familyCount; f++) {
         for(int j=0;j<trans.pPheno[f].Length();j++) {
            ped[trans.pPheno[f][j]].traits[traitNum] = residuals_unrelated[index];
            index++;
         }
      }
      FastFit::useCovariates=false;
   }

   if(inverseNormal)
      InverseNormalTransform(ped,traitNum);

   trans.GetFixedEffectsCount(ped,useCovariates);
   beta.Clear();
   betaHat.Clear();
   beta.Dimension(trans.fixedEffects);
   betaHat.Dimension(trans.fixedEffects);
}

void FastFit::FitModels(IFILE SCOREoutput,IFILE SCOREoutput_rec,IFILE SCOREoutput_dom,Pedigree & ped,double tol,FastTransform & trans,KinshipEmp & kin_emp,FILE * log)
{
   heritability=0.0;
   heritabilityAuto=0.0;
   Xheritability = 0.0;
   //if N>5000, choose to fit one VC based on kinshipX only.
   if(!separateX && FastTransform::empKin && AutoFit::fitX && trans.persons>10000)
   {
      printf("Warning: sample size is >10000. --separateX option is activated.\n");
      separateX=true;
   }

   if(makeResiduals)
      FastFit::useCovariates=true;

   PreFit(SCOREoutput,SCOREoutput_rec,SCOREoutput_dom,ped, tol, trans,log);

   if(unrelated) {
      printf("  Start fitting linear model ... ");
      fprintf(log,"  Start fitting linear model ... ");
      fflush(stdout);
   }
   else
   {
      printf("  Start fitting linear mixed model ... ");
      fprintf(log,"  Start fitting linear mixed model ... ");
      fflush(stdout);
   }
   if(unrelated)
   {
      if (binary)
         FitSimpleLogisticModels(ped, trans, log);
      else
         FitSimpleLinearModels(ped,tol,trans,kin_emp,log);
   }
   else {
      //when fitting the basic variance component model, use fast fit algorithm
      FastFitPolyGenicModels(ped,tol,trans,kin_emp,log,false);
   }

   //Calculate Cov(beta_hat)
   //double sigma_g2 = engine.sigma_g2Hat;
   //double sigma_e2 = engine.sigma_e2Hat;
   if(!unrelated)
   {
	  Vector sigma2(trans.persons);
	  for(int i=0;i<trans.persons;i++)
         sigma2[i] = sigma_g2Hat*trans.D[i]+sigma_e2Hat;

	 inv.Dimension(trans.UX.cols,trans.UX.cols);
	 Eigen::MatrixXd X;
	 Eigen::MatrixXd left;
	 Eigen::MatrixXd tmp;
	 Eigen::MatrixXd inv_;
	 X.resize(trans.UX.rows,trans.UX.cols);
	 left.resize(trans.UX.cols,trans.UX.rows);
	 inv_.resize(trans.UX.cols,trans.UX.cols);
	 for(int i=0;i<trans.UX.rows;i++)
	    for(int j=0;j<trans.UX.cols;j++)
	       X(i,j)=trans.UX[i][j];
	 left = X.transpose();

	 for(int i=0;i<left.rows();i++)
	    for(int j=0;j<left.cols();j++)
	       left(i,j) /= sigma2[j]; 
	 tmp = left*X;
	 inv_ = tmp.inverse();
	 for(int i=0;i<inv.rows;i++)
	    for(int j=0;j<inv.cols;j++)
	       inv[i][j] = inv_(i,j);
   }
   else {
      inv.Dimension(trans.X.cols,trans.X.cols);
      Eigen::MatrixXd X;
      Eigen::MatrixXd tmp;
      Eigen::MatrixXd inv_;
      X.resize(trans.X.rows,trans.X.cols);
      for(int i=0;i<trans.X.rows;i++)
         for(int j=0;j<trans.X.cols;j++)
            X(i,j)=trans.X[i][j];
      tmp.resize(trans.persons,trans.persons);
      inv_.resize(trans.X.cols,trans.X.cols);
      tmp = X.transpose()*X;
      inv_ = tmp.inverse(); 
      inv_ *= sigma2Hat;
      for(int i=0;i<inv.rows;i++)
       for(int j=0;j<inv.cols;j++)
          inv[i][j] = inv_(i,j);
   }
   ifprintf(SCOREoutput,"## - NullModelEstimates\n");
   ifprintf(SCOREoutput,"## - Name\tBetaHat\tSE(BetaHat)\n");
   ifprintf(SCOREoutput,"## - Intercept\t%g\t%g\n",betaHat[0],sqrt(inv[0][0]));
   if(betaHat.Length()>1) {
      for(int b=1;b<betaHat.Length();b++)
         ifprintf(SCOREoutput,"## - %s\t%g\t%g\n",ped.covariateNames[b-1].c_str(),betaHat[b],sqrt(inv[b][b]));
   }

   Vector y;
   for (int f = 0; f < ped.familyCount; f++) {
      for(int j=0;j<trans.pPheno[f].Length();j++)
         y.Push(ped[trans.pPheno[f][j]].traits[traitNum]);
   }
   y.Sort();
   int first,second,middle;
   first = (int) y.Length()/4;
   second = first*3;
   middle = (int) (first+second)/2;
   ifprintf(SCOREoutput,"##AnalyzedTrait\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",y.Min(),y[first],y[middle],y[second],y.Max(),y.Average(),y.Var());

   if(unrelated)
      ifprintf(SCOREoutput,"##Sigma_e2_Hat\t%g\n",sigma2Hat);
   else
   {
      ifprintf(SCOREoutput,"##Sigma_g2_Hat\t%g\n",sigma_g2Hat);
      ifprintf(SCOREoutput,"##Sigma_e2_Hat\t%g\n",sigma_e2Hat);
   }

   printf("  done.\n");
   fprintf(log,"  done.\n");

//house keeping
   if(CleanKin)
   {
      kin_emp.CleanUpAuto();
      if(AutoFit::fitX)
         kin_emp.CleanUpX();
   }
   //trans.CleanUp();


   //   printf("sigma_g2Hat,sigma_e2Hat,sigma_gXHat,sigma_s2Hat are: %g, %g, %g, %g\n",sigma_g2Hat,sigma_e2Hat,sigma_gXHat,sigma_s2Hat);
}

void FastFit::FastFitPolyGenicModels(Pedigree & ped,double tol,FastTransform & trans,KinshipEmp & kin_emp, FILE * log, bool forX)
{
   if(forX)
   {
      beta.Clear();
      betaHat.Clear();
      beta.Dimension(trans.fixedEffects);
      betaHat.Dimension(trans.fixedEffects);
   }
   trans.Prepare(ped,traitNum,useCovariates,log,true);

   /*
   //TS
   if(forX)
   {
   printf("Y's are:\n");
   for(int i=0;i<trans.persons;i++)
   printf("%g\t",trans.Y[i]);
   printf("\n");
   printf("X's are:\n");
   for(int i=0;i<trans.persons;i++)
   {
   for(int j=0;j<trans.X.cols;j++)
   printf("%g\t",trans.X[i][j]);
   }
   printf("\n");
   }
   */

   if(FastTransform::pedKin) {
      if(forX)
         trans.TransformPedkinSepX(ped,traitNum,useCovariates,log);
      else
         trans.TransformPedkin(ped,traitNum,useCovariates,log);
   }
   else
   {
      if(forX)
         trans.TransformEmpkin(ped,traitNum,useCovariates,kin_emp,kin_emp.allPairsX,log);
      else
         trans.TransformEmpkin(ped,traitNum,useCovariates,kin_emp,kin_emp.allPairs,log);
   }
   OptimizeDelta(0.0001,trans);

   //Save the parameter estiamtes
   GetBeta(deltaHat,trans);
   betaHat = beta;
   if(forX)
   {
      sigma_gXHat = GetSigma_g(deltaHat,betaHat,trans);
      sigma_g2Hat= sigma_gXHat;
      sigma_e2Hat = sigma_gXHat * deltaHat;
   }
   else
   {
      sigma_g2Hat = GetSigma_g(deltaHat,betaHat,trans);
      sigma_e2Hat = sigma_g2Hat * deltaHat;
   }

   if(forX)
      Xheritability = sigma_gXHat/(sigma_gXHat+sigma_e2Hat);
   else
      heritability = sigma_g2Hat/(sigma_g2Hat+sigma_e2Hat);
   //printf("Model fitting is finished.\n");
   fprintf(log,"Fitting linear mixed model completed.\n");

   //Fill in UDX, UD, and UDY
   /*
      if(forX)
      trans.FinalizeTransform(sigma_gXHat, sigma_e2Hat);
      else
      trans.FinalizeTransform(sigma_g2Hat, sigma_e2Hat);
   //house keeping
    */
   if(forX)
   {
      if(CleanKin)
      {
	 kin_emp.CleanUpX();
	 kin_emp.CleanUpAuto();
      }
   }
   /*
   //TS:
   printf("betaHat is:\n");
   for(int i=0;i<betaHat.Length();i++)
   printf("%g ",betaHat[i]);
   printf("\n");
   //TS
    */
   //printf("sigma_g2Hat,sigma_e2Hat,sigma_gXHat,sigma_s2Hat are: %g, %g, %g, %g\n",sigma_g2Hat,sigma_e2Hat,sigma_gXHat,sigma_s2Hat);
}

void FastFit::AutoFitLinearMixModels2(Pedigree & ped,double tol,FastTransform & trans,KinshipEmp & kin_emp,FILE * log,bool forX)
{
   //TS
   if(forX)
   {
      beta.Clear();
      betaHat.Clear();
      beta.Dimension(trans.fixedEffects);
      betaHat.Dimension(trans.fixedEffects);
   }
   AutoFit2 fitEngine;
   fitEngine.useCovariates = useCovariates;
   trans.Prepare(ped,traitNum,fitEngine.useCovariates,log,true);
   fitEngine.FitPolygenicModels(ped,kin_emp,trans,traitNum,true);

   sigma_g2Hat = fitEngine.heritability * fitEngine.variance;
   sigma_e2Hat = fitEngine.sigma_e2;
   sigma_gXHat = fitEngine.heritabilityX * fitEngine.variance;
   betaHat.Copy(fitEngine.betaHat);

   heritability = sigma_g2Hat/(sigma_g2Hat+sigma_e2Hat);
   if(forX)
   {
      heritabilityAuto = sigma_g2Hat/fitEngine.variance;
      Xheritability = sigma_gXHat/fitEngine.variance;
   }

   //printf("Model fitting is finished.\n");
   if(forX)
   {
      printf("    Fitting linear mixed model for chromosome X completed.\n");
      fprintf(log,"    Fitting linear mixed model for chromosome X completed.\n");
   }
   else
   {
      printf("    Fitting linear mixed model completed.\n");
      fprintf(log,"    Fitting linear mixed model completed.\n");
   }
   //SVD var-cov matrix and do the transformations
   trans.TransformEmpkinX(fitEngine.covMatrix);
   //TS
   //printf("sigma_g2Hat,sigma_e2Hat,sigma_gXHat,sigma_s2Hat are: %g, %g, %g, %g\n",sigma_g2Hat,sigma_e2Hat,sigma_gXHat,sigma_s2Hat);

   //Fill in UDX, UD, and UDY
   //trans.FinalizeTransform(sigma_g2Hat, sigma_e2Hat);
}

void FastFit::AutoFitLinearMixModels(Pedigree & ped,double tol,FastTransform & trans,KinshipEmp & kin_emp,FILE * log, bool forX)
{
   if(forX)
   {
      beta.Clear();
      betaHat.Clear();
      beta.Dimension(trans.fixedEffects);
      betaHat.Dimension(trans.fixedEffects);
   }
   AutoFit fitEngine;
   fitEngine.useCovariates = useCovariates;
   trans.Prepare(ped,traitNum,fitEngine.useCovariates,log,true);
   fitEngine.FitPolygenicModels(ped,traitNum,true);

   sigma_g2Hat = fitEngine.heritability * fitEngine.variance;
   sigma_e2Hat = fitEngine.sigma_e2 * fitEngine.variance;
   sigma_gXHat = fitEngine.sigma_gX * fitEngine.variance;
   sigma_s2Hat = fitEngine.sharedEnv * fitEngine.variance;
   betaHat.Copy(fitEngine.betaHat);

   //printf("Model fitting is finished.\n");
   if(forX)
   {
      printf("    Fitting linear mixed model for chromosome X completed.\n");
      fprintf(log,"    Fitting linear mixed model for chromosome X completed.\n");
   }
   else
   {
      printf("    Fitting linear mixed model completed.\n");
      fprintf(log,"    Fitting linear mixed model completed.\n");
   }
   //SVD var-cov matrix and do the transformations
   trans.TransformPedkinX(fitEngine, ped);
   heritability = sigma_g2Hat/fitEngine.variance;
   if(forX)
   {
      heritabilityAuto = sigma_g2Hat/fitEngine.variance;
      Xheritability = sigma_gXHat/fitEngine.variance;
   }

   //Fill in UDX, UD, and UDY
   //trans.FinalizeTransform(sigma_g2Hat, sigma_e2Hat);
   //printf("sigma_g2Hat,sigma_e2Hat,sigma_gXHat,sigma_s2Hat are: %g, %g, %g, %g\n",sigma_g2Hat,sigma_e2Hat,sigma_gXHat,sigma_s2Hat);
}

void FastFit::FitSimpleLinearModels(Pedigree & ped,double tol,FastTransform & trans, KinshipEmp & kin_emp,FILE * log)
{
   trans.Prepare(ped,traitNum,useCovariates,log,true);

   /* debug
   printf("\nX size = %d,%d\n", trans.X.rows, trans.X.cols);
   printf("X:");
   for(int i=0; i<trans.X.rows; i++) {
      for(int j=0; j<trans.X.cols; j++)
         printf(",%g",trans.X[i][j]);
      printf("\n");
   }*/

   Matrix inv;
   Matrix transX;
   transX.Transpose(trans.X);

   Vector Y;
   for (int f = 0; f < ped.familyCount; f++)
   {
      for(int j=0;j<trans.pPheno[f].Length();j++)
         Y.Push(ped[trans.pPheno[f][j]].traits[traitNum]);
   }
   bool is_possible_binary = isBinaryY(Y);
   if (is_possible_binary)
      printf("Warning: only 2 values found in phenotype. Should we try --binary instead?");

   Vector tmp;
   tmp.Dimension(trans.fixedEffects);
   for(int j=0;j<trans.fixedEffects;j++)
      tmp[j] = Y.InnerProduct(transX[j]);

   inv.Product(transX,trans.X);
   SVD svd;
   svd.InvertInPlace(inv);

   for(int i=0;i<trans.fixedEffects;i++)
      betaHat[i] = inv[i].InnerProduct(tmp);

   residuals_unrelated.Dimension(trans.persons);
   for(int i=0;i<trans.persons;i++)
      residuals_unrelated[i] = Y[i] - betaHat.InnerProduct(trans.X[i]);

   //         for(int i=0;i<trans.persons;i++)
   //printf("%g\n",residuals_unrelated[i]);
   sigma_g2Hat=0.0;
   sigma2Hat = residuals_unrelated.SumSquares()/trans.persons;
   //   for(int i=0;i<trans.persons;i++)
   //     residuals_unrelated[i] /= sigma2Hat;
}

void FastFit::FitSimpleLogisticModels(Pedigree & ped, FastTransform & trans,FILE * log)
{
   trans.Prepare(ped,traitNum,useCovariates,log,true);

   Vector Y;
   for (int f = 0; f < ped.familyCount; f++) {
      for(int j=0;j<trans.pPheno[f].Length();j++) {
//printf("%g\n",ped[trans.pPheno[f][j]].traits[traitNum]);
         Y.Push(ped[trans.pPheno[f][j]].traits[traitNum]);
      }
   }

   // fit model
   int Yplus1, Yminus1;
   convertBinaryTrait(Y, Yplus1, Yminus1);
   GetLogitCoeff(betaHat, trans.X, Y);

   residuals_unrelated.Dimension(trans.persons);
   for(int i=0;i<trans.persons;i++)
      residuals_unrelated[i] = Y[i] - sigmoid(betaHat.InnerProduct(trans.X[i]));

   sigma_g2Hat=0.0;
   sigma2Hat = residuals_unrelated.SumSquares()/trans.persons;
}


int FastFit::GetTraitID( Pedigree & ped, const char * name,FILE * log)
{
   int idx = ped.traitLookup.Integer(name);
   if (idx != -1)
      return idx;
   else {
      fprintf(log,"Error! Trait not found. Please check .dat file for the correct trait name.\n");
      error("Error! Trait not found. Please check .dat file for the correct trait name.\n");
      return -1;
   }
}

double FastFit::Brent(double tol,FastTransform & trans)
{
   double temp;

   if (a > c) {
      SHIFT(temp, a, c, temp);
      SHIFT(temp, fa, fc, temp);
   }

   min = b; fmin = fb;
   double w = b, v = b;
   double fw = fb, fv = fb;

   double delta = 0.0;         // distance moved in step before last
   double u, fu, d = 0.0;      // Initializing d is not necesary, but avoids warnings
   for (int iter = 1; iter <= ITMAX; iter++) {
      double middle = 0.5 * (a + c);
      double tol1 = tol * fabs(min) + ZEPS;
      double tol2 = 2.0 * tol1;

      if (fabs(min - middle) <= (tol2 - 0.5 * (c - a)))
         return fmin;

      if (fabs(delta) > tol1) {// Try a parabolic fit
         double r = (min - w) * (fmin - fv);
         double q = (min - v) * (fmin - fw);
         double p = (min - v) * q - (min - w) * r;

         q = 2.0 * (q - r);
         if (q > 0.0) p = -p;
            q = fabs(q);

         temp = delta;
         delta = d;

         if (fabs(p) >= fabs(0.5 * q * temp) || p <= q * (a - min) || p >= q * (c - min)) {
            // parabolic doesn't look like a good idea
            delta = min >= middle ? a - min : c - min;
            d = CGOLD * delta;
         }
         else {
         // parabolic fit is the way to go
            d = p / q;
            u = min + d;
            if (u - a < tol2 || c - u < tol2)
               d = sign(tol1, middle - min);
         }
      }
      else { // Golden ratio for first step
         delta = min >= middle ? a - min : c - min;
         d = CGOLD * delta;
      }

      // Don't take steps smaller than tol1
      u = fabs(d) >= tol1 ? min + d : min + sign(tol1, d);
      fu = -1.0*Evaluate(u,trans,false);

      if (fu <= fmin) {
         if (u >= min)
            a = min;
         else
            c = min;
         SHIFT(v, w, min, u);
         SHIFT(fv, fw, fmin, fu);
      }
      else {
         if (u < min)
            a = u;
         else
            c = u;
         if (fu <= fw || w == min) {
            SHFT3(v, w, u);
            SHFT3(fv, fw, fu);
         }
         else if (fu <= fv || v == min || v == w) {
            v = u;
            fv = fu;
         }
      }
   }
   numerror("ScalarMinimizer::Brent got stuck");
   return fmin;
}

