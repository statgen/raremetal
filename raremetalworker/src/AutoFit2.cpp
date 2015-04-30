#include "AutoFit2.h"
#include "FastFit.h"
#include "MathNormal.h"
#include "KinshipEmp.h"
#include "TransformResiduals.h"

#include <math.h>


//bool AutoFit2::useProbands = false;
bool AutoFit2::useCovariates = false;

//double mean =0.0;
//double heritability =0.0;
//double heritabilityX = 0.0;
//double variance = 0.0;
//double sigma_e2=0.0;

AutoFit2::AutoFit2() { }

void AutoFit2::FitPolygenicModels(Pedigree & ped, KinshipEmp & kin_emp,FastTransform & trans, int traitNum,bool quiet)
{
   Matrix allPairs,allPairsX;
   allPairs.Dimension(trans.persons,trans.persons);
   if(AutoFit::fitX)
      allPairsX.Dimension(trans.persons,trans.persons);

   //Subset empirical kinship matrix
   if(FastTransform::readInEmp!="")
   {
      printf("\n    Matching IDs in kinship matrix ... ");
      fflush(stdout);
      IntArray selectedSample;
      //Checking if all individuals are included in read-in kinship.
      for (int f = 0; f < ped.familyCount; f++)
      {
	 for(int j=0;j<trans.pPheno[f].Length();j++)
	 {
	    int p = kin_emp.IDFromEmp.Integer(ped[trans.pPheno[f][j]].pid);
	    if(p==-1)
	    {
	       error("Fatal Error: person %s in ped file is not included in read-in kinship. The program will stop here. Please fix the problem and re-run.\n",ped[trans.pPheno[f][j]].pid.c_str());
	       //fprintf(log,"Fatal Error: person %s in ped file is not included in read-in kinship. The program will stop here. Please fix the problem and re-run.\n",ped[trans.pPheno[f][j]].pid.c_str());
	    }
	    else
	    {
	       selectedSample.Push(p);
	    }
	 }
	 printf("done.\n");
      }

      //subsetting
      for(int r=0;r<selectedSample.Length();r++)
      {
	 for(int c=0;c<selectedSample.Length();c++)
	 {
	    allPairs[r][c] = kin_emp.allPairs[selectedSample[r]][selectedSample[c]];
	 }
      }
      if(AutoFit::fitX)
      {
	 for(int r=0;r<selectedSample.Length();r++)
	 {
	    for(int c=0;c<selectedSample.Length();c++)
	    {
	       allPairsX[r][c] = kin_emp.allPairsX[selectedSample[r]][selectedSample[c]];
	    }
	 }
      }
   }
   else
   {
      trans.SubSet(kin_emp.allPairs,ped,allPairs);
      if(AutoFit::fitX)
	 trans.SubSet(kin_emp.allPairsX,ped,allPairsX);
   }

   if (!quiet)
   {
      // Print header line
      printf("FITTED MODELS%s\n",
	    useCovariates ? " (for covariate adjusted residuals)" : "");
      printf("======================================================\n"
	    "%15s %12s %12s %12s\n",
	    "Trait", "Mean", "Variance", "Heritability");
   }

   // Loop through traits in the pedigree
   //for (int t = 0; t < ped.traitCount; t++)
   //   {
   // First list phenotyped individuals 
   // If we are using covariates, we only consider individuals
   // for which all covariates have been recorded

   int vc_count;
   if(AutoFit::fitX)
      vc_count = 3;
   else
      vc_count = 2;

   // The normal set class is our workhorse
   NormalSet mvn;
   mvn.Dimension(1, vc_count);

   // Index to the normal equations within NormalSet
   int index = 0;

   // Build model with polygenes and environment only
   int count = trans.persons;

   // Phenotypes
   mvn[index].scores.Dimension(count);
   int idx=0;
   for (int f = 0; f < ped.familyCount; f++)
   {
      for(int j=0;j<trans.pPheno[f].Length();j++)
      {
	 /*
	    if(PreMeta::genoFromVCF || PreMeta::dosage)
	    {
	    int p = trans.sampleVCFIDHash.Integer(ped[pPheno[f][j]].pid);
	    if(p==-1)
	    {
	    error("Fatal Error: person %s in ped file is not included in the VCF. The program will abort.\n",ped[trans.pPheno[f][j]].pid.c_str());
	 //fprintf(log,"Fatal Error: person %s in ped file is not included in read-in kinship. The program will stop here. Please fix the problem and re-run.\n",ped[trans.pPheno[f][j]].pid.c_str());
	 }
	 }
	 else if(PreMeta::genoFromPed)
	 {
	 int p = samplePEDIDHash.Integer(ped[i].pid);
	 if(p==-1)
	 {
	 error("Fatal Error: person %s in ped file is not included in the PED. The program will abort.\n",ped[trans.pPheno[f][j]].pid.c_str());
	 //fprintf(log,"Fatal Error: person %s in ped file is not included in read-in kinship. The program will stop here. Please fix the problem and re-run.\n",ped[trans.pPheno[f][j]].pid.c_str());
	 }
	 }
	  */
	 mvn[index].scores[idx] = ped[trans.pPheno[f][j]].traits[traitNum];
	 idx++;
      }
   }

   // Total number of fixed effects
   int fixedEffects = 1 + (useCovariates ? ped.covariateCount : 0);

   // Setup matrix of fixed effects for each family
   mvn[index].linearModel.Dimension(count, fixedEffects);
   for (int i = 0; i < count; i++)
   {
      // Constant for regressing grand mean
      mvn[index].linearModel[i][0] = 1.0;
   }

   // User specified covariates
   if (useCovariates)
   {
      int idx=0;
      for (int f = 0; f < ped.familyCount; f++)
      {
	 for(int j=0;j<trans.pPheno[f].Length();j++)
	 {
	    /*
	       if(PreMeta::genoFromVCF || PreMeta::dosage)
	       {
	       int p = trans.sampleVCFIDHash.Integer(ped[pPheno[f][j]].pid);
	       if(p==-1)
	       {
	       error("Fatal Error: person %s in ped file is not included in the VCF. The program will abort.\n",ped[trans.pPheno[f][j]].pid.c_str());
	    //fprintf(log,"Fatal Error: person %s in ped file is not included in read-in kinship. The program will stop here. Please fix the problem and re-run.\n",ped[trans.pPheno[f][j]].pid.c_str());
	    }
	    }
	    else if(PreMeta::genoFromPed)
	    {
	    int p = samplePEDIDHash.Integer(ped[i].pid);
	    if(p==-1)
	    {
	    error("Fatal Error: person %s in ped file is not included in the PED. The program will abort.\n",ped[trans.pPheno[f][j]].pid.c_str());
	    //fprintf(log,"Fatal Error: person %s in ped file is not included in read-in kinship. The program will stop here. Please fix the problem and re-run.\n",ped[trans.pPheno[f][j]].pid.c_str());
	    }
	    }
	     */
	    for (int c = 1; c <= ped.covariateCount; c++)
	       mvn[index].linearModel[idx][c] = ped[trans.pPheno[f][j]].covariates[c-1];
	    idx++;
	 }
      }
   }

   // Setup non-shared variances
   mvn[index].varComponents[0].Dimension(count, count);
   mvn[index].varComponents[0].Identity();

   // Setup polygenic variances
   mvn[index].varComponents[1].Dimension(count, count);

   for (int i = 0; i < count; i++)
      for (int j = i; j < count; j++)
	 mvn[index].varComponents[1][i][j] =
	    mvn[index].varComponents[1][j][i] =
	    allPairs[i][j];

   if(AutoFit::fitX)
   {
      // Setup polygenic variances for X chromosome
      mvn[index].varComponents[2].Dimension(count, count);
      for (int i = 0; i < count; i++)
	 for (int j = i; j < count; j++)
	    mvn[index].varComponents[2][i][j] =
	       mvn[index].varComponents[2][j][i] =
	       allPairsX[i][j];
   }

   // Fit polygenic model
   mvn.Solve();

   mean = useCovariates ? 0.0 : mvn.means[0];
   variance = mvn.variances.Sum();
   heritability = mvn.variances[1] / variance;
   sigma_e2=mvn.variances[0];
   betaHat.Copy(mvn[0].means);

   if(AutoFit::fitX)
      heritabilityX = mvn.variances[2]/variance;

   if (!quiet)
      printf("%15s %12.5f %12.5f %12.5f\n",
	    (const char *) ped.traitNames[traitNum], mean, variance, heritability);

   /*
      if (!useCovariates) continue;

   // Evaluate model to update residuals for each family
   mvn.Evaluate();

   index = 0;

   // First wipe out all the phenotypes, to ensure that we
   // don't get residuals and covariates mixed up ...
   for (int i = 0; i < ped.count; i++)
   ped[i].traits[traitNum] = _NAN_;

   // Replace phenotypes with residuals
   idx=0;
   for (int f = 0; f < ped.familyCount; f++)
   {
   for(int j=0;j<trans.pPheno[f].Length();j++)
   {
   ped[trans.pPheno[f][j]].traits[traitNum] = mvn[index].residuals[idx];
   idx++;
   }
   }
   ped[pheno[i]].traits[traitNum] = mvn[index].residuals[i];
   //}
    */

   //generate covariance matrix and save
   covMatrix.Dimension(trans.persons,trans.persons);
   if(AutoFit::fitX)
   {
      for(int i=0;i<trans.persons;i++)
      {
	 covMatrix[i][i]=allPairs[i][i]*mvn.variances[1]+allPairsX[i][i]*mvn.variances[2] + mvn.variances[0];
	 for(int j=i+1;j<trans.persons;j++)
	    covMatrix[j][i]=covMatrix[i][j] = allPairs[i][j]*mvn.variances[1]+allPairsX[i][j]*mvn.variances[2];
      }
   }
   else
   {
      for(int i=0;i<trans.persons;i++)
      {
	 covMatrix[i][i]=allPairs[i][i]*mvn.variances[1]+mvn.variances[0];
	 for(int j=i+1;j<trans.persons;j++)
	    covMatrix[j][i]=covMatrix[i][j] = allPairs[i][j]*mvn.variances[1];
      }
   }
}
