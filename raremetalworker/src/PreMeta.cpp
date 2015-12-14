// PreMeta.cpp 
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

#include "Constant.h"
#include "math.h"
#include <string>
#include <cstring>
#include <iostream>
#include "PreMeta.h"
#include "snpHWE.h"
#include "WritePDF.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "InputFile.h"
#include "GenomeSequence.h"
//#include "CheckRef.h"

bool PreMeta::recessive = false;
bool PreMeta::dominant = false;
bool PreMeta::additive = true;
bool PreMeta::checkRef = false;
String PreMeta::vcfInput = ""; //input annotated vcf file name here
String PreMeta::dosageFlag = "DS";
int PreMeta::Xstart = 2699520;
int PreMeta::Xend = 154931044;
bool PreMeta::dosage = false;
bool PreMeta::genoFromPed = false;
bool PreMeta::genoFromVCF = false;
bool PreMeta::FounderFreq = false;
bool PreMeta::zipOutput=false;
String PreMeta::outputFile = "";
String PreMeta::xLabel = "X";
int PreMeta::window = 1000000;
int PreMeta::maleLabel = 1;
int PreMeta::femaleLabel =2;
bool PreMeta::correctGC =false;
bool PreMeta::calculateOR=false;
bool PreMeta::Simplify=false; // print %.2e for cov
String PreMeta::Region = "";
String PreMeta::varListName = "";

PreMeta::PreMeta(){}

PreMeta::~PreMeta()
{}

void PreMeta::CalculateProjectionMatrix(FastTransform & trans,FastFit & engine,Vector & sigma2)
{
   projectionMat.Dimension(trans.persons,trans.persons);
   Matrix tmp,tmp2,trans_UX;
   if(FastFit::unrelated)
   {
      tmp.Product(trans.X,engine.inv);
      trans_UX.Transpose(trans.X);
      tmp2.Product(tmp,trans_UX);

      for(int i=0;i<projectionMat.rows;i++)
      {
	 projectionMat[i][i] =1.0-tmp2[i][i];
	 for(int j=i+1;j<projectionMat.cols;j++)
	 {
	    projectionMat[i][j] = projectionMat[j][i] = -tmp2[i][j];
	 }
      }

      for(int i=0;i<projectionMat.rows;i++)
	 for(int j=0;j<projectionMat.cols;j++)
	    projectionMat[i][j] /=engine.sigma2Hat;
   }
   else
   {
      tmp.Product(trans.UX,engine.inv);
      trans_UX.Transpose(trans.UX);
      tmp2.Product(tmp,trans_UX);

      for(int i=0;i<projectionMat.rows;i++)
      {
	 projectionMat[i][i] =1.0-tmp2[i][i];
	 for(int j=i+1;j<projectionMat.cols;j++)
	 {
	    projectionMat[i][j] = projectionMat[j][i] = -tmp2[i][j];
	 }
      }

      for(int i=0;i<projectionMat.rows;i++)
	 for(int j=0;j<projectionMat.cols;j++)
	    projectionMat[i][j] *= 1.0/sigma2[i];
   }
}

double PreMeta::CalculateCov(FastFit & engine,FastTransform & trans, Pedigree & ped, Matrix & genotypeAll,Vector & sigma2,int m)
{
   double cov = 0.0;
/*
   if(!FastFit::useCovariates)
   {
*/
      if(FastFit::unrelated)
      {    
	 cov = genotypeAll[0].InnerProduct(genotypeAll[m]);
	 cov /= engine.sigma2Hat;
      }    
      else 
      {    
	 double scale;
	 for(int p=0;p<trans.persons;p++)
	 {     
	    scale = 1.0/sigma2[p];
	    cov += genotypeAll[0][p] * genotypeAll[m][p]*scale;
	 }     
      }   
/*   }
   else
   {
      Vector tmp,tmp2;
      tmp.Dimension(trans.persons,0.0);
      for(int i=0;i<trans.persons;i++)
	 for(int j=0;j<trans.persons;j++)
	    tmp[i] += genotypeAll[0][i]*projectionMat[j][i];
      for(int i=0;i<trans.persons;i++)
	 cov += genotypeAll[m][i]*tmp[i];
   }
 */
  //printf("cov is: %g\n",cov); 
   cov /= trans.persons;
   return cov;
}

void PreMeta::CalculateAssocStats(double & effSize,double &pvalue,double &numerator,double &denominator,double &chisq,Eigen::VectorXf & transGeno, FastFit & engine,FastTransform & trans,Vector & sigma2)
{
   numerator = 0.0;
   denominator = 0.0;

   for(int i=0;i<trans.persons;i++)
   {
      double tmp = transGeno[i]/sigma2[i];
      numerator += residuals[i]*tmp;
      denominator += transGeno[i]*tmp;
   }
   if(denominator==0.0)
   {
      pvalue=_NAN_;
      effSize = _NAN_;
      return;
   }
   /*
      printf("numerator and denominator are: %g,%g.\n",numerator,denominator);
      printf("residuals are:\n");
      for(int i=0;i<trans.persons;i++)
      printf("%g\n",trans.Y[i]-trans.X[i].InnerProduct(engine.betaHat));

      printf("transGeno is:\n");
      for(int i=0;i<trans.persons;i++)
      printf("%g\n",transGeno[i]);
      printf("transUY is:\n");
      for(int i=0;i<trans.persons;i++)
      printf("%g\n",trans.UY[i]);
      printf("sigma2 is:\n");
      for(int i=0;i<trans.persons;i++)
      printf("%g\n",sigma2[i]);
      printf("\n");
    */
   chisq = numerator*numerator/denominator;
   pvalue = pchisq(chisq, 1,0,0);
   effSize = numerator/denominator;
}

void PreMeta::RelatedAssoc(IFILE SCOREoutput, IFILE SCOREoutput_rec,IFILE SCOREoutput_dom, Pedigree & ped,FastFit & engine,FastTransform & trans,Vector & sigma2)
{
   //get transGeno
   if(trans.pedKin)
   {
      for(int i=0;i<trans.persons;i++)
      {
	 transGeno[i] = 0.0;
	 if(recessive)
	    transGeno_rec[i] = 0.0;
	 if(dominant)
	    transGeno_dom[i] = 0.0;
      }
      int index=0;
      for (int f = 0; f < ped.familyCount; f++)
      {
	 if (trans.pPheno[f].Length() > 0)
	 {
	    int count = trans.pPheno[f].Length();
	    int block = index;
	    for(int i=0;i<count;i++)
	    {
	       for(int c=0;c<count;c++)
		  transGeno[index] += trans.transU(block+i,block+c)*genotype[block+c];
	       if(recessive)
	       {
		  for(int c=0;c<count;c++)
		     transGeno_rec[index] += trans.transU(block+i,block+c)*genotype_rec[block+c];
	       }
	       if(dominant)
	       {
		  for(int c=0;c<count;c++)
		     transGeno_dom[index] += trans.transU(block+i,block+c)*genotype_dom[block+c];
	       }
	       index++;
	    }
	 }
      }
   }
   else
   {
      transGeno = trans.transU*genotype;
      if(recessive)
	 transGeno_rec = trans.transU*genotype_rec;
      if(dominant)
	 transGeno_dom = trans.transU*genotype_dom;
   }
   //calculate score statistics and output
   double recessive_effSize, recessive_pvalue = _NAN_, dominant_effSize, dominant_pvalue=_NAN_, chisq, pvalue = _NAN_,effSize,numerator,denominator;
   CalculateAssocStats(effSize,pvalue,numerator,denominator,chisq,transGeno,engine,trans,sigma2);
   if(pvalue!=_NAN_)
   {
      if(FounderFreq)
      {
	 maf_GC.Push(founderMaf);
      }
      else
      {
	 maf_GC.Push(markerMaf);
      }
      pvalueAll.Push(pvalue);
      chr_plot.Push(chr);
      pos_plot.Push(pos);
      chisq_before_GCcorrect.Push(chisq);
   }

   if(pvalue==_NAN_)
   {
      if(hwe_pvalue != _NAN_ || hwe_pvalue != 6.6666e-66)
      {
	 ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%g\t%g\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2,numerator,sqrt(denominator));
      }
      else
      {
	 ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d\t%g\t%g\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2,numerator,sqrt(denominator));
      }
   }
   else
   {
      if(hwe_pvalue!=_NAN_ || hwe_pvalue != 6.66666e-66)
	 ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%g\t%g\t%g\t%g",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2,numerator,sqrt(denominator),effSize,pvalue);
      else
	 ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d\t%g\t%g\t%g\t%g",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2,numerator,sqrt(denominator),effSize,pvalue);
   }
   ifprintf(SCOREoutput,"\n");

   if(recessive)
   {
      if(hwe_pvalue!=_NAN_ || hwe_pvalue != 6.66666e-66)
	 ifprintf(SCOREoutput_rec,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2);
      else
	 ifprintf(SCOREoutput_rec,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2);
   }
   if(dominant)
   {
      if(hwe_pvalue!=_NAN_  || hwe_pvalue != 6.66666e-66)
	 ifprintf(SCOREoutput_dom,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2);
      else
	 ifprintf(SCOREoutput_dom,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2);
   }
   //recessive model
   if(recessive)
   {
      if(hom2==0 || hom1+het==0)
	 ifprintf(SCOREoutput_rec,"\tNA\tNA\tNA\tNA");
      else
      {
	 CalculateAssocStats(recessive_effSize,recessive_pvalue,numerator,denominator,chisq,transGeno_rec,engine,trans,sigma2);
	 if(recessive_pvalue==_NAN_)
	    ifprintf(SCOREoutput_rec,"\t%g\t%g\tNA\tNA",numerator,sqrt(denominator));
	 else
	 {
	    ifprintf(SCOREoutput_rec,"\t%g\t%g\t%g\t%g",numerator,sqrt(denominator),recessive_effSize,recessive_pvalue);
	 }
      }
      //if(recessive_pvalue != _NAN_  || hwe_pvalue != 6.66666e-66)
      if(recessive_pvalue != _NAN_)
      {
	 chisq_before_GCcorrect_rec.Push(chisq);
	 pvalueAll_rec.Push(recessive_pvalue);
	 if(FounderFreq)
	    maf_GC_rec.Push(founderMaf);
	 else
	    maf_GC_rec.Push(markerMaf);
	 chr_plot_rec.Push(chr);
	 pos_plot_rec.Push(pos);
      }
      ifprintf(SCOREoutput_rec,"\n");
   }
   //dominant model
   if(dominant)
   {
      if(hom1==0 || hom2+het==0)
	 ifprintf(SCOREoutput_dom,"\tNA\tNA\tNA\tNA");
      else
      {
	 CalculateAssocStats(dominant_effSize,dominant_pvalue,numerator,denominator,chisq,transGeno_dom,engine,trans,sigma2);
	 if(dominant_pvalue==_NAN_)
	    ifprintf(SCOREoutput_dom,"\t%g\t%g\tNA\tNA",numerator,sqrt(denominator));
	 else
	    ifprintf(SCOREoutput_dom,"\t%g\t%g\t%g\t%g",numerator,sqrt(denominator),dominant_effSize,dominant_pvalue);
      }
      if(dominant_pvalue != _NAN_)
      {
	 chisq_before_GCcorrect_dom.Push(chisq);
	 pvalueAll_dom.Push(dominant_pvalue);
	 if(FounderFreq)
	    maf_GC_dom.Push(founderMaf);
	 else
	    maf_GC_dom.Push(markerMaf);
	 chr_plot_dom.Push(chr);
	 pos_plot_dom.Push(pos);
      }
      ifprintf(SCOREoutput_dom,"\n");
   }
}

void PreMeta::UnrelatedAssoc(IFILE SCOREoutput, IFILE SCOREoutput_rec,IFILE SCOREoutput_dom, Pedigree & ped,FastFit & engine,FastTransform & trans,Vector & sigma2)
{
   double effSize,pvalue=_NAN_,chisq,numerator,denominator,recessive_effSize,recessive_pvalue = _NAN_,dominant_effSize,dominant_pvalue = _NAN_;
   //Additive model
   CalculateUnrelatedAssocStats(effSize,pvalue,chisq,numerator,denominator,genotype,engine,trans,sigma2);

   if(pvalue!=_NAN_)
   {
      if(FounderFreq)
      {
	 maf_GC.Push(founderMaf);
      }
      else
      {
	 maf_GC.Push(markerMaf);
      }
      pvalueAll.Push(pvalue);
      chr_plot.Push(chr);
      pos_plot.Push(pos);
      //calculate standard error here
      chisq_before_GCcorrect.Push(chisq);
   }

   if(pvalue==_NAN_)
   {
      if(hwe_pvalue!=_NAN_  || hwe_pvalue != 6.66666e-66)
	 ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%g\t%g\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2,numerator,sqrt(denominator));
      else
	 ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d\t%g\t%g\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2,numerator,sqrt(denominator));
   }
   else
   {
      if(hwe_pvalue!=_NAN_  || hwe_pvalue != 6.66666e-66)
	 ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%g\t%g\t%g\t%g",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2,numerator,sqrt(denominator),effSize,pvalue);
      else
	 ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d\t%g\t%g\t%g\t%g",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2,numerator,sqrt(denominator),effSize,pvalue);
   }
   ifprintf(SCOREoutput,"\n");
   //recessive model
   if(recessive)
   {
      if(hwe_pvalue!=_NAN_  || hwe_pvalue != 6.66666e-66)
	 ifprintf(SCOREoutput_rec,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2);
      else
	 ifprintf(SCOREoutput_rec,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2);

      if(hom2==0 || hom1+het==0)
	 ifprintf(SCOREoutput_rec,"\tNA\tNA\tNA\tNA");
      else
      {
	 CalculateUnrelatedAssocStats(recessive_effSize,recessive_pvalue,chisq,numerator,denominator,genotype_rec,engine,trans,sigma2);
	 if(recessive_pvalue==_NAN_)
	    ifprintf(SCOREoutput_rec,"\t%g\t%g\tNA\tNA",numerator,sqrt(denominator));
	 else
	 {
	    ifprintf(SCOREoutput_rec,"\t%g\t%g\t%g\t%g",numerator,sqrt(denominator),recessive_effSize,recessive_pvalue);
	 }
      }
      if(recessive_pvalue != _NAN_)
      {
	 chisq_before_GCcorrect_rec.Push(chisq);
	 pvalueAll_rec.Push(recessive_pvalue);
	 if(FounderFreq)
	    maf_GC_rec.Push(founderMaf);
	 else
	    maf_GC_rec.Push(markerMaf);
	 chr_plot_rec.Push(chr);
	 pos_plot_rec.Push(pos);
      }
      ifprintf(SCOREoutput_rec,"\n");
   }
   //dominant model
   if(dominant)
   {
      if(hwe_pvalue!=_NAN_  || hwe_pvalue != 6.66666e-66)
	 ifprintf(SCOREoutput_dom,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2);
      else
	 ifprintf(SCOREoutput_dom,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2);

      if(hom1==0 || hom2+het==0)
	 ifprintf(SCOREoutput_dom,"\tNA\tNA\tNA\tNA");
      else
      {
	 CalculateUnrelatedAssocStats(dominant_effSize,dominant_pvalue,chisq,numerator,denominator,genotype_dom,engine,trans,sigma2);
	 if(dominant_pvalue==_NAN_)
	    ifprintf(SCOREoutput_dom,"\t%g\t%g\tNA\tNA",numerator,sqrt(denominator));
	 else
	 {
	    ifprintf(SCOREoutput_dom,"\t%g\t%g\t%g\t%g",numerator,sqrt(denominator),dominant_effSize,dominant_pvalue);
	 }
      }
      if(dominant_pvalue != _NAN_)
      {
	 chisq_before_GCcorrect_dom.Push(chisq);
	 pvalueAll_dom.Push(dominant_pvalue);
	 if(FounderFreq)
	    maf_GC_dom.Push(founderMaf);
	 else        
	    maf_GC_dom.Push(markerMaf); 
	 chr_plot_dom.Push(chr);
	 pos_plot_dom.Push(pos);
      }
      ifprintf(SCOREoutput_dom,"\n");
   }
}

void PreMeta::CalculateUnrelatedAssocStats(double & effSize,double & pvalue, double & chisq, double & numerator,double & denominator, Eigen::VectorXf & genotype,FastFit & engine,FastTransform & trans,Vector & sigma2)
{
   numerator = 0.0;
   denominator = 0.0;

   for(int i=0;i<trans.persons;i++)
   {
      //troubleshooting
      //printf("%g\n",engine.residuals_unrelated[i]);
      numerator += engine.residuals_unrelated[i]*genotype[i];
      denominator += genotype[i]*genotype[i];
   }
   if(denominator==0)
   {
      effSize = _NAN_;
      pvalue = _NAN_;
      return;
   }
   effSize = numerator/denominator;
   //denominator *= engine.sigma2Hat;
   double inv_var = 1.0/engine.sigma2Hat;
   numerator *= inv_var; 
   denominator *= inv_var;

   chisq = numerator*numerator/denominator;
   pvalue = pchisq(chisq, 1,0,0);
}

void PreMeta::Run(IFILE SCOREoutput, IFILE SCOREoutput_rec,IFILE SCOREoutput_dom, IFILE SCOREcov, IFILE SCOREcov_rec,IFILE SCOREcov_dom,Pedigree & ped, FastTransform & trans,FastFit & engine,GroupFromAnnotation & group, FILE * log,SanityCheck & checkData,KinshipEmp & kin_emp)
{
	if ( PreMeta::varListName != "" ) {
		printf("\nReading variant list...\n");
		setVarList();
		printf("  done.\n\n ");
	}

   warnings = 0;
   bool fitXStatus = false;
   //setup the reference genome
   //CheckRef checkref;
   //if(checkRef)
   //   checkref.Setup();

   double averageAF=0.0;
   int marker_count=0;
   malehwewarning=0;
   /*
      if(genoFromVCF)
      {
   //This is to make sure VCF file is sorted
   VCFSanityCheck();
   }
    */

   PrintHeader(SCOREoutput,engine,trans,ped);
   ifprintf(SCOREcov,"##ProgramName=RareMetalWorker\n");
   ifprintf(SCOREcov,"##Version=%s\n",VERSION);
   ifprintf(SCOREcov,"#CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES\n");

   if(recessive)
   {
      PrintHeader(SCOREoutput_rec,engine,trans,ped);
      ifprintf(SCOREcov_rec,"##ProgramName=RareMetalWorker\n");
      ifprintf(SCOREcov_rec,"##Version=%s\n",VERSION);
      ifprintf(SCOREcov_rec,"#CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES\n");
   }
   if(dominant)
   {
      PrintHeader(SCOREoutput_dom,engine,trans,ped);
      ifprintf(SCOREcov_dom,"##ProgramName=RareMetalWorker\n");
      ifprintf(SCOREcov_dom,"##Version=%s\n",VERSION);
      ifprintf(SCOREcov_dom,"#CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES\n");
   }

   Matrix genotypeAll,genotypeAll_rec,genotypeAll_dom; //this matrix has the all marker genotype within a window
   genotypeAll.Dimension(0,0);
   if(recessive)
      genotypeAll_rec.Dimension(0,0);
   if(dominant)
      genotypeAll_dom.Dimension(0,0);

   IntArray markerInLD;

   Vector sigma2;
   sigma2.Dimension(trans.persons);
   transGeno.resize(trans.persons);
   if(recessive)
      transGeno_rec.resize(trans.persons);
   if(dominant)
      transGeno_dom.resize(trans.persons);

   double sigma_g2 = engine.sigma_g2Hat;
   double sigma_e2 = engine.sigma_e2Hat;

   if(!FastFit::unrelated)
   {
      for(int i=0;i<trans.persons;i++)
      {
	 sigma2[i] = sigma_g2*trans.D[i]+sigma_e2;
      }
      residuals.Dimension(trans.persons);
      for(int i=0;i<trans.persons;i++)
	 residuals[i]=trans.UY[i]-trans.UX[i].InnerProduct(engine.betaHat);
   }
   else
   {
      for(int i=0;i<trans.persons;i++)
	 sigma2[i] = engine.sigma2Hat;
   }

   if(FastFit::useCovariates)
      CalculateProjectionMatrix(trans,engine,sigma2);

   if(genoFromPed)
   {
      //TODO: this assumes all markers are sorted by chromosome and position
      StringArray name;
      String initial_chr;
      int initial_pos = 0;
      int startingMarker=_NAN_;
      for(int i=0;i<ped.markerCount;i++)
      {
	 name.AddTokens(ped.markerNames[i],":");
	 /*
	    if(name.Length()<=1)
	    error("ERROR! Marker names in DAT file has to be in CHR:POS format.\n");
	  */
	 if(name.Length()<=1)
	 {
	    printf("WARNING: marker name %s in DAT file is not in CHR:POS format and will be skipped from analysis.\n",ped.markerNames[i].c_str());
	    fprintf(log,"WARNING: marker name %s in DAT file is not in CHR:POS format and will be skipped from analysis.\n",ped.markerNames[i].c_str());
	    continue;
	 }
	 initial_chr = name[0];
	 initial_pos = atoi(name[1].c_str());
	 startingMarker=i;
	 break;
      }

      for(int m=startingMarker; m<ped.markerCount; m++)
      {
	 int status = GetGenotypeVectorPED(trans,ped,m,log);

	 if(status==2)
	    continue;

	 if(chr==PreMeta::xLabel && (AutoFit::fitX || FastFit::separateX) && !fitXStatus && !FastFit::unrelated)
	 {
	    printf("\n  Analyzing chromosome X ... \n");
	    fprintf(log,"\n  Analyzing chromosome X ... \n");
	    double tol = 0.0001;

	    //setup kinshipX if have not done so
	    if(FastTransform::empKin)
	    {
	       kin_emp.SetupEmpKinX(ped,trans.genotypedSamplePED,trans.genotypedSampleVCF,trans.samplePEDIDHash,checkData.skippedSNPs,log);
	    }

	    //Then use different engines to fit LMM model.
	    if(FastFit::separateX)
	    {
	       if(FastTransform::pedKin)
	       {
		  engine.FastFitPolyGenicModels(ped,tol,trans,kin_emp,log,true);
	       }
	       else if(FastTransform::empKin)
	       {
		  engine.FastFitPolyGenicModels(ped,tol,trans,kin_emp,log,true);
	       }
	    }
	    else
	    {
	       if(FastTransform::pedKin)
	       {
		  //when fitting variance component model with sigmag2 and sigmag2X, together with pedigree kinship, use AutoFit.
		  engine.AutoFitLinearMixModels(ped,tol,trans,kin_emp,log,true);
	       }
	       else
	       {
		  //when fitting variance component model with sigmag2 and sigmag2X, together with empirical kinship, use AutoFit2.
		  printf("    Fitting variance component model with three components ...\n      Friendly reminder: if your sample size is large, this might take long.\n");
		  engine.AutoFitLinearMixModels2(ped,tol,trans,kin_emp,log,true);
		  printf("    done.\n");
	       }
	    }

	    //refill sigma2 vector
	    if(FastFit::separateX)
	    {
	       for(int i=0;i<trans.persons;i++)
	       {
		  sigma2[i] = engine.sigma_gXHat*trans.D[i]+engine.sigma_e2Hat;
	       }
	    }
	    else
	    {
	       for(int i=0;i<trans.persons;i++)
	       {
		  sigma2[i] = trans.D[i];
	       }
	    }
	    residuals.Dimension(trans.persons);
	    for(int i=0;i<trans.persons;i++)
	       residuals[i]=trans.UY[i]-trans.UX[i].InnerProduct(engine.betaHat);
	    //change status
	    fitXStatus=true;
	 }

	if(status==1)
//	 if(status==1 && ( current_score_var < 0 || current_score_var == pos ))
	 {
	    ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\tNA\tNA\t%d\t%d\t%d\tNA\tNA\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,hom1,het,hom2);
	    ifprintf(SCOREoutput,"\n");

	    if(recessive)
	    {
	       ifprintf(SCOREoutput_rec,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\tNA\tNA\t%d\t%d\t%d\tNA\tNA\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,hom1,het,hom2);
	       ifprintf(SCOREoutput_rec,"\n");
	    }

	    if(dominant)
	    {
	       ifprintf(SCOREoutput_dom,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d\tNA\tNA\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2);
	       ifprintf(SCOREoutput_dom,"\n");
	    }
	    continue;
	 }
	 averageAF += markerMaf;
	 marker_count++;
	 //output score statistics
	 if(FastFit::unrelated)
	 {
	    UnrelatedAssoc(SCOREoutput,SCOREoutput_rec,SCOREoutput_dom,ped,engine,trans,sigma2);
	 }
	 else
	 {
	    RelatedAssoc(SCOREoutput,SCOREoutput_rec,SCOREoutput_dom,ped,engine,trans,sigma2);
	 }

	 if(chr != initial_chr)
	 {
	    //output cov for the rest of the markers saved in genotypeAll
	    int markers = genotypeAll.rows;
	    for(int i=0;i<markers;i++)
	    {
	       initial_pos = markerInLD[0];
	       //       printf("current position: %d; initial position: %d,markerInLD length %d\t",pos,initial_pos,markerInLD.Length());
	       ifprintf(SCOREcov,"%s\t%d\t",initial_chr.c_str(),initial_pos);
	       for(int j=0;j<markerInLD.Length();j++)
	       {
		  ifprintf(SCOREcov,"%d,",markerInLD[j]);
	       }
	       ifprintf(SCOREcov,"\t");
	       for(int m=0;m<genotypeAll.rows;m++)
	       {
		  double cov = CalculateCov(engine,trans,ped,genotypeAll,sigma2,m);
		  ifprintf(SCOREcov,"%g,",cov);
	       }
	       ifprintf(SCOREcov,"\n");

	       if(recessive)
	       {
		  ifprintf(SCOREcov_rec,"%s\t%d\t",initial_chr.c_str(),initial_pos);
		  for(int j=0;j<markerInLD.Length();j++)
		  {
		     ifprintf(SCOREcov_rec,"%d,",markerInLD[j]);
		  }
		  ifprintf(SCOREcov_rec,"\t");
		  for(int m=0;m<genotypeAll_rec.rows;m++)
		  {
		     double cov = CalculateCov(engine,trans,ped,genotypeAll_rec,sigma2,m);
		     ifprintf(SCOREcov_rec,"%g,",cov);
		  }
		  ifprintf(SCOREcov_rec,"\n");
	       }

	       if(dominant)
	       {
		  ifprintf(SCOREcov_dom,"%s\t%d\t",initial_chr.c_str(),initial_pos);
		  for(int j=0;j<markerInLD.Length();j++)
		  {
		     ifprintf(SCOREcov_dom,"%d,",markerInLD[j]);
		  }
		  ifprintf(SCOREcov_dom,"\t");
		  for(int m=0;m<genotypeAll_dom.rows;m++)
		  {  
		     double cov = CalculateCov(engine,trans,ped,genotypeAll_dom,sigma2,m);
		     ifprintf(SCOREcov_dom,"%g,",cov);
		  }
		  ifprintf(SCOREcov_dom,"\n");
	       }

	       //reorganize the trackers
	       if(markerInLD.Length()>0)
	       {
		  markerInLD.Delete(0);
		  genotypeAll.DeleteRow(0);
		  if(recessive)
		     genotypeAll_rec.DeleteRow(0);
		  if(dominant)
		     genotypeAll_dom.DeleteRow(0);
	       }
	       //    printf("markerInID length after %d\n",markerInLD.Length());
	    }

	    //Start a new sequence
	    initial_chr = chr;
	    initial_pos = pos;
	    markerInLD.Clear();
	    genotypeAll.Dimension(0,0);
	    if(recessive)  
	       genotypeAll_rec.Dimension(0,0);
	    if(dominant)
	       genotypeAll_dom.Dimension(0,0);
	 }
	 else
	 {
	    while(pos - initial_pos > window && markerInLD.Length() >0 )
	    {
	       //output cov for marker at initial_pos
	       //      printf("current position: %d; initial position: %d,markerInLD length %d\t",pos,initial_pos,markerInLD.Length());
	       ifprintf(SCOREcov,"%s\t%d\t",chr.c_str(),initial_pos);
	       for(int i=0;i<markerInLD.Length();i++)
	       {
		  ifprintf(SCOREcov,"%d,",markerInLD[i]);
	       }
	       ifprintf(SCOREcov,"\t");
	       for(int m=0;m<genotypeAll.rows;m++)
	       {
		  double cov = CalculateCov(engine,trans,ped,genotypeAll,sigma2,m);
		  if ( PreMeta::Simplify ) {
		  	ifprintf(SCOREcov,"%.2e,",cov);
		  }
		  else
		  ifprintf(SCOREcov,"%g,",cov);
	       }
	       ifprintf(SCOREcov,"\n");

	       if(recessive)
	       {
		  ifprintf(SCOREcov_rec,"%s\t%d\t",chr.c_str(),initial_pos);
		  for(int i=0;i<markerInLD.Length();i++)
		  {
		     ifprintf(SCOREcov_rec,"%d,",markerInLD[i]);
		  }
		  ifprintf(SCOREcov_rec,"\t");
		  for(int m=0;m<genotypeAll_rec.rows;m++)
		  {
		     double cov = CalculateCov(engine,trans,ped,genotypeAll_rec,sigma2,m);
		     ifprintf(SCOREcov_rec,"%g,",cov);
		  }
		  ifprintf(SCOREcov_rec,"\n");
	       }

	       if(dominant)
	       {
		  ifprintf(SCOREcov_dom,"%s\t%d\t",chr.c_str(),initial_pos);
		  for(int i=0;i<markerInLD.Length();i++)
		  {
		     ifprintf(SCOREcov_dom,"%d,",markerInLD[i]);
		  }
		  ifprintf(SCOREcov_dom,"\t");
		  for(int m=0;m<genotypeAll_dom.rows;m++)
		  {
		     double cov = CalculateCov(engine,trans,ped,genotypeAll_dom,sigma2,m);
		     ifprintf(SCOREcov_dom,"%g,",cov);
		  }
		  ifprintf(SCOREcov_dom,"\n");
	       }

	       markerInLD.Delete(0);
	       genotypeAll.DeleteRow(0);
	       if(recessive)
		  genotypeAll_rec.DeleteRow(0);
	       if(dominant)
		  genotypeAll_dom.DeleteRow(0);
	       initial_pos = markerInLD[0];
	    }
	    //	    printf("markerInID length after %d\n",markerInLD.Length());
	 }
	 //if not, then update current sequence  
	 markerInLD.Push(pos);
	 initial_pos = markerInLD[0];
	 genotypeAll.GrowTo(genotypeAll.rows+1,trans.persons,0.0);
	 if(recessive)
	    genotypeAll_rec.GrowTo(genotypeAll_rec.rows+1,trans.persons,0.0);
	 if(dominant)
	    genotypeAll_dom.GrowTo(genotypeAll_dom.rows+1,trans.persons,0.0);
	 for(int i=0;i<trans.persons;i++)
	 {
	    if(FastFit::unrelated)
	       genotypeAll[genotypeAll.rows-1][i] = genotype[i];
	    else
	       genotypeAll[genotypeAll.rows-1][i] = transGeno[i];
	    if(recessive)
	    {
	       if(FastFit::unrelated)
		  genotypeAll_rec[genotypeAll_rec.rows-1][i] = genotype_rec[i];
	       else
		  genotypeAll_rec[genotypeAll_rec.rows-1][i] = transGeno_rec[i];
	    }
	    if(dominant)
	    {
	       if(FastFit::unrelated)
		  genotypeAll_dom[genotypeAll_dom.rows-1][i] = genotype_dom[i];
	       else
		  genotypeAll_dom[genotypeAll_dom.rows-1][i] = transGeno_dom[i];
	    }
	 }
      }
      //output cov for the rest of the markers saved in genotypeAll
      int markers = genotypeAll.rows;
      for(int i=0;i<markers;i++)
      {
	 initial_pos = markerInLD[0];
	 ifprintf(SCOREcov,"%s\t%d\t",chr.c_str(),initial_pos);
	 for(int j=0;j<markerInLD.Length();j++)
	 {
	    ifprintf(SCOREcov,"%d,",markerInLD[j]);
	 }
	 ifprintf(SCOREcov,"\t");
	 for(int m=0;m<genotypeAll.rows;m++)
	 {
	    double cov = CalculateCov(engine,trans,ped,genotypeAll,sigma2,m);
	    if (PreMeta::Simplify)
	    	ifprintf(SCOREcov,"%.2e,",cov);
	    else
	    	ifprintf(SCOREcov,"%g,",cov);
	 }
	 ifprintf(SCOREcov,"\n");

	 if(recessive)
	 {
	    ifprintf(SCOREcov_rec,"%s\t%d\t",chr.c_str(),initial_pos);
	    for(int j=0;j<markerInLD.Length();j++)
	    {
	       ifprintf(SCOREcov_rec,"%d,",markerInLD[j]);
	    }
	    ifprintf(SCOREcov_rec,"\t");
	    for(int m=0;m<genotypeAll_rec.rows;m++)
	    {
	       double cov = CalculateCov(engine,trans,ped,genotypeAll_rec,sigma2,m);
	       ifprintf(SCOREcov_rec,"%g,",cov);
	    }
	    ifprintf(SCOREcov_rec,"\n");
	 }

	 if(dominant)
	 {
	    ifprintf(SCOREcov_dom,"%s\t%d\t",chr.c_str(),initial_pos);
	    for(int j=0;j<markerInLD.Length();j++)
	    {
	       ifprintf(SCOREcov_dom,"%d,",markerInLD[j]);
	    }
	    ifprintf(SCOREcov_dom,"\t");
	    for(int m=0;m<genotypeAll_dom.rows;m++)
	    {
	       double cov = CalculateCov(engine,trans,ped,genotypeAll_dom,sigma2,m);
	       ifprintf(SCOREcov_dom,"%g,",cov);
	    }
	    ifprintf(SCOREcov_dom,"\n");
	 }

	 //reorganize the trackers
	 if(markerInLD.Length()>0)
	 {
	    markerInLD.Delete(0);
	    genotypeAll.DeleteRow(0);
	    if(recessive)
	       genotypeAll_rec.DeleteRow(0);
	    if(dominant)
	       genotypeAll_dom.DeleteRow(0);
	 }
      }
   }

	if(genoFromVCF)
	{
		StringArray chromosome;
		reader.open(vcfInput,header);
		reader.readVcfIndex();
		const Tabix* indexPtr = reader.getVcfIndex();
		String range_chr;
		int range_start;
		int range_end;
		if ( Region == "" ) {
			for(int i = 0; i < indexPtr->getNumRefs(); i++) {
				chromosome.Push(indexPtr->getRefName(i));
			}
		}
		else {
			printf("Restrict RMW to region %s!\n", Region.c_str());
			StringArray tf;
			tf.AddTokens(Region, ":-");
			range_chr = tf[0];
			range_start = tf[1].AsInteger();
			range_end = tf[2].AsInteger();
			if ( range_end <= range_start )
				error("Invalid range: %s!\n", Region.c_str());
			chromosome.Push( range_chr );
    	}
		reader.close();

      //loop through chromosomes by reading tabix data
		for(int i=0;i<chromosome.Length();i++)
		{
			reader.open(vcfInput,header);
			reader.readVcfIndex();
		// see if need to specify region
			if ( Region == "" )
				reader.setReadSection(chromosome[i].c_str());
			else {
				bool section_status = reader.set1BasedReadSection(chromosome[i].c_str(), range_start, range_end);
				if (!section_status)
	 				error("Unable to set vcf at range: %s\n", Region.c_str());
			}
		// read from varList if needed
			int current_score_var = -1;
			int current_cov_var = -1;
			int current_score_index = -1;
			int current_cov_index = -1;
			std::string chr_str = std::string(chromosome[i].c_str());
			if ( varListName != "" ) {
				if ( varList.find(chr_str) == varList.end() ) // no var in this chr
					continue;
				bool section_status = 0;
				for( int ii=0; ii<(int)varList[chr_str].size(); ii++ ) {
					current_score_var = varList[chr_str][0];
					current_cov_var = current_score_var;
					int last_end = varList[chr_str][ (int)varList[chr_str].size()-1 ] + window + 1;
					section_status = reader.set1BasedReadSection(chr_str.c_str(), current_score_var, last_end);
					if (section_status) {
						current_score_index = ii;
						current_cov_index = ii;
						break;
					}
	 			}
	 			if (!section_status) // no available variant in this chr
	 				continue;
			}
			
			int initial_pos = 0;

			while(reader.readRecord(record))
			{
				bool status;
				if(dosage)
					status = GetGenotypeVectorVCFDosage(trans,ped);
				else
					status = GetGenotypeVectorVCF(trans,ped,checkData,log);

	    /*
	    //TODO:TS
	    printf("genotype is:\n");
	    for(int i=0;i<trans.persons;i++)
	    printf("%g ",genotype[i]);
	    printf("\n");

	    printf("genotype_dom is:\n");
	    for(int i=0;i<trans.persons;i++)
	    printf("%g ",genotype_dom[i]);
	    printf("\n");

	    printf("residual is:\n");
	    for(int i=0;i<trans.persons;i++)
	    printf("%g ",engine.residuals_unrelated[i]);
	    printf("\n");

	    printf("sigma2 is: %g.\n",engine.sigma2Hat);
	     */

				if(chr==PreMeta::xLabel && !FastFit::unrelated && (AutoFit::fitX || FastFit::separateX) && !fitXStatus)
				{
	       printf("\n  Analyzing chromosome X ... \n");
	       fprintf(log,"\n  Analyzing chromosome X ... \n");
	       double tol = 0.0001;

	       //setup kinshipX if have not done so
	       if(FastTransform::empKin)
	       {
		  kin_emp.SetupEmpKinX(ped,trans.genotypedSamplePED,trans.genotypedSampleVCF,trans.samplePEDIDHash,checkData.skippedSNPs,log);
		  /*
		     if(AutoFit::fitX)
		     {
		     kin_emp.SetupEmpKin(ped,trans.genotypedSamplePED,trans.genotypedSampleVCF,trans.samplePEDIDHash,checkData.skippedSNPs,checkData.chromosomeVCF,log);
		     }
		   */
	       }

	       //Then use different engines to fit LMM model.
	       if(FastFit::separateX)
	       {
		  if(FastTransform::pedKin)
		  {
		     engine.FastFitPolyGenicModels(ped,tol,trans,kin_emp,log,true);
		  }
		  else if(FastTransform::empKin)
		  {
		     engine.FastFitPolyGenicModels(ped,tol,trans,kin_emp,log,true);
		  }
	       }
	       else
	       {
		  if(FastTransform::pedKin)
		  {
		     //when fitting variance component model with sigmag2 and sigmag2X, together with pedigree kinship, use AutoFit.
		     engine.AutoFitLinearMixModels(ped,tol,trans,kin_emp,log,true);
		  }
		  else
		  {
		     //when fitting variance component model with sigmag2 and sigmag2X, together with empirical kinship, use AutoFit2.
		     printf("    Fitting variance component model with three components ...\n      Friendly reminder: if your sample size is large, this might take long.\n");
		     engine.AutoFitLinearMixModels2(ped,tol,trans,kin_emp,log,true);
		     printf("    done.\n");
		  }
	       }
	       //refill sigma2 vector
	       if(FastFit::separateX)
	       {
		  for(int i=0;i<trans.persons;i++)
		  {
		     sigma2[i] = engine.sigma_gXHat*trans.D[i]+engine.sigma_e2Hat;
		  }
	       }
	       else
	       {
		  for(int i=0;i<trans.persons;i++)
		  {
		     sigma2[i] = trans.D[i];
		  }
	       }
	       residuals.Dimension(trans.persons);
	       for(int i=0;i<trans.persons;i++)
		  residuals[i]=trans.UY[i]-trans.UX[i].InnerProduct(engine.betaHat);
	       //change status
	       fitXStatus=true;
	    }

			if ( pos > current_score_var && current_score_var > 0 ) {
				printf("\nWarning: cannot find variants in vcf: %d", current_score_var);
				bool found = 0;
				for( int jj=current_score_index+1; jj<(int)varList[chr_str].size(); jj++ ) {
					if ( pos <= varList[chr_str][jj] ) {
						current_score_index = jj;
						current_cov_index = jj;
						current_score_var = varList[chr_str][jj];
						current_cov_var = varList[chr_str][jj];
						found = 1;
						break;
					}
					else
						printf(" %d", varList[chr_str][jj]);
				}
				printf(". These variants are skipped!\n");
				if ( !found ) {
					current_score_var = 2000000000;
					current_cov_var = 2000000000;
				}
			}

				if(!status)
				{
					if ( current_score_var < 0 || current_score_var == pos ) {
						if(hwe_pvalue != _NAN_ || hwe_pvalue != 6.66666e-66)
							ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\tNA\tNA\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2);
						else
		  					ifprintf(SCOREoutput,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d\tNA\tNA\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2);

	       				ifprintf(SCOREoutput,"\n");
	       				if(recessive)
	       				{
		  					if(hwe_pvalue != _NAN_ || hwe_pvalue != 6.66666e-66)
		     					ifprintf(SCOREoutput_rec,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\tNA\tNA\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2);
		  					else
		     					ifprintf(SCOREoutput_rec,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d\tNA\tNA\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2);
		  					ifprintf(SCOREoutput_rec,"\n");
	       				}

	       				if(dominant)
	       				{
		  					if(hwe_pvalue != _NAN_ || hwe_pvalue != 6.66666e-66)
								ifprintf(SCOREoutput_dom,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\tNA\tNA\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hwe_pvalue,hom1,het,hom2);
					  		else
					    		ifprintf(SCOREoutput_dom,"%s\t%d\t%s\t%s\t%d\t%g\t%g\t%g\t%g\tNA\t%d\t%d\t%d\tNA\tNA\tNA\tNA",chr.c_str(),pos,refAllele.c_str(),rareAllele.c_str(),trans.persons,founderMaf,markerMaf,mac,callRate,hom1,het,hom2);
							ifprintf(SCOREoutput_dom,"\n");
						}
					}
					if (current_score_var == pos) {
						updateScoreVar( chr_str, current_score_var, current_score_index);
					}
					continue;
				}
	    		marker_count++;
	    		averageAF += markerMaf;

	    		if(marker_count==1) 
	       			initial_pos = pos;

	    //Output score statistics here
	    		if ( current_score_var < 0 || pos == current_score_var ) {
		    		if(FastFit::unrelated)
		    		{
	    	   			UnrelatedAssoc(SCOREoutput,SCOREoutput_rec,SCOREoutput_dom,ped,engine,trans,sigma2);
	    			}
			    	else
			    	{
			       		RelatedAssoc(SCOREoutput,SCOREoutput_rec,SCOREoutput_dom,ped,engine,trans,sigma2);
			    	}
			    	if ( current_score_var == pos )
						updateScoreVar( chr_str, current_score_var, current_score_index);
			    }

	    //check if position has passed the window
			    while(pos - initial_pos > window && markerInLD.Length() >0 )
			    {
			      while ( current_cov_var > 0 && initial_pos > current_cov_var ) { // this variant might  be skipped above
			      	updateCovVar( chr_str, current_cov_var, current_cov_index );
			      }
			      if (current_cov_var < 0 || initial_pos == current_cov_var) {
			       //output cov for marker at initial_pos
			       ifprintf(SCOREcov,"%s\t%d\t",chr.c_str(),initial_pos);
			       for(int i=0;i<markerInLD.Length();i++)
			       {
					  ifprintf(SCOREcov,"%d,",markerInLD[i]);
			       }
			       ifprintf(SCOREcov,"\t");
			       for(int m=0;m<genotypeAll.rows;m++)
			       {
						double cov = CalculateCov(engine,trans,ped,genotypeAll,sigma2,m);
						if (PreMeta::Simplify)
							ifprintf(SCOREcov,"%.2e,",cov);
						else
		  					ifprintf(SCOREcov,"%g,",cov);
	       			}
			    	ifprintf(SCOREcov,"\n");

					if(recessive)
					{
						ifprintf(SCOREcov_rec,"%s\t%d\t",chr.c_str(),initial_pos);
						for(int i=0;i<markerInLD.Length();i++)
						{
							ifprintf(SCOREcov_rec,"%d,",markerInLD[i]);
						}
						ifprintf(SCOREcov_rec,"\t");
		  				for(int m=0;m<genotypeAll_rec.rows;m++)
		  				{
		     				double cov = CalculateCov(engine,trans,ped,genotypeAll_rec,sigma2,m);
		     				ifprintf(SCOREcov_rec,"%g,",cov);
		  				}
		  				ifprintf(SCOREcov_rec,"\n");
	       			}

	       			if(dominant)
	       			{
		  				ifprintf(SCOREcov_dom,"%s\t%d\t",chr.c_str(),initial_pos);
		  				for(int i=0;i<markerInLD.Length();i++)
		  				{
		     				ifprintf(SCOREcov_dom,"%d,",markerInLD[i]);
		 				}
		  				ifprintf(SCOREcov_dom,"\t");
		  				for(int m=0;m<genotypeAll_dom.rows;m++)
		  				{
		     				double cov = CalculateCov(engine,trans,ped,genotypeAll_dom,sigma2,m);
		     				ifprintf(SCOREcov_dom,"%g,",cov);
		  				}
		  				ifprintf(SCOREcov_dom,"\n");
	       			}
	       			if ( initial_pos == current_cov_var )
						updateCovVar( chr_str, current_cov_var, current_cov_index );
	       		  }

	       //reorganize the trackers
	      			markerInLD.Delete(0);
	       			initial_pos = markerInLD[0];
	       			genotypeAll.DeleteRow(0);
	       			if(recessive)
		  				genotypeAll_rec.DeleteRow(0);
	       			if(dominant)
		 				genotypeAll_dom.DeleteRow(0);
	    		}

	    //push the new marker in
	    		markerInLD.Push(pos);
	    		initial_pos = markerInLD[0];
	    		genotypeAll.GrowTo(genotypeAll.rows+1,trans.persons,0.0);

	    		for(int i=0;i<trans.persons;i++)
	    		{
	       			if(FastFit::unrelated)
		  				genotypeAll[genotypeAll.rows-1][i] = genotype[i];
	       			else
		  				genotypeAll[genotypeAll.rows-1][i] = transGeno[i];
	    		}
	    		if(recessive)
	    		{
	       			genotypeAll_rec.GrowTo(genotypeAll_rec.rows+1,trans.persons,0.0);

	       			for(int i=0;i<trans.persons;i++)
	       			{
		  				if(FastFit::unrelated)
		     				genotypeAll_rec[genotypeAll_rec.rows-1][i] = genotype_rec[i];
		  				else
		     				genotypeAll_rec[genotypeAll_rec.rows-1][i] = transGeno_rec[i];
	       			}
	    		}
	    		if(dominant)
	    		{
	       			genotypeAll_dom.GrowTo(genotypeAll_dom.rows+1,trans.persons,0.0);

	       			for(int i=0;i<trans.persons;i++)
	       			{
		  				if(FastFit::unrelated)
		     				genotypeAll_dom[genotypeAll_dom.rows-1][i] = genotype_dom[i];
		  				else
		     				genotypeAll_dom[genotypeAll_dom.rows-1][i] = transGeno_dom[i];
	       			}
	    		}
	 		}
	 //output the LD matrix for the rest of the markers in genotypeAll
	 		int markers = genotypeAll.rows;

			for(int i=0;i<markers;i++)
	 		{
	 		  while ( current_cov_var > 0 && initial_pos > current_cov_var ) { // this variant might  be skipped above
				updateCovVar( chr_str, current_cov_var, current_cov_index );
			  }
	 		  if (current_cov_var < 0 || initial_pos == current_cov_var) {
	    		ifprintf(SCOREcov,"%s\t%d\t",chr.c_str(),initial_pos);
	    		for(int j=0;j<markerInLD.Length();j++)
	    		{
	       			ifprintf(SCOREcov,"%d,",markerInLD[j]);
	    		}
	    		ifprintf(SCOREcov,"\t");
	    		for(int m=0;m<genotypeAll.rows;m++)
	    		{
	       			double cov = CalculateCov(engine,trans,ped,genotypeAll,sigma2,m);
	       			if (PreMeta::Simplify)
	       				ifprintf(SCOREcov,"%.2e,",cov);
	       			else
	       				ifprintf(SCOREcov,"%g,",cov);
	    		}
	    		ifprintf(SCOREcov,"\n");
	    	  }
	    		genotypeAll.DeleteRow(0);

	    		if(recessive)
	    		{
	    		  if (current_cov_var < 0 || initial_pos == current_cov_var) {
	       			ifprintf(SCOREcov_rec,"%s\t%d\t",chr.c_str(),initial_pos);
	       			for(int j=0;j<markerInLD.Length();j++)
	       			{
		  				ifprintf(SCOREcov_rec,"%d,",markerInLD[j]);
	       			}
	       			ifprintf(SCOREcov_rec,"\t");
	       			for(int m=0;m<genotypeAll_rec.rows;m++)
	       			{
		  				double cov = CalculateCov(engine,trans,ped,genotypeAll_rec,sigma2,m);
		  				ifprintf(SCOREcov_rec,"%g,",cov);
	       			}
	       			ifprintf(SCOREcov_rec,"\n");
	       		  }
	       			genotypeAll_rec.DeleteRow(0);
	    		}

	    		if(dominant)
	    		{
	    		  if (current_cov_var < 0 || initial_pos == current_cov_var) {
	       			ifprintf(SCOREcov_dom,"%s\t%d\t",chr.c_str(),initial_pos);
	       			for(int j=0;j<markerInLD.Length();j++)
	       			{
		  				ifprintf(SCOREcov_dom,"%d,",markerInLD[j]);
	       			}
	       			ifprintf(SCOREcov_dom,"\t");
	       			for(int m=0;m<genotypeAll_dom.rows;m++)
	       			{
		  				double cov = CalculateCov(engine,trans,ped,genotypeAll_dom,sigma2,m);
		  				ifprintf(SCOREcov_dom,"%g,",cov);
	       			}
	       			ifprintf(SCOREcov_dom,"\n");
	       		  }
	       			genotypeAll_dom.DeleteRow(0);
	    		}
	    		
	    		if ( initial_pos == current_cov_var )
	    			updateCovVar( chr_str, current_cov_var, current_cov_index );
	
	    	//reorganize the trackers
	    		markerInLD.Delete(0);
	    		initial_pos = markerInLD[0];
	 		}
	 		reader.close();
		}
	}

   if(averageAF/marker_count>0.5)
   {
      printf("\n\nWARNING: Please check for allele flips, the reference allele frequency in your sample is <0.5 on average!\n\n");
      fprintf(log,"\n\nWARNING: Please check for allele flips, the reference allele frequency in your sample is <0.5 on average!\n\n");
   }
   if(chisqred.Length()>0)
   {
      chisqred.Sort();
      //lambda = chisqred[0.5]/0.4549364;
   }

   String filename,filename_rec,filename_dom;
   GetRealPrefix(filename);

   if(recessive)
      filename_rec =filename + "recessive.plots.pdf";
   if(dominant)
      filename_dom =filename + "dominant.plots.pdf";
   filename += "additive.plots.pdf";

   String model;
   model = "additive";
   GeneratePlots(filename,pvalueAll,chr_plot,pos_plot,group,SCOREoutput,chisq_before_GCcorrect,maf_GC,model);

   if(recessive)
   {
      model = "recessive";
      GeneratePlots(filename_rec,pvalueAll_rec,chr_plot_rec,pos_plot_rec,group,SCOREoutput_rec,chisq_before_GCcorrect_rec,maf_GC_rec,model);
   }

   if(dominant)
   {
      model = "dominant";
      GeneratePlots(filename_dom,pvalueAll_dom,chr_plot_dom,pos_plot_dom,group,SCOREoutput_dom,chisq_before_GCcorrect_dom,maf_GC_dom,model);
   }
   if(malehwewarning>10)
      printf("\n  WARNING: Morer than 10 warnings were generated for HWE pvalue calculation. Please check log for details.\n");
   if(warnings>10)
      printf("\n  WARNING: More than 10 warnings were generated for male heterozygous genotyes in nonPAR region. Please check log for details.\n");

   if(FastFit::separateX)
   {
      ifprintf(SCOREoutput,"#Xheritabilit is %.2f.\n",engine.Xheritability);
   }
   if(AutoFit::fitX && !FastFit::separateX)
   {
      ifprintf(SCOREoutput,"#Autosomal heritabilit is %.2f.\n",engine.heritabilityAuto);
      ifprintf(SCOREoutput,"#Xheritabilit is %.2f.\n",engine.Xheritability);
   }

   printf("completed.\n\n");
   fprintf(log,"completed.\n\n");
}

void PreMeta::GeneratePlots(String & filename_,Vector & pvalueAll_,StringArray & chr_plot_,Vector & pos_plot_,GroupFromAnnotation & group,IFILE & SCOREoutput_,Vector & chisq_before_GCcorrect_,Vector & maf_GC_,String & model)
{
   Vector pvalue1, pvalue5;
   PDF pdf;
   WritePDF writepdf;
   StringArray geneLabel;
   if(group.labelHits)
   {
      String target_chr="";
      int target_pos=0;
      int target=0;
      double target_pvalue = _NAN_;
      for(int lb=0;lb<pvalueAll_.Length();lb++)
      {
	 if(pvalueAll_[lb] < 0.05/pvalueAll_.Length())
	 {
	    //if this locus is annotated already
	    if(target_chr==chr_plot_[lb] && target_pos>pos_plot_[lb]-1000000)
	    {
	       //if this is the higher peak then update the annotation of this locus
	       if(pvalueAll_[lb]<target_pvalue)
	       {
		  String current_anno = group.AnnotateSingleVar(chr_plot_[lb],pos_plot_[lb]);
		  int distance = 1000;
		  while(current_anno=="" && distance<=1000000)
		  {
		     current_anno = group.AnnotateSingleVar(chr_plot_[lb],pos_plot_[lb]+distance);
		     if(current_anno=="")
			current_anno = group.AnnotateSingleVar(chr_plot_[lb],pos_plot_[lb]-distance);
		     distance += 1000;
		  }
		  if(geneLabel[target]!="" && current_anno!="" && geneLabel[target].Find(current_anno)==-1)
		     current_anno = geneLabel[target] + "/" + current_anno;
		  if(geneLabel[target]!="" && current_anno=="")
		     current_anno = geneLabel[target];

		  geneLabel.Push(current_anno);
		  geneLabel[target]="";
		  target = geneLabel.Length()-1;
		  target_pvalue = pvalueAll_[lb];
		  target_pos = pos_plot_[lb];
		  target_chr = chr_plot_[lb];
	       }
	       //otherwise, leave the original annotation of this locus
	       else
		  geneLabel.Push("");
	    }
	    else
	    {
	       geneLabel.Push(group.AnnotateSingleVar(chr_plot_[lb],pos_plot_[lb]));
	       target_chr = chr_plot_[lb];
	       target_pos = pos_plot_[lb];
	       target_pvalue = pvalueAll_[lb];
	       target = geneLabel.Length()-1;
	    }
	 }
	 else
	    geneLabel.Push("");
      }
   }
   //calculate GC
   Vector tmp;
   tmp.Copy(chisq_before_GCcorrect_);
   tmp.Sort();
   double GC=_NAN_;
   if(tmp.Length()>0)
      GC = tmp[0.5]/0.456;
   else
      printf("  WARNING: GC for %s model is not valid. One possible reason is that all of your variants are monomorphic.\n",model.c_str());

   pdf.OpenFile(filename_);
   String title = FastFit::traitName;
   String demo = "GC = ";
   demo+= GC;
   if(tmp.Length()>0)
   {
      printf("  Genomic Control for %s model is: %g\n",model.c_str(),GC);
      ifprintf(SCOREoutput_,"#Genomic control for %s model is: %g\n",model.c_str(),GC);
   }

   for(int i=0;i<pvalueAll_.Length();i++)
   {
      if(maf_GC_[i]<0.01)
	 pvalue1.Push(pvalueAll_[i]);
      if(maf_GC_[i]<0.05)
	 pvalue5.Push(pvalueAll_[i]);
   }
   if(group.labelHits)
      writepdf.Draw(pdf,geneLabel,pvalueAll_,pvalue1,pvalue5,chr_plot_,pos_plot_,title,demo,true);
   else
      writepdf.Draw(pdf,geneLabel,pvalueAll_,pvalue1,pvalue5,chr_plot_,pos_plot_,title,demo,false);

   if(correctGC)
   {
      Vector pvalue_GCcorrected,pvalue1_GCcorrected,pvalue5_GCcorrected;
      chisq_before_GCcorrect_.Multiply(1.0/GC);
      for(int i=0;i<chisq_before_GCcorrect_.Length();i++)
      {
	 double p = pchisq(chisq_before_GCcorrect_[i],1,0,0);
	 pvalue_GCcorrected.Push(p);
	 if(maf_GC_[i]<0.01)
	    pvalue1_GCcorrected.Push(p);
	 if(maf_GC_[i]<0.05)
	    pvalue5_GCcorrected.Push(p);
      }
      title += " (GC corrected)";
      if(group.labelHits)
	 writepdf.Draw(pdf,geneLabel,pvalue_GCcorrected,pvalue1_GCcorrected,pvalue5_GCcorrected,chr_plot_,pos_plot_,title,demo,true);
      else
	 writepdf.Draw(pdf,geneLabel,pvalue_GCcorrected,pvalue1_GCcorrected,pvalue5_GCcorrected,chr_plot_,pos_plot_,title,demo,false);
   }
   pdf.CloseFile();
}

bool PreMeta::GetGenotypeVectorVCFDosage(FastTransform & trans, Pedigree & ped)
{
   bool status = false;

   founderMaf = markerMaf = 0.0;
   rareAllele = record.getAltStr();
   refAllele = record.getRefStr();
   chr = record.getChromStr();
   pos = record.get1BasedPosition();

   if(refAllele=="."||rareAllele==".")
      return status;

   int n=0,nf=0;
   double fmaf = 0.0, maf=0.0;
   founderMaf = markerMaf = 0.0;
   VcfRecordGenotype & genoInfo = record.getGenotypeInfo();
   /*
      bool XStatus = false;
      if(chr==xLabel && pos>=Xstart && pos<=Xend)
      XStatus=true;
    */
   //printf("Now working on marker %d.\n",markerIdx);
   genotype.resize(trans.persons);
   for(int i=0;i<trans.persons;i++)
      genotype[i]=-1.0;

   mac=0.0;

   //Loop through samples according to trans.pPheno 
   //read in dosage into genotype vector
   int idx=0,nmiss=0;
   double mean = 0.0;
   for (int f = 0; f < ped.familyCount; f++)
   {
      for(int j=0;j<trans.pPheno[f].Length();j++)
      {
	 String sample;

	 if(trans.mergedVCFID)
	    sample = ped[trans.pPheno[f][j]].famid + "_" + ped[trans.pPheno[f][j]].pid;
	 else
	    sample = ped[trans.pPheno[f][j]].pid;

	 //get the sample number by sample ID from hash
	 int s = trans.sampleVCFIDHash.Integer(sample);
	 if(s==-1)
	    error("ERROR! Check sample from VCF and PED file.\n",sample.c_str());
	 if(s!=-1)
	 {
	    const std::string * geno = genoInfo.getString(dosageFlag.c_str(),s);
	    if(strcmp((*geno).c_str(),".")>0)
	    {
	       float dose = atof((*geno).c_str());
	       genotype[idx] = dose;
	       /*
		  if(XStatus && ped[trans.samplePEDIDHash.Integer(sample)].sex==maleLabel && dose >0)
		  {
		  genotype[idx] = 2;
		  }
		*/
	       n++;
	       maf += genotype[idx];
	       int founder = trans.foundersHash.Integer(sample);
	       if(founder!=-1)
	       {
		  fmaf+=dose;
		  nf++;
	       }
	    }
	    else
	    {
	       nmiss++;
	    }
	 }
	 idx++;
      }
   }
   mac = maf;
   maf /= 2.0*n;
   fmaf /= 2.0*nf;
   if(fmaf==0.0 && maf>0.0)
      fmaf = 1.0/(2.0*nf);
   founderMaf = fmaf;
   markerMaf = maf;

   mean = 2.0 * markerMaf;
   if(maf==0.0 || 1.0-maf==0.0)
      return status;

   if(recessive)
   {
      double avg =0.0;
      genotype_rec.resize(trans.persons);
      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype[g]>=1.5)
	 {
	    genotype_rec[g]=1;
	    avg += 1;
	 }
	 else if(genotype[g]==-1)
	    genotype_rec[g]=-1;
	 else
	    genotype_rec[g]=0;
      }
      //center the genotypes
      avg /= n;
      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype_rec[g]==-1)
	    genotype_rec[g] = 0.0;
	 else
	    genotype_rec[g] -= avg;
      }
   }

   if(dominant)
   {
      double avg =0.0;
      genotype_dom.resize(trans.persons);
      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype[g]>0.5)
	 {
	    genotype_dom[g]=1;
	    avg += 1;
	 }
	 else if(genotype[g]==-1)
	    genotype_dom[g]=-1;
	 else
	    genotype_dom[g]=0;
      }
      //center the genotypes
      avg /= n;
      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype_dom[g]==-1)
	    genotype_dom[g] = 0.0;
	 else
	    genotype_dom[g] -= avg;
      }
   }

   status = true;
   for(int i=0;i<trans.persons;i++)
   {
      if(genotype[i]==-1.0)
	 genotype[i] = 0.0;
      else
	 genotype[i] -= mean;
   }
   callRate = 1.0-(double) nmiss/trans.persons;
   hwe_pvalue = 1.0;
   het=0; hom1=0; hom2=0;

   return status;
}

bool PreMeta::GetGenotypeVectorVCF(FastTransform & trans, Pedigree & ped,SanityCheck & checkData,FILE * log)
{
   bool status = false;

   founderMaf = markerMaf = 0.0;
   rareAllele = record.getAltStr();
   refAllele = record.getRefStr();
   chr = record.getChromStr();
   pos = record.get1BasedPosition();

   if(checkData.skippedSNPs.Integer(chr + ":" + pos)!=-1)
      return status;

   //Loop through samples according to trans.pPheno 
   int hwe_het=0,hwe_hom1=0,hwe_hom2=0;
   het=0;hom1=0;hom2=0;
   int numMiss = 0;
   double mean = 0.0;
   int n=0,nf=0;
   double fmaf=0.0,maf=0.0;
   mac=0.0;
   //hwe_pvalue=_NAN_;

   bool XStatus = false;
   if(chr==xLabel && pos>=Xstart && pos<=Xend)
      XStatus=true;

   genotype.resize(trans.persons);
   for(int g=0;g<genotype.size();g++)
      genotype[g] = 0.0;

   int idx=0;
   for (int f = 0; f < ped.familyCount; f++)
   {
      for(int j=0;j<trans.pPheno[f].Length();j++)
      {
	 String sample;

	 if(trans.mergedVCFID)
	    sample = ped[trans.pPheno[f][j]].famid + "_" + ped[trans.pPheno[f][j]].pid;
	 else
	    sample = ped[trans.pPheno[f][j]].pid;

	 bool miss_stat = false;

	 //get the sample number by sample ID from hash
	 int s = trans.sampleVCFIDHash.Integer(sample);
	 if(s==-1)
	    error("ERROR! Check sample %s in VCF and PED file.\n",sample.c_str());
	 int numGTs = record.getNumGTs(s);
	 int sum = 0;
	 int founder = trans.foundersHash.Integer(sample);
	 if(XStatus && ped[trans.samplePEDIDHash.Integer(sample)].sex==maleLabel)
	 {
	    if(numGTs==1)
	    {
	       int a = record.getGT(s,j);
	       if(a==VcfGenotypeSample::MISSING_GT)
		  miss_stat=true;
	       else
	       {
		  n++;
		  if(founder!=-1)
		  {
		     if(a>1)
			fmaf++;
		     nf++;
		  }
		  if(a>1)
		  {
		     maf++;
		  }
		  sum+=a;
	       }
	    }
	    else
	    {
	       int malemissing=0;
	       int a;
	       for(int j = 0; j < numGTs; j++)
	       {
		  a=record.getGT(s,j);
		  if(a==VcfGenotypeSample::MISSING_GT)
		     malemissing++;
		  else
		     sum += a;
	       }
	       if(malemissing == numGTs)
	       {
		  miss_stat=true;
	       }
	       if(!miss_stat)
	       {
		  if(founder!=-1)
		  {
		     nf++;
		     if(sum>1)
			fmaf+=1;
		  }
		  n++;
		  if(sum>1)
		     maf++;
	       }

	       //if this nonPAR variant in this male is heterozygous, then set the genotype to missing.
	       if(sum==1)
	       {
		  miss_stat=true;
		  warnings++;
		  fprintf(log,"WARNING: variant %s:%d for sample %s is heterozygous and treated as missing.\n",chr.c_str(),pos,sample.c_str());
		  if(warnings<10)
		  {
		     printf("WARNING: variant %s:%d for sample %s is heterozygous and treated as missing.\n",chr.c_str(),pos,sample.c_str());
		  }
	       }
	    }
	    if(miss_stat)
	    {
	       numMiss++;
	       genotype[idx] = -1.0;
	       idx++;
	    }
	 }
	 else
	 {
	    for(int j = 0; j < numGTs; j++)
	    {
	       int a = record.getGT(s,j); 
	       if(a==VcfGenotypeSample::MISSING_GT)
	       {
		  numMiss++;
		  genotype[idx] = -1.0;
		  //printf("%s\t%g\n",sample.c_str(),genotype[idx]);
		  idx++;
		  miss_stat=true;
		  break;
	       }
	       sum += a;
	    }
	    if(!miss_stat)
	    {
	       n+=2;
	       maf+=sum;
	       if(founder!=-1) 
	       {
		  fmaf+=sum;
		  nf+=2;
	       }
	    }
	 }
	 if(miss_stat)
	 {
	    continue;
	 }

	 genotype[idx] = sum;
	 if(sum==1)
	    het++;
	 else if(sum==0)
	    hom1++;
	 else if(sum==2)
	    hom2++;

	 if(XStatus && ped[trans.samplePEDIDHash.Integer(sample)].sex==femaleLabel)
	 {
	    if(sum==1)
	       hwe_het++;
	    if(sum==0)
	       hwe_hom1++;
	    if(sum==2)
	       hwe_hom2++;
	 }
	 idx++;
      }
   }

   mac = maf;
   maf /= n;
   fmaf /= nf;
   if(fmaf==0.0 && maf>0.0)
      fmaf = 1.0/nf;

   callRate = 1.0-(double) numMiss/trans.persons;

   hwe_pvalue = 0.0;
   //Get hwe_pvalue using the entire sample
   if(het==0 && hom1==0 && hom2==0)
   {
      //printf("marker %s:%d has problem.\n",chr.c_str(),pos);
      refAllele = rareAllele = ".";
      return status;
   }

   if(XStatus)
   {
      if(hwe_het + hwe_hom1 + hwe_hom2 ==0)
      {
	 if(malehwewarning<10)
	    printf("HWE pvalue is NA for variant %s:%d, because no female sample is genotyped at this site, which belongs to nonPAR region in chr X.\n",chr.c_str(),pos);
	 fprintf(log, "HWE pvalue is NA for variant %s:%d, because no female sample is genotyped at this site, which belongs to nonPAR region in chr X.\n",chr.c_str(),pos);
	 hwe_pvalue = _NAN_;
	 malehwewarning++;
      }
      else
      {
	 hwe_pvalue = SNPHWE(hwe_het,hwe_hom1,hwe_hom2);
      }
   }
   else
   {
      hwe_pvalue = SNPHWE(het,hom1,hom2);
   }
   founderMaf = fmaf;
   markerMaf = maf;

   if(maf==0.0 || 1.0-maf==0.0)
      return status;

   status = true;

   mean = 2.0 * markerMaf;

   n /= 2;

   if(recessive)
   {
      double avg =0.0;
      genotype_rec.resize(trans.persons);
      if(XStatus)
      {
	 int n_=0;
	 for(int g=0;g<trans.persons;g++)
	 {
	    if(genotype[g]==2)
	    {
	       genotype_rec[g]=1;
	       avg += 1;
	       n_++;
	    }
	    else if(genotype[g]==-1)
	       genotype_rec[g]=-1;
	    else
	    {
	       genotype_rec[g]=0;
	       n_++;
	    }
	 }
	 avg /= n_;
      }
      else
      {
	 for(int g=0;g<trans.persons;g++)
	 {
	    if(genotype[g]==2)
	    {
	       genotype_rec[g]=1;
	       avg += 1;
	    }
	    else if(genotype[g]==-1)
	       genotype_rec[g]=-1;
	    else
	       genotype_rec[g]=0;
	 }
	 avg /= n;
      }
      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype_rec[g]==-1)
	    genotype_rec[g] = 0.0;
	 else
	    genotype_rec[g] -= avg;
      }
   }

   if(dominant)
   {
      genotype_dom.resize(trans.persons);
      double avg =0.0;
      if(XStatus)
      {
	 int n_=0;
	 for(int g=0;g<trans.persons;g++)
	 {
	    if(genotype[g]>=1)
	    {
	       genotype_dom[g]=1;
	       avg += 1;
	       n_++;
	    }
	    else if(genotype[g]==-1)
	       genotype_dom[g]=-1;
	    else
	    {
	       genotype_dom[g]=0;
	       n_++;
	    }
	 }
	 //center the genotypes
	 avg /= n_;
      }
      else
      {
	 for(int g=0;g<trans.persons;g++)
	 {
	    if(genotype[g]>=1)
	    {
	       genotype_dom[g]=1;
	       avg += 1;
	    }
	    else if(genotype[g]==-1)
	       genotype_dom[g]=-1;
	    else
	       genotype_dom[g]=0;
	 }
	 //center the genotypes
	 avg /= n;
      }
      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype_dom[g]==-1)
	    genotype_dom[g] = 0.0;
	 else
	    genotype_dom[g] -= avg;
      }
   }

   if(XStatus)
   {
      mean =0.0;n=0;
      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype[g]!=-1.0)
	 {
	    mean += genotype[g];
	    n++;
	 }
      }
      mean /= n;
   }

   for(int g=0;g<trans.persons;g++)
   {
      if(genotype[g]==-1.0)
	 genotype[g] = 0.0;
      else
	 genotype[g] -= mean;
   }
   //printf("centered sum is: %g\n",genotype.Sum());
   return status;
}

int PreMeta::GetGenotypeVectorPED(FastTransform & trans, Pedigree & ped,int markerid,FILE * log)
{
   //status=0: polymorphic
   //status=1: monomorphic or allele count greater than 2 
   //status=2: marker name from DAT file is incorrect and skip this marker
   int status = 1;

   StringArray markerName;
   markerName.AddTokens(ped.markerNames[markerid],":");
   if(markerName.Length()<=1)
   {
      printf("WARNING: marker name %s is not in CHR:POS format; this marker will be skipped from analysis.\n",ped.markerNames[markerid].c_str());
      fprintf(log,"WARNING: marker name %s is not in CHR:POS format; this marker will be skipped from analysis.\n",ped.markerNames[markerid].c_str());
      status=2;
      return status;
   }

   chr = markerName[0];
   pos = atoi(markerName[1].c_str());
   markerName.Clear();

   founderMaf = markerMaf = 0.0;
   refAllele = rareAllele=".";

   genotype.resize(trans.persons);
   for(int g=0;g<genotype.size();g++)
      genotype[g] = 0.0;
   het=0;hom1=0;hom2=0;
   int hwe_het=0,hwe_hom1=0,hwe_hom2=0;
   int numMiss=0;
   double mean =0.0;
   int n=0,nf=0;
   double maf=0.0,fmaf=0.0;
   callRate = 0.0;
   hwe_pvalue = 0.0;
   mac =0.0;

   int alleleCount = CountAlleles(ped, markerid);

   if(alleleCount==1)
   {
      rareAllele = ped.GetMarkerInfo(markerid)->alleleLabels[1];
      refAllele = rareAllele;
      return status;
   }

   if(alleleCount!=2)
      return status;

   refAllele = ped.GetMarkerInfo(markerid)->alleleLabels[1];
   rareAllele = ped.GetMarkerInfo(markerid)->alleleLabels[2];

   int idx=0;
   bool XStatus = false;
   if(chr==xLabel && pos>=Xstart && pos<=Xend)
      XStatus=true;

   for(int f=0;f<ped.familyCount;f++)
   {
      for(int p=0;p<trans.pPheno[f].Length();p++)
      {
	 int founder = trans.foundersHash.Integer(ped[trans.pPheno[f][p]].pid);

	 if(!ped[trans.pPheno[f][p]].isGenotyped(markerid))
	 {
	    numMiss++;
	    genotype[idx]=-1.0;
	    idx++;
	    continue;
	 }
	 else if(ped[trans.pPheno[f][p]].markers[markerid].isHeterozygous())
	 {
	    if(XStatus && ped[trans.pPheno[f][p]].sex==maleLabel)
	    {
	       numMiss++;
	       genotype[idx]=-1.0;
	    }
	    else
	    {
	       het++;
	       genotype[idx]=1.0 ;
	       maf++;
	       n += 2;
	       if(founder != -1)
	       {
		  fmaf++;
		  nf +=2;
	       }
	    }
	 }
	 else if(ped[trans.pPheno[f][p]].markers[markerid].isHomozygous())
	 {
	    if(ped[trans.pPheno[f][p]].markers[markerid].hasAllele(2))
	    {
	       hom2++;
	       genotype[idx]=2.0;
	       if(XStatus && ped[trans.pPheno[f][p]].sex==maleLabel)
	       {
		  maf++;
		  n++;
		  if(founder != -1)
		  {
		     fmaf++;
		     nf++;
		  }
	       }
	       else
	       {
		  maf+=2;
		  n+=2;
		  if(founder != -1)
		  {
		     fmaf+=2;
		     nf+=2;
		  }
	       }
	    }
	    else if(ped[trans.pPheno[f][p]].markers[markerid].hasAllele(1))
	    {
	       hom1++;
	       genotype[idx]=0.0;
	       if(XStatus && ped[trans.pPheno[f][p]].sex==maleLabel)
	       {
		  n++;
		  if(founder != -1)
		  {
		     nf++;
		  }
	       }
	       else
	       {
		  n+=2;
		  if(founder != -1)
		  {
		     nf +=2;
		  }
	       }
	    }
	 }
	 if(XStatus && ped[trans.pPheno[f][p]].sex==femaleLabel)
	 {
	    if(ped[trans.pPheno[f][p]].markers[markerid].isHeterozygous())
	       hwe_het++;
	    if(ped[trans.pPheno[f][p]].markers[markerid].isHomozygous())
	    {
	       if(ped[trans.pPheno[f][p]].markers[markerid].hasAllele(2))
		  hwe_hom2++;
	       if(ped[trans.pPheno[f][p]].markers[markerid].hasAllele(1))
		  hwe_hom1++;
	    }
	 }

	 idx++;
      }
   }

   mac = maf;

   maf /= n;
   fmaf /= nf;

   if(fmaf==0.0 && maf>0.0)
      fmaf = 1.0/nf;

   callRate = 1.0-(double) numMiss/trans.persons;
   hwe_pvalue=0.0;
   //Get hwe_pvalue using the entire sample
   if(het==0 && hom1==0 && hom2==0)
   {
      //printf("marker %s:%d has problem.\n",chr.c_str(),pos);
      refAllele = rareAllele = ".";
      status=1;
      return status;
   }
   if(XStatus)
   {
      if(hwe_het + hwe_hom1 + hwe_hom2 ==0)
      {
	 if(malehwewarning<10)
	 {
	    printf("HWE pvalue is NA for variant %s:%d, because no female sample is genotyped at this site, which belongs to nonPAR region in chr X.\n",chr.c_str(),pos);
	    fprintf(log, "HWE pvalue is NA for variant %s:%d, because no female sample is genotyped at this site, which belongs to nonPAR region in chr X.\n",chr.c_str(),pos);
	 }
	 hwe_pvalue = _NAN_;
	 malehwewarning++;
      }
      else
	 hwe_pvalue = SNPHWE(hwe_het,hwe_hom1,hwe_hom2);
   }
   else
      hwe_pvalue = SNPHWE(het,hom1,hom2);

   founderMaf = fmaf;
   markerMaf = maf;

   if(maf==0.0 || 1.0-maf==0.0)
      return status;

   status = 0;

   /*
      if(markerMaf > 0.5)
      {
      String tmp;
      tmp=refAllele;
      refAllele = rareAllele;
      rareAllele = tmp;
      markerMaf = 1.0-markerMaf;
      founderMaf = 1.0-founderMaf;
      for(int g=0;g<trans.persons;g++)
      {
      if(genotype[g]!=-1.0)
      genotype[g] = 2.0 -genotype[g];
      }
      int tmp_hom2 = hom2;
      hom2 = hom1;
      hom1 = tmp_hom2;
      hwe_pvalue = SNPHWE(het,hom1,hom2);
      }
    */
   mean = 2.0 * markerMaf;
   n /= 2;

   if(recessive)
   {
      double avg =0.0;
      genotype_rec.resize(trans.persons);
      if(XStatus)
      {
	 int n_=0;
	 for(int g=0;g<trans.persons;g++)
	 {
	    if(genotype[g]==2)
	    {
	       genotype_rec[g]=1;
	       avg += 1;
	       n_++;
	    }
	    else if(genotype[g]==-1)
	       genotype_rec[g]=-1;
	    else
	    {
	       genotype_rec[g]=0;
	       n_++;
	    }
	 }
	 avg /= n_;
      }
      else
      {
	 for(int g=0;g<trans.persons;g++)
	 {
	    if(genotype[g]==2)
	    {
	       genotype_rec[g]=1;
	       avg += 1;
	    }
	    else if(genotype[g]==-1)
	       genotype_rec[g]=-1;
	    else
	       genotype_rec[g]=0;
	 }
	 avg /= n;
      }

      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype_rec[g]==-1)
	    genotype_rec[g] = 0.0;
	 else
	    genotype_rec[g] -= avg;
      }
   }
   if(dominant)
   {
      genotype_dom.resize(trans.persons);
      double avg =0.0;
      if(XStatus)
      {
	 int n_=0;
	 for(int g=0;g<trans.persons;g++)
	 {
	    if(genotype[g]>=1)
	    {
	       genotype_dom[g]=1;
	       avg += 1;
	       n_++;
	    }
	    else if(genotype[g]==-1)
	       genotype_dom[g]=-1;
	    else
	    {
	       genotype_dom[g]=0;
	       n_++;
	    }
	 }
	 //center the genotypes
	 avg /= n_;
      }
      else
      {
	 for(int g=0;g<trans.persons;g++)
	 {
	    if(genotype[g]>=1)
	    {
	       genotype_dom[g]=1;
	       avg += 1;
	    }
	    else if(genotype[g]==-1)
	       genotype_dom[g]=-1;
	    else
	       genotype_dom[g]=0;
	 }
	 //center the genotypes
	 avg /= n;
      }
      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype_dom[g]==-1)
	    genotype_dom[g] = 0.0;
	 else
	    genotype_dom[g] -= avg;
      }
   }
   if(XStatus)
   {
      mean =0.0;n=0;
      for(int g=0;g<trans.persons;g++)
      {
	 if(genotype[g]!=-1.0)
	 {
	    mean += genotype[g];
	    n++;
	 }
      }
      mean /= n;
   }

   for(int g=0;g<trans.persons;g++)
   {
      if(genotype[g]==-1.0)
	 genotype[g] = 0.0;
      else
	 genotype[g] -= mean;
   }

   return status;
}

void PreMeta::PrintHeader(IFILE output,FastFit & engine,FastTransform & trans,Pedigree & ped)
{
   double totalVar=engine.sigma_g2Hat + engine.sigma_gXHat + engine.sigma_s2Hat + engine.sigma_e2Hat;
   double h2 = engine.sigma_g2Hat/totalVar;
   //double h2X = 100.0*engine.sigma_gXHat/totalVar;
   double sharedEnv =engine.sigma_s2Hat/totalVar;
   //printf("simga_g2 is: %g,sigma_e2 is: %g,betaHat is: %g\n",sigma_g2,engine.sigma_e2Hat,engine.betaHat[0]);

   /*
      if(!AutoFit::fitSharedEnvironment && !AutoFit::fitX)
      {
      fprintf(output,"\n%s mean=%6.2f, variance=%6.2f, heritability=%6.2f \n",engine.traitName.c_str(),trans.traitMean,trans.traitVar,h2);
      }
      else if(AutoFit::fitSharedEnvironment && AutoFit::fitX)
      {
      fprintf(output,"\n%s mean=%6.2f, variance=%6.2f, heritability=%6.2f, Xheritability = %6.2f, SharedEnv = %6.2f \n",engine.traitName.c_str(),trans.traitMean,trans.traitVar,h2,h2X,sharedEnv);
      }
      else if(AutoFit::fitSharedEnvironment)
      {
      fprintf(output,"\n%s mean=%6.2f, variance=%6.2f, heritability=%6.2f, SharedEnv = %6.2f \n",engine.traitName.c_str(),trans.traitMean,trans.traitVar,h2,sharedEnv);
      }
      else if(AutoFit::fitX)
      {
      fprintf(output,"\n%s mean=%6.2f, variance=%6.2f, heritability=%6.2f, Xheritability = %6.2f \n",engine.traitName.c_str(),trans.traitMean,trans.traitVar,h2,h2X);
      }
      fprintf(output,"===================================================================\n");
    */
   //fprintf(output,"\nCHR\tPOS\t#FOUNDERS\tFOUNDER_AF\tALL_AF\tINFORMATIVE_AC\tREF_ALLELE\tALT_ALLELE\tALT_ALLELE_EFFSIZE\tPVALUE\n");
   if(!FastFit::unrelated)
   {
      ifprintf(output,"##Heritability=%.2f\n",h2);
      if(AutoFit::fitSharedEnvironment)
	 ifprintf(output,"##SharedEnv=%.2f%\n",sharedEnv/100);
   }
   ifprintf(output,"#CHROM\tPOS\tREF\tALT\tN_INFORMATIVE\tFOUNDER_AF\tALL_AF\tINFORMATIVE_ALT_AC\tCALL_RATE\tHWE_PVALUE\tN_REF\tN_HET\tN_ALT\tU_STAT\tSQRT_V_STAT\tALT_EFFSIZE\tPVALUE");
   if(calculateOR)
      ifprintf(output,"\tOddsRatio");
   ifprintf(output,"\n");
}

void PreMeta::GetRealPrefix(String & realPrefix)
{
   if(outputFile != "")
   {
      if(outputFile.Last()=='.' || outputFile.Last()=='/')
	 realPrefix = outputFile + FastFit::traitName + ".";
      else
	 realPrefix = outputFile + "." + FastFit::traitName + ".";
   }
   else
      realPrefix = FastFit::traitName + ".";
}


void PreMeta::setVarList()
{
	IFILE infile = ifopen( varListName, "r" );
	int prev_pos = 0;
	std::string prev_chr = "";
	while( !ifeof(infile) ) {
		String buffer;
		buffer.ReadLine(infile);
		StringArray tokens;
		tokens.AddTokens( buffer, "\t");
		if ( tokens.Length() != 2 )
			error("Incorrect format of input variant list!\n");
		std::string chr_name = std::string( tokens[0].c_str() );
		int current_pos = tokens[1].AsInteger();
		if ( chr_name.compare(prev_chr) != 0 )
			prev_pos = 0;
		else {
			if (current_pos < prev_pos)
				error("Input variant list is not sorted by coordinates!\n");
		}
		varList[chr_name].push_back(current_pos);
	}
	ifclose(infile);
}

void PreMeta::updateScoreVar( std::string & chr_str, int & current_score_var, int & current_score_index)
{
	current_score_index++;
	if (current_score_index >= (int)varList[chr_str].size())
		current_score_var = 2000000000;
	else
		current_score_var = varList[chr_str][current_score_index];
}

void PreMeta::updateCovVar( std::string & chr_str, int & current_cov_var, int & current_cov_index)
{
	current_cov_index++;
	if ( current_cov_index >= (int)varList[chr_str].size())
	    current_cov_var = 2000000000;
	else
		current_cov_var = varList[chr_str][current_cov_index];
}



