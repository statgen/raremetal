////////////////////////////////////////////////////////////////////// 
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

#ifndef __PREMETA_H__
#define __PREMETA_H__

#include "FastFit.h"
#include "VcfRecord.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include "InputFile.h"
#include "GroupFromAnnotation.h"
#include "DataQC.h"
#include <vector>
#include <map>
#include <string>

class PreMeta
{
public:
      PreMeta( FILE * logfile,IFILE & SCOREoutput1,IFILE & SCOREoutput_rec1,IFILE & SCOREoutput_dom1,IFILE & SCOREcov1,IFILE & SCOREcov_rec1,IFILE & SCOREcov_dom1);
      ~PreMeta();

      FILE * log;

      static bool vcfAnnotated;
      static bool zipOutput;
      static bool checkRef;

      //Input/Output options	
      static String vcfInput;
      static String outputFile;
      static int Xstart,Xend;
      static String xLabel;
      static String dosageFlag;
      static bool dosage;
      static bool genoFromPed;
      static bool genoFromVCF;
      static bool FounderFreq;
      static bool correctGC;
      static int window;
      static int maleLabel;
      static int femaleLabel;
      static bool recessive;
      static bool dominant;
      static bool additive;
      static bool calculateOR;
      static bool Simplify;
      static String Region;
      static String varListName;
      static bool splitUV;
      static bool newFormat;

      int warnings;
      int numFounders;
      int malehwewarning;
      VcfFileReader reader;
      VcfHeader header;
      VcfRecord record;

      //this is the position of each marker in a gene for output
      int pos,hom1,het,hom2;
      double founderMaf,markerMaf,callRate;
      String chr,refAllele,rareAllele;
      double mac,hwe_pvalue;
      Eigen::VectorXf genotype, genotype_rec, genotype_dom;
      Eigen::VectorXf transGeno, transGeno_rec,transGeno_dom;
      Vector residuals;
      Vector chisq_before_GCcorrect_rec, chisq_before_GCcorrect_dom, chisq_before_GCcorrect;
      Vector pvalueAll,pvalueAll_rec,pvalueAll_dom;
      Vector pos_plot,pos_plot_rec,pos_plot_dom;
      StringArray chr_plot,chr_plot_rec,chr_plot_dom;
      Vector maf_GC,maf_GC_rec,maf_GC_dom;
      IFILE SCOREoutput,SCOREoutput_rec,SCOREoutput_dom;
      IFILE SCOREcov,SCOREcov_rec,SCOREcov_dom;
      Matrix genotypeAll,genotypeAll_rec,genotypeAll_dom; //this matrix has the all marker genotype within a window
      IntArray markerInLD;

      Vector chisqred;
      double lambda;
      Matrix projectionMat,inv;
      std::map< std::string, std::vector<int> > varList; // variant list if needed

      void CalculateProjectionMatrix(FastTransform & trans,FastFit & engine,Vector & sigma2);

      void Run(Pedigree & ped, FastTransform & trans,FastFit & engine,GroupFromAnnotation & group,SanityCheck & checkData,KinshipEmp & kin_emp);
      void runGenoFromPed(Pedigree & ped, FastTransform & trans, FastFit & engine, SanityCheck & checkData, KinshipEmp & kin_emp, Vector & sigma2);
      void runGenoFromVcf(Pedigree & ped, FastTransform & trans, FastFit & engine, SanityCheck & checkData, KinshipEmp & kin_emp, Vector & sigma2);
      void smallSanityCheck( double averageAF, int marker_count );


      void VCFSanityCheck();
      bool GetGenotypeVectorVCFDosage(FastTransform & trans, Pedigree & ped);
      bool GetGenotypeVectorVCF(FastTransform & trans, Pedigree & ped,SanityCheck & checkData,FILE * log);
      int GetGenotypeVectorPED(FastTransform & trans, Pedigree & ped,int markerid,FILE *log);
      //bool GetGenotypeVector(FastTransform & trans, Pedigree & ped,int markerid);
      void PrintScoreFileHeader(IFILE & output,FastFit & engine,FastTransform & trans,Pedigree & ped);
      void PrintCovFileHeader( IFILE & file );
      void PrintToCovOutput( IFILE & file, Matrix & gt, String & chr, int pos, Pedigree & ped, FastTransform & trans, FastFit & engine,Vector & sigma2);
      double CalculateCovWithSigma(FastFit & engine,FastTransform & trans, Pedigree & ped, Matrix & genotypeAl,Vector & sigma2,int m);
      double CalculatePlainCov( Matrix & gt, Vector & sigma2, int m, FastTransform & trans );

      void CalculateUnrelatedAssocStats(double & effSize,double & pvalue, double & chisq, double & numerator,double & denominator, Eigen::VectorXf & genotype,FastFit & engine,FastTransform & trans,Vector & sigma2);
      void CalculateAssocStats(double & effSize,double &pvalue,double &numerator,double &denominator,double &chisq, Eigen::VectorXf & transGeno,FastFit & engine,FastTransform & trans,Vector & sigma2);

      void RelatedAssoc(IFILE SCOREoutput, IFILE SCOREoutput_rec,IFILE SCOREoutput_dom, Pedigree & ped,FastFit & engine,FastTransform & trans,Vector & sigma2);
      void UnrelatedAssoc(Pedigree & ped,FastFit & engine,FastTransform & trans,Vector & sigma2);

      void GeneratePlots(String & filename,Vector & pvalueAll,StringArray & chr_plot,Vector & pos_plot,GroupFromAnnotation & group,IFILE & SCOREoutput,Vector & chisq_before_GCcorrect,Vector & maf_GC,String & model);
      void GetRealPrefix(String & file);
      void setVarList();
      void updateScoreVar( std::string & chr_str, int & current_score_var, int & current_score_index);
      void updateCovVar( std::string & chr_str, int & current_cov_var, int & current_cov_index);
};

#endif
