////////////////////////////////////////////////////////////////////// 
// KinshipEmp.cpp 
// (c) 2012-2013 Shuang Feng, Dajiang Liu, Goncalo Abecasis
// 
// This file is distributed as part of the RareMetalWorker source code package   
// and may not be redistributed in any form, without prior written    

// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile RareMetalWorker.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Wednesday November 28, 2012
// 

#include "KinshipEmp.h"
#include "MathMatrix.h"
#include "PedigreeAlleleFreq.h"
#include "OutputKin.h"
#include "PreMeta.h"
#include <string>
#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "math.h"
#include <algorithm>
#include "InputFile.h"
#include "StringHash.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "time.h"

double KinshipEmp::q = 0.05;
double KinshipEmp::miss = 0.05;
int KinshipEmp::cpus = 1;

void KinshipEmp::CleanUpAuto()
{
   allPairs.Dimension(0,0);
}

void KinshipEmp::CleanUpX()
{
   allPairsX.Dimension(0,0);
}

void KinshipEmp::ReadEmpKin(FILE * log)
{
   StringArray tmp;

   //Get the pids from kinship matrix
   IFILE empIn;
   empIn = ifopen(FastTransform::readInEmp.c_str(),"r");

   if(empIn==NULL)
   {
      fprintf(log,"FATAL ERROR! Estiamted kinship can not be obtained from file: %s.\n",FastTransform::readInEmp.c_str());
      error("Estiamted kinship can not be obtained from file: %s.\n",FastTransform::readInEmp.c_str());
   }
   printf("  Loading empirical kinship from file %s ...",FastTransform::readInEmp.c_str());
   fprintf(log,"  Loading empirical kinship from file %s ...",FastTransform::readInEmp.c_str());
   StringArray IDs;
   while (!ifeof(empIn))
   {
      String buffer;
      buffer.ReadLine(empIn);
      tmp.Clear();
      tmp.AddTokens(buffer, SEPARATORS);
      if(tmp.Length()>0)
      {
	 for(int i=0;i<tmp.Length();i++)
	 {
	    IDs.Push(tmp[i]);
	    IDFromEmp.SetInteger(tmp[i],i);
	 }
	 break;
      }
   }
   ifclose(empIn);

   //Fill in the empirical kinship matrix
   allPairs.Dimension(tmp.Length(),tmp.Length(),0.0);

   empIn = ifopen(FastTransform::readInEmp.c_str(),"r");
   int idx=0;
   while (!ifeof(empIn))
   {
      idx++;
      String buffer;
      buffer.ReadLine(empIn);
      if(idx==1)
	 continue;
      if(buffer.Find("nan")!=-1)
      {
	 error("Estiamted kinship is not all numerical. Please check to make sure your file does not contain \"nan\"s.\n");
      }
      tmp.Clear();
      tmp.AddTokens(buffer, SEPARATORS);
      for(int i=0;i<tmp.Length();i++)
	 allPairs[idx-2][i] = tmp[i].AsDouble();
   }
   ifclose(empIn);

   //Fill in the other half of the matrix
   for(int i=0;i<allPairs.rows;++i)
      for(int j=i+1;j<allPairs.cols;++j)
	 allPairs[i][j]=allPairs[j][i];
   printf(" done.\n\n");
}

void KinshipEmp::ReadEmpKinX(FILE * log)
{
   StringArray tmp;

   //Get the pids from kinship matrix
   IFILE empIn;

   empIn = ifopen(FastTransform::readInEmpX.c_str(),"r");
   if(empIn==NULL)
      FastTransform::readInEmpX += ".gz";

   empIn = ifopen(FastTransform::readInEmpX.c_str(),"r");

   if(empIn==NULL)
   {
      fprintf(log,"Cannot open file %s for kinship X. Please check file name.\n",FastTransform::readInEmpX.c_str());
      error("Error! Cannot open file %s for kinship X. Please check file name.\n",FastTransform::readInEmpX.c_str());
   }
   printf("    Loading empirical kinshipX from file %s ...",FastTransform::readInEmpX.c_str());
   fprintf(log,"    Loading empirical kinshipX from file %s ...",FastTransform::readInEmpX.c_str());

   while (!ifeof(empIn))
   {
      String buffer;
      buffer.ReadLine(empIn);
      tmp.Clear();
      tmp.AddTokens(buffer, SEPARATORS);
      if(tmp.Length()>0)
      {
	 for(int i=0;i<tmp.Length();i++)
	    IDFromEmpX.SetInteger(tmp[i],i);
	 break;
      }
   }
   ifclose(empIn);

   //Fill in the empirical kinship matrix
   allPairsX.Dimension(tmp.Length(),tmp.Length());

   empIn = ifopen(FastTransform::readInEmpX.c_str(),"r");
   int idx=0;
   while (!ifeof(empIn))
   {
      idx++;
      String buffer;
      buffer.ReadLine(empIn);
      if(idx==1)
	 continue;
      tmp.Clear();
      tmp.AddTokens(buffer, SEPARATORS);
      for(int i=0; i<tmp.Length(); i++)
      {
	 allPairsX[idx-2][i] = tmp[i].AsDouble();
      }
   }
   ifclose(empIn);
   //Fill in the other half of the matrix
   for(int i=0;i<allPairsX.rows;++i)
      for(int j=i+1;j<allPairsX.cols;++j)
	 allPairsX[i][j]=allPairsX[j][i];

   printf(" done.\n");
   fprintf(log," done.\n");
}

void KinshipEmp::SetupEmpKin(Pedigree & ped, IntArray & genotypedSamplePED, IntArray & genotypedSampleVCF, StringIntHash & samplePEDIDHash, StringIntHash & skippedSNPs,StringArray & chromosomeVCF, FILE * log)
{
#ifdef _OPENMP
   omp_set_num_threads(cpus);
#endif
   warnings=0;
   if(FastTransform::readInEmp != "")
   {
      ReadEmpKin(log);
   }
   else
   {
      if(PreMeta::genoFromPed)
      {
	 SetupPEDAuto(ped,genotypedSamplePED,log);
      }
      else
      {
	 SetupVCFAuto(ped,genotypedSampleVCF,chromosomeVCF,log);
      }
   }
   if(warnings>10)
   {
      printf("\nThere are more than 10 warnings generated when estimating empirical kinship. Please refer to the log file for details.\n");
   }
}

void KinshipEmp::SetupEmpKinX(Pedigree & ped, IntArray & genotypedSamplePED, IntArray & genotypedSampleVCF, StringIntHash & samplePEDIDHash, StringIntHash & skippedSNPs,FILE * log)
{
   warnings = 0;
   if(FastTransform::readInEmpX != "")
   {
      ReadEmpKinX(log);
   }
   else
   {
      if(PreMeta::genoFromPed)
      {
	 SetupPEDX(ped,genotypedSamplePED,log);
      }
      else
      {
	 SetupVCFX(ped,genotypedSampleVCF,samplePEDIDHash,skippedSNPs,log);
      }
   }
   if(warnings>10)
   {
      printf("\nThere are more than 10 warnings generated when estimating empirical kinshipX. Please refer to the log file for details.\n");
   }
}

void KinshipEmp::SetupPEDAuto(Pedigree & ped,IntArray & genotypedSamplePED, FILE * log)
{
   printf("Calculating empirical kinship matrix ...");
   fflush(stdout);
   fprintf(log,"Calculating empirical kinship matrix ...");

   int N=0;
   int total_n = genotypedSamplePED.Length();

   allPairs.Dimension(total_n, total_n);
   for(int r=0;r<total_n;r++)
   {
      for(int c=r;c<total_n;c++)
	 allPairs[r][c] = 0.0;
   }

   float * X = new float [total_n];
   for(int i=0;i<total_n;i++)
      X[i] = -1.0;

   String chr;
   int pos;

   double mean=_NAN_;
   double var_inv=_NAN_;
   double maf = _NAN_;
   int n_=0,nmiss=0;
   int NMISS=total_n*miss;
   for(int i=0;i<ped.markerCount;i++)
   {
      StringArray markerName;
      markerName.AddTokens(ped.markerNames[i],":");

      if(markerName.Length()<=1)
      {
	 printf("WARNING: marker name %s in DAT file is not in CHR:POS format and will be skipped from calculating empirical kinsihp.\n",ped.markerNames[i].c_str());
	 fprintf(log,"WARNING: marker name %s in DAT file is not in CHR:POS format and will be skipped from calculating empirical kinship.\n",ped.markerNames[i].c_str());
	 continue;
      }

      chr = markerName[0];
      pos = atoi(markerName[1].c_str());

      if(chr==PreMeta::xLabel && pos>=PreMeta::Xstart && pos<=PreMeta::Xend)
	 continue;

      //X.Clear();
      bool missingStatus=false;
      int alleleCount = CountAlleles(ped, i);
      //if no one is genotyped, mornomorphic, or marker is not biallelic
      if(alleleCount!=2)
	 continue;

      mean=0.0;

      //now fill in the genotypes
      n_=0;nmiss=0;
      //only generate matrix for genotyped individuals
      for(int p=0;p<total_n;p++)
      {
	 if(!ped[genotypedSamplePED[p]].isGenotyped(i))
	 {
	    X[p] = -1.0;
	    nmiss++;
	    if(nmiss>NMISS)
	    {
	       missingStatus=true;
	       break;
	    }
	 }
	 else if(ped[genotypedSamplePED[p]].markers[i].isHeterozygous())
	 {
	    X[p]=1.0;
	    mean += X[p];
	    n_++;
	 }
	 else if(ped[genotypedSamplePED[p]].markers[i].isHomozygous())
	 {
	    if(ped[genotypedSamplePED[p]].markers[i].hasAllele(2))
	    {
	       X[p]=2.0;
	       mean += X[p];
	       n_++;
	    }
	    else if(ped[genotypedSamplePED[p]].markers[i].hasAllele(1))
	    {
	       X[p]=0.0;
	       n_++;
	    }
	 }
	 else
	 {
	    X[p]=0.0;
	    n_++;
	 }
      }

      if(missingStatus)
	 continue;

      maf = mean/(2.0*n_);

      if(nmiss>NMISS || maf<q || 1.0-maf<q)
	 continue;

      N++;
      mean /= n_;
      var_inv = 1.0/(2.0*maf*(1.0-maf));

      for(int p=0;p<total_n;p++)
      {
	 if(X[p]==-1.0)
	    X[p] =0.0;
	 else
	    X[p] -= mean;
      }
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int r=0;r<total_n;r++)
      {
	 for(int c=r;c<total_n;c++)
	 {
	    allPairs[r][c] += X[r]*X[c]*var_inv;
	 }
      }
   }

   if(X) delete [] X;

   if(N==0)
      error("ERROR! No qualified marker was included in calculating empirical kinsihp matrix.\n");

   printf("  %d markers were included to calculate empirical kinship.\n",N);
   fprintf(log,"  %d markers were included to calculate empirical kinship.\n",N);


   double scale = 1.0/N;

   for(int i=0;i<total_n;i++)
   {
      allPairs[i][i] *= scale;
      for(int j=i+1;j<total_n;j++)
      {
	 allPairs[i][j] *= scale;
	 allPairs[j][i] = allPairs[i][j];
      }
   }

   printf("completed. \n ");
   fprintf(log,"completed. \n ");
   printf("  %d markers were included to calculate empirical kinship.\n",N);
   fprintf(log,"  %d markers were included to calculate empirical kinship.\n",N);

   if(OutputKin::outputKin)
   {
      WriteKinship(ped,allPairs,genotypedSamplePED,true,false,log);
   }
   printf("\n");
   fprintf(log,"\n");
}

void KinshipEmp::SetupVCFAuto(Pedigree & ped, IntArray & genotypedSampleVCF,StringArray & chromosomeVCF,FILE * log)
{
   printf("Estimating kinship matrix ...\n");
   fflush(stdout);
   fprintf(log,"Estimating kinship matrix ...\n");

   int chr_count = chromosomeVCF.Length();
   int N=0;
   int total_n = genotypedSampleVCF.Length();
   int NMISS = total_n * miss;
   allPairs.Dimension(total_n,total_n);

   for(int r=0;r<total_n;r++)
   {
      for(int c=r;c<total_n;c++)
	 allPairs[r][c] = 0.0;
   }
   VcfFileReader reader;
   VcfHeader header;
   VcfRecord record;

   int markerCount = 0;
   //#pragma omp parallel for
   for(int chr_idx=0;chr_idx<chr_count;chr_idx++)
   {
      reader.open(PreMeta::vcfInput,header);
      reader.readVcfIndex();
      reader.setReadSection(chromosomeVCF[chr_idx].c_str());
printf("  processing chromosome %s\n",chromosomeVCF[chr_idx].c_str());
fflush(stdout);
fprintf(log,"  processing chromosome %s\n",chromosomeVCF[chr_idx].c_str());

      float * genotype = new float [total_n];
      String chr;
      int pos;

      if(PreMeta::dosage)
      {
	 double mean = 0.0, var_inv = 0.0, maf=0.0;
	 int n_=0, nmiss=0;
	 while(reader.readRecord(record))
	 {
	    ++markerCount;
	    chr = record.getChromStr();
	    pos = record.get1BasedPosition();

	    if(chr==PreMeta::xLabel)
	       if(pos>=PreMeta::Xstart && pos<=PreMeta::Xend)
		  continue;

	    mean = 0.0;
	    n_=0; nmiss=0;
	    VcfRecordGenotype & genoInfo = record.getGenotypeInfo();
	    for (int p = 0; p < genotypedSampleVCF.Length(); p++)
	    {
	       int s = genotypedSampleVCF[p];
	       //fill in genotype vector for this sample
	       const std::string * geno = genoInfo.getString("DS",s);
	       if(*geno == ".")
	       {
		  genotype[p] = _NAN_;
		  nmiss++;
	       }
	       else
	       {
		  genotype[p] = atof((*geno).c_str());
		  mean += genotype[p];
		  n_++;
	       }
	       if(nmiss>NMISS)
		  break;
	    }
	    if(nmiss>NMISS)
	       continue;
	    maf = mean/(2.0*n_);
	    if( maf<q || 1.0-maf<q)
	       continue;

	    N++;
	    mean /= n_;
	    var_inv = 1.0/(2.0*maf*(1.0-maf));
	    //standardize genotype vector
	    for(int i=0;i<total_n;i++)
	    {
	       if(genotype[i]==_NAN_)
		  //make sure if genotype is missing, then it is not contributing
		  genotype[i]=0.0;
	       else
		  genotype[i] -= mean;
	    }
	    //add up to kinship
#ifdef _OPENMP
#pragma omp parallel for
#endif
	    for(int r=0;r<total_n;r++)
	    {
	       for(int c=r;c<total_n;c++)
	       {
		  allPairs[r][c] += genotype[r]*genotype[c]*var_inv;
	       }
	    }
	 }
	 reader.close();
      }
      else
      {
	 double var_inv=_NAN_,mean = 0.0,maf=0.0;
	 int n_=0,nmiss=0,a=0;
	 bool skipSNP = false,skip=false;
	 while(reader.readRecord(record))
	 {
	    ++markerCount;
	    chr = record.getChromStr();
	    pos = record.get1BasedPosition();

	    if(chr==PreMeta::xLabel)
	       if(pos>=PreMeta::Xstart && pos<=PreMeta::Xend)
		  continue;

	    mean = 0.0; maf=0.0;
	    n_=0; nmiss=0;
	    for(int g=0;g<total_n;g++)
	       genotype[g] = 0.0;

	    skipSNP=false;
	    skip=false;
	    for(int s=0;s<total_n;s++)
	    {
	       skip=false;
	       int i= genotypedSampleVCF[s];
	       int numGTs = record.getNumGTs(i);
	       for(int j = 0; j < numGTs; j++)
	       {
		  a = record.getGT(i,j);
		  if(a==VcfGenotypeSample::MISSING_GT)
		  {
		     nmiss++;
		     if(nmiss>NMISS)
		     {
			skipSNP=true;
			break;
		     }
		     skip =true;
		     genotype[s]=_NAN_;
		     break;
		  }
		  genotype[s] += a;
	       }
	       if(skipSNP)
		  break;
	       if(skip)
		  continue;
	       mean+=genotype[s];
	       n_++;
	    }

	    if(skipSNP)
	       continue;
	    mean /= n_;
	    maf = mean/2.0;

	    if(nmiss > NMISS || maf<q || 1.0-maf<q)
	       continue;

	    N++;
	    var_inv = 1.0/((1.0-maf)*2.0*maf);

	    for(int i=0;i<total_n;i++)
	    {
	       if(genotype[i]==_NAN_)
		  genotype[i]=0.0;
	       else
		  genotype[i] -= mean;
	    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
	    for(int r=0;r<total_n;r++)
	    {
	       for(int c=r;c<total_n;c++)
	       {
		  allPairs[r][c] += genotype[r]*genotype[c]*var_inv;
	       }
	    }
	 }
      }
      if(genotype) delete [] genotype;
      reader.close();
   }
printf("  processed %d markers in total.\n",markerCount);
fprintf(log,"  processed %d markers in total.\n",markerCount);
printf("  %d markers were included to calculate empirical kinship.\n",N);
fprintf(log,"  %d markers were included to calculate empirical kinship.\n",N);

   double scale = 1.0/N;

   for(int i=0;i<total_n;i++)
   {
      allPairs[i][i] *= scale;
      for(int j=i+1;j<total_n;j++)
      {
	 allPairs[i][j] *= scale;
	 allPairs[j][i] = allPairs[i][j];
      }
   }

   if(N==0)
      error("ERROR! No qualified marker was included in calculating empirical kinsihp matrix.\n");

   printf("completed.\n\n");
   fprintf(log,"completed.\n\n");

time_t now;
    time(&now);
    printf("Relationship inference completed at: %s\n", ctime(&now));
    fprintf(log,"Relationship inference completed at: %s\n", ctime(&now));

   if(OutputKin::outputKin)
   {
      WriteKinship(ped,allPairs,genotypedSampleVCF,true,true,log);
   }
   printf("\n");
   fprintf(log,"\n");
}

void KinshipEmp::WriteKinship(Pedigree & ped,Matrix & allPairs, IntArray & genotypedSample, bool AUTO,bool VCF, FILE * log)
{
  
   IFILE file = NULL;
   String filename;

   if(AUTO)
   {
      if(PreMeta::outputFile != "")
      {
	 if(PreMeta::outputFile.Last()=='.' || PreMeta::outputFile.Last()=='/')
	    filename.printf("%sEmpirical.Kinship.gz",(const char *) PreMeta::outputFile);
	 else
	    filename.printf("%s.Empirical.Kinship.gz",(const char *) PreMeta::outputFile);
      }
      else
	 filename.printf("Empirical.Kinship.gz");
   }
   else
   {
      if(PreMeta::outputFile != "")
      {  
	 if(PreMeta::outputFile.Last()=='.' || PreMeta::outputFile.Last()=='/')
	    filename.printf("%sEmpirical.KinshipX.gz",(const char *) PreMeta::outputFile);
	 else
	    filename.printf("%s.Empirical.KinshipX.gz",(const char *) PreMeta::outputFile);
      }  
      else
	 filename.printf("Empirical.KinshipX.gz");
   }
   int n = genotypedSample.Length();

   if(VCF)
   {
      file = ifopen(filename, "wt", InputFile::GZIP);
      VcfFileReader reader;
      VcfHeader header;
      reader.open(PreMeta::vcfInput,header);
      for (int i = 0; i < n; i++)
      {
	 const char * sample = header.getSampleName(genotypedSample[i]);
	 ifprintf(file,"%s ",sample);
      }
      ifprintf(file,"\n");
      reader.close();
      for(int r=0; r<n; r++)
      {
	 for(int c=0; c<=r; c++)
	 {
	    ifprintf(file,"%g ",allPairs[r][c]);
	 }
	 ifprintf(file,"\n");
      }
      ifclose(file);
   }
   else
   {
      file = ifopen(filename, "wt",InputFile::GZIP);
      for(int p=0;p<n;p++)
      {
	 ifprintf(file,"%s ",ped[genotypedSample[p]].pid.c_str());
      }
      ifprintf(file,"\n");

      for(int r=0;r<n;r++)
      {
	 for(int c=0;c<n;c++)
	 {
	    ifprintf(file,"%g ",allPairs[r][c]);
	 }
	 ifprintf(file,"\n");
      }
      ifclose(file);
   }

   if(AUTO)
   {
      printf("  Estimated kinship saved in %s.\n",filename.c_str());
      fprintf(log,"  Estimated kinship saved in %s.\n",filename.c_str());
   }
   else
   {
      printf("  Estiamted kinshipX saved in %s.\n",filename.c_str());
      fprintf(log,"  Estiamted kinshipX saved in %s.\n",filename.c_str());
   }
}

//this function calculates relationship matrix from autosomal markers
void KinshipEmp::SetupPEDX(Pedigree & ped, IntArray & genotypedSamplePED, FILE * log)
{
   printf("Calculating empirical kinshipX matrix ...");
   fflush(stdout);
   fprintf(log,"Calculating empirical kinshipX matrix ...");
   int N=0;
   int total_n = genotypedSamplePED.Length();

   allPairsX.Dimension(total_n, total_n);
   for(int r=0;r<total_n;r++)
   {
      for(int c=r;c<total_n;c++)
	 allPairsX[r][c] = 0.0;
   }

   float * X = new float [total_n];
   for(int i=0;i<total_n;i++)
      X[i] = -1.0;

   String chr;
   int pos;
   int NMISS=total_n*miss;

   for(int i=0;i<ped.markerCount;i++)
   {
      StringArray markerName;
      markerName.AddTokens(ped.markerNames[i],":");

      if(markerName.Length()<=1)
      {
	 printf("WARNING: marker name %s in DAT file is not in CHR:POS format and will be skipped from calculating empirical kinsihp.\n",ped.markerNames[i].c_str());
	 fprintf(log,"WARNING: marker name %s in DAT file is not in CHR:POS format and will be skipped from calculating empirical kinship.\n",ped.markerNames[i].c_str());
	 continue;
      }

      chr = markerName[0];
      pos = atoi(markerName[1].c_str());

      if(chr!=PreMeta::xLabel)
	 continue;

      int alleleCount = CountAlleles(ped, i);
      //if no one is genotyped, mornomorphic, or marker is not biallelic
      if(alleleCount!=2)
	 continue;

      double mean=0.0;
      double variance = 0.0;
      double mean_sqr = 0.0;

      //now fill in the genotypes
      int n_=0,nmiss=0;
      double maf = 0.0;
      bool missingStatus=false;
      //only generate matrix for genotyped individuals
      /*
	 bool maleGenotypedWrong=false;
	 if(pos>=PreMeta::Xstart && pos<=PreMeta::Xend)
	 {
	 for(int p=0;p<total_n;p++)
	 {
      //check if this sample is male and heterozygous then skip this marker
      if(ped[genotypedSamplePED[p]].sex==PreMeta::maleLabel && ped[genotypedSamplePED[p]].markers[i].isHeterozygous())
      {
      maleGenotypedWrong=true;
      break;
      }
      }
      }
      if(maleGenotypedWrong)
      {
      fprintf(log,"marker %s:%d is excluded from calculating reltionship matrix due observation of heterozygous X chromosome in at least one male sample.\n",chr.c_str(),pos);
      continue;
      }
       */

      for(int p=0;p<total_n;p++)
      {
	 if(!ped[genotypedSamplePED[p]].isGenotyped(i))
	 {
	    X[p] = -1.0;
	    nmiss++;
	    if(nmiss>NMISS)
	    {
	       missingStatus=true;
	       break;
	    }
	 }
	 else if(ped[genotypedSamplePED[p]].markers[i].isHeterozygous())
	 {
	    if(pos>=PreMeta::Xstart && pos<=PreMeta::Xend && ped[genotypedSamplePED[p]].sex==PreMeta::maleLabel)
	       X[p]=2.0;
	    else
	       X[p]=1.0;
	    mean += X[p];
	    mean_sqr += X[p]*X[p];
	    n_++;
	 }
	 else if(ped[genotypedSamplePED[p]].markers[i].isHomozygous())
	 {
	    if(ped[genotypedSamplePED[p]].markers[i].hasAllele(2))
	    {
	       X[p]=2.0;
	       mean += X[p];
	       mean_sqr += X[p]*X[p];
	       n_++;
	    }
	    else if(ped[genotypedSamplePED[p]].markers[i].hasAllele(1))
	    {
	       X[p]=0.0;
	       n_++;
	    }
	 }
	 else
	 {
	    X[p]=0.0;
	    n_++;
	 }
      }

      if(missingStatus)
      {
	 fprintf(log,"marker %s:%d is excluded from calculating reltionship matrix due to missing rate > %g.\n",chr.c_str(),pos,miss);
	 continue;
      }
      maf = mean/(2.0*n_);

      if(nmiss>NMISS || maf<q || 1.0-maf<q)
	 continue;

      N++;
      mean /= n_;
      variance = mean_sqr/n_ - mean*mean;

      for(int p=0;p<total_n;p++)
      {
	 if(X[p]==-1.0)
	    X[p] =0.0;
	 else
	    X[p] -= mean;
      }
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int r=0;r<total_n;r++)
      {
	 for(int c=r;c<total_n;c++)
	 {
	    allPairsX[r][c] += X[r]*X[c]/variance;
	 }
      }
   }

   if(X) delete []X;

   if(N==0)
      error("ERROR! No qualified marker was included in calculating empirical kinsihp matrix of X chromosome.\n");

   for(int i=0;i<total_n;i++)
      for(int j=i+1;j<total_n;j++)
	 allPairsX[j][i] = allPairsX[i][j];

   allPairsX.Multiply(1.0/N);

   printf("completed\n ");
   fprintf(log,"completed\n ");

   printf("  %d markers were included to calculate empirical kinship of X chromosome.\n",N);
   fprintf(log,"  %d markers were included to calculate empirical kinship of X chromosome.\n",N);
   if(OutputKin::outputKin)
   {
      WriteKinship(ped,allPairsX,genotypedSamplePED,false,false,log);
   }
   printf("\n");
   fprintf(log,"\n");
}

void KinshipEmp::SetupVCFX(Pedigree & ped, IntArray & genotypedSampleVCF, StringIntHash & samplePEDIDHash, StringIntHash & skippedSNPs,FILE * log)
{
   printf("Calculating empirical kinshipX matrix ...");
   fflush(stdout);
   fprintf(log,"Calculating empirical kinshipX matrix ...");
   VcfFileReader reader;
   VcfHeader header;
   VcfRecord record;
   reader.open(PreMeta::vcfInput,header);
   reader.readVcfIndex();
   reader.setReadSection(PreMeta::xLabel.c_str());

   int N=0;
   int total_n = genotypedSampleVCF.Length();
   allPairsX.Dimension(total_n,total_n);
   for(int r=0;r<total_n;r++)
   {
      for(int c=r;c<total_n;c++)
	 allPairsX[r][c] = 0.0;
   }

   float * genotype = new float [total_n];

   String chr;
   int pos;
   int NMISS = total_n *miss;

   IntArray counts;
   counts.Dimension(2);

   double mean = 0.0,maf=_NAN_,var_inv=_NAN_;
   int n_=0, nmiss=0;
   if(PreMeta::dosage)
   {
      while(reader.readRecord(record))
      {
	 chr = record.getChromStr();
	 pos = record.get1BasedPosition();

	 if(pos<PreMeta::Xstart || pos>PreMeta::Xend)
	    continue;
	 mean = 0.0;
	 n_=0; nmiss=0;
	 VcfRecordGenotype & genoInfo = record.getGenotypeInfo();
	 for (int p = 0; p < total_n; p++)
	 {   
	    int s = genotypedSampleVCF[p];

	    //fill in genotype vector for this sample
	    const std::string * geno = genoInfo.getString("DS",s);
	    if(*geno == ".")
	    {
	       genotype[p] = _NAN_;
	       nmiss++;
	    }
	    else
	    {
	       genotype[p] = atof((*geno).c_str());
	       mean += genotype[p];
	       n_++;
	    }
	    if(nmiss>NMISS)
	       break;
	 }
	 if(nmiss>NMISS)
	    continue;
	 maf = mean/(2.0*n_);
	 if( maf<q || 1.0-maf<q)
	    continue;
	 N++;
	 mean = maf*2.0;
	 var_inv = 1.0/(2.0*maf*(1.0-maf));
	 //standardize genotype vector
	 for(int i=0;i<total_n;i++)
	 {
	    if(genotype[i]==_NAN_)
	       //make sure if genotype is missing, then it is not contributing
	       genotype[i]=0.0;
	    else
	       genotype[i] -= mean;
	 }
	 //add up to kinship
#ifdef _OPENMP
#pragma omp parallel for
#endif
	 for(int r=0;r<total_n;r++)
	 {
	    for(int c=r;c<total_n;c++)
	    {
	       allPairsX[r][c] += genotype[r]*genotype[c]*var_inv;
	    }
	 }
      }
   }
   else
   {
      while(reader.readRecord(record))
      {
	 pos = record.get1BasedPosition();
	 if(pos<PreMeta::Xstart || pos>PreMeta::Xend)
	    continue;
	 chr = record.getChromStr();

	 String SNPname = chr + ":"+ pos;
/*
	    if(skippedSNPs.Integer(SNPname)>0)
	       continue;
*/
	 mean = 0.0;
	 n_=0;
	 nmiss=0;

	 counts.Zero();

	 bool skipSNP=false;
	 for(int s=0;s<genotypedSampleVCF.Length();s++)
	 {
	    bool skip=false;
	    int i= genotypedSampleVCF[s];
	    int numGTs = record.getNumGTs(i);

	    const char * sample = header.getSampleName(i);
	    int ped_idx = samplePEDIDHash.Integer(sample);

	    if(ped[ped_idx].sex==PreMeta::maleLabel 
		  && pos>=PreMeta::Xstart && pos<=PreMeta::Xend) 
	    {
	       int * geno_tmp = new int [numGTs];
	       for(int j = 0; j < numGTs; j++)
	       {
		  //if marker is not biallelic, skip this marker
		  int a = record.getGT(i,j);
		  if(a>1)
		  {
		     skipSNP=true;
		     warnings++;
		     fprintf(log,"Warning: vairant %s is skipped because male has allele that is not 0 or 1.\n",SNPname.c_str());
		     if(warnings<20)
		     {
			printf("Warning: vairant %s is skipped because male has allele that is not 0 or 1.\n",SNPname.c_str());
		     }
		     break;
		  }
		  geno_tmp[j] = a;
	       }

	       if(skipSNP)
		  break;
	       if(numGTs==1)
	       {
		  if(geno_tmp[0]==VcfGenotypeSample::MISSING_GT)
		  {
		     nmiss++;
		     if(nmiss>NMISS)
		     {
			skipSNP=true;
			warnings++;
			fprintf(log,"Warning: variant %s is skipped when calculating kinshipX, because genotype missing rate is >%g.\n",SNPname.c_str(),miss);
			if(warnings<20)
			{
			   printf("Warning: variant %s is skipped when calculating kinshipX, because genotype missing rate is >%g.\n",SNPname.c_str(),miss);
			}
			break;
		     }
		     skip =true;
		     genotype[s]=_NAN_;
		  }
		  else
		  {
		     counts[geno_tmp[0]] += 2;
		     genotype[s] = geno_tmp[0]*2;
		  }
	       }
	       else
	       {
		  if(geno_tmp[0] != geno_tmp[1])
		  {
		     nmiss++;
		     skip =true;
		     warnings++;
		     fprintf(log,"Warning: genotype of sample %s of variant %s is set to be missing because male genotype on non-pseudo-autosomal region of chr X can not be heterozygous.\n",sample,SNPname.c_str());
		     if(warnings<20)
		     {
			printf("Warning: genotype of sample %s of variant %s is set to be missing because male genotype on non-pseudo-autosomal region of chr X can not be heterozygous.\n",sample,SNPname.c_str());
		     }
		     if(nmiss>NMISS)
		     {
			skipSNP=true;
			warnings++;
			fprintf(log,"Warning: variant %s is skipped when calculating kinshipX, because genotype missing rate is >%g.\n",SNPname.c_str(),miss);
			if(warnings<20)
			{
			   printf("Warning: variant %s is skipped when calculating kinshipX, because genotype missing rate is >%g.\n",SNPname.c_str(),miss);
			}
			break;
		     }
		     genotype[s]=_NAN_;
		  }
		  else if(geno_tmp[0]==VcfGenotypeSample::MISSING_GT)
		  {
		     nmiss++;
		     skip =true;
		     if(nmiss>NMISS)
		     {
			skipSNP=true;
			warnings++;
			fprintf(log,"Warning: variant %s is skipped when calculating kinshipX, because genotype missing rate is >%g.\n",SNPname.c_str(),miss);
			if(warnings<20)
			{
			   printf("Warning: variant %s is skipped when calculating kinshipX, because genotype missing rate is >%g.\n",SNPname.c_str(),miss);
			}
			break;
		     }
		     genotype[s]=_NAN_;
		  }
		  else
		  {
		     genotype[s]=0;
		     for(int j=0;j<numGTs;j++)
			genotype[s] += geno_tmp[j];
		     counts[geno_tmp[0]]++;
		     counts[geno_tmp[1]]++;
		  }
	       }
	       if(skipSNP)
		  break;
	       if(skip)
		  continue;
	    }
	    else
	    {
	       for(int j = 0; j < numGTs; j++)
	       {  
		  int a = record.getGT(i,j);
		  if(a==VcfGenotypeSample::MISSING_GT)
		  {  
		     nmiss++;
		     if(nmiss>NMISS)
		     {  
			skipSNP=true;
			warnings++;
			fprintf(log,"Warning: variant %s is skipped when calculating kinshipX, because genotype missing rate is >%g.\n",SNPname.c_str(),miss);
			if(warnings<20)
			{
			   printf("Warning: variant %s is skipped when calculating kinshipX, because genotype missing rate is >%g.\n",SNPname.c_str(),miss);
			}
			break;
		     }  
		     skip =true;
		     genotype[s]=_NAN_;
		     break;
		  }  
		  //if marker is not biallelic, skip this marker
		  if(a>1)
		  {
		     skipSNP=true;
		     warnings++;
		     fprintf(log,"Warning: variant %s is skipped when calculating kinshipX, because marker is not bi-allelic.\n",SNPname.c_str());
		     if(warnings<20)                                 {
			printf("Warning: variant %s is skipped when calculating kinshipX, because marker is not bi-allelic.\n",SNPname.c_str());
		     }
		     break;
		  }
		  genotype[s] += a;
		  counts[a]++;
	       }
	    }
	    if(skipSNP)
	       break;
	    if(skip)
	       continue;
	    mean+=genotype[s];
	    n_++;
	 }

	 if(skipSNP)
	    continue;

	 maf = (double)counts[1]/(counts[0]+counts[1]);

	 if(nmiss > NMISS || maf<q || 1.0-maf<q)
	    continue;

	 N++;
	 mean /= n_;
	 var_inv = 1/(2.0*maf*(1.0-maf));

	 for(int i=0;i<total_n;i++)
	 {
	    if(genotype[i]==_NAN_)
	       genotype[i]=0.0;
	    else
	       genotype[i] -= mean;
	 }
#ifdef _OPENMP
#pragma omp parallel for
#endif
	 for(int r=0;r<total_n;r++)
	 {
	    for(int c=r;c<total_n;c++)
	    {
	       allPairsX[r][c] += genotype[r]*genotype[c]*var_inv;
	    }
	 }
      }
   }

   for(int i=0;i<total_n;i++)
      for(int j=i+1;j<total_n;j++)
	 allPairsX[j][i] = allPairsX[i][j];

   if(genotype) delete [] genotype;

   if(N==0)
      error("ERROR! No qualified marker was included in calculating empirical kinsihp matrix.\n");

   allPairsX.Multiply(1.0/N);

   printf("completed\n ");
   fprintf(log,"completed\n ");

   printf("  %d markers were included to calculate empirical kinship of X chromosome.\n",N);
   fprintf(log,"  %d markers were included to calculate empirical kinship of X chromosome.\n",N);

   if(OutputKin::outputKin)
   {
      WriteKinship(ped,allPairsX,genotypedSampleVCF,false,true,log);
   }
   printf("\n");
   fprintf(log,"\n");
}
