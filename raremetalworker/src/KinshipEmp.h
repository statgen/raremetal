////////////////////////////////////////////////////////////////////// 
// KinshipEmp.h 
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

#ifndef __EMPKINSHIP_H__
#define __EMPKINSHIP_H__

#include "Kinship.h"
#include "StringHash.h"


class KinshipEmp : public Kinship
{
   public:

      KinshipEmp() {};
      ~KinshipEmp() {};

      static double q;
      static double miss;
      static int cpus;

      StringIntHash IDFromEmp,IDFromEmpX;

      Matrix allPairs;
      Matrix allPairsX;
      int warnings;

      void ReadEmpKin(FILE * log);
      void ReadEmpKinX(FILE * log);
      void CleanUpAuto();
      void CleanUpX();

      void SetupPEDAuto(Pedigree & ped, IntArray & genotypedSamplePED, FILE * log);
      void SetupPEDX(Pedigree & ped, IntArray & genotypedSamplePED, FILE * log);
      void SetupVCFAuto(Pedigree & ped, IntArray & genotypedSampleVCF, StringArray & chromosomeVCF,FILE * log);
      void SetupVCFX(Pedigree & ped, IntArray & genotypedSampleVCF, StringIntHash & samplePEDIDHash, StringIntHash & skippedSNPs,FILE * log);
      void SetupEmpKin(Pedigree & ped, IntArray & genotypedSamplePED, IntArray & genotypedSampleVCF, StringIntHash & samplePEDIDHash, StringIntHash & skippedSNPs,StringArray & chromosomeVCF,FILE * log);
      void SetupEmpKinX(Pedigree & ped, IntArray & genotypedSamplePED, IntArray & genotypedSampleVCF, StringIntHash & samplePEDIDHash, StringIntHash & skippedSNPs,FILE * log);
void WriteKinship(Pedigree & ped,Matrix & allPairs,IntArray & genotypedSample,bool AUTO,bool VCF,FILE * log);
};

#endif
