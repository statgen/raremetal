////////////////////////////////////////////////////////////////////// 
// WriteLog.h 
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

#ifndef __WRITELOG_H__
#define __WRITELOG_H__

#include <stdio.h>
#include "StringBasics.h"
#include "WritePDF.h"


void WriteLog(String & pedfile,String & datfile,bool noeof,FILE * log)
{
   fprintf(log,"\nRare-Metal-Worker handled all individuals as %s.\n\n",FastFit::unrelated?"unrelated":"related");
   fprintf(log,"The following parameters are in effect:\n\n");
   fprintf(log,"\nInput Files:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--ped [%s]\n",pedfile.c_str());
   fprintf(log,"--dat [%s]\n",datfile.c_str());
   fprintf(log,"--vcf [%s]\n",PreMeta::vcfInput.c_str());
   fprintf(log,"--dosage [%s]\n",PreMeta::dosage ? "true":"false");
   fprintf(log,"--noeof [%s]\n", noeof ? "true":"false");
   fprintf(log,"\nOutput Files:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--prefix [%s]\n",PreMeta::outputFile.c_str());
   fprintf(log,"--LDwindow [%d]\n",PreMeta::window);
   fprintf(log,"--zip [%s]\n", PreMeta::zipOutput ? "true":"false");
   fprintf(log,"--thin [%s]\n",WritePDF::thinPoints ? "true":"false");
   fprintf(log,"--labelHits [%s]\n",GroupFromAnnotation::labelHits ? "true":"false");
   fprintf(log,"\nVC Options:\n");
   fprintf(log,"============================\n");
//   fprintf(log,"--vcShared [%s]\n",AutoFit::fitSharedEnvironment ? "true":"false");
   fprintf(log,"--vcX [%s]\n",AutoFit::fitX?"true":"false");
   fprintf(log,"--separateX [%s]\n",FastFit::separateX?"true":"false");
 //  fprintf(log,"--useCovariates [%s]\n\n",FastFit::makeResiduals || FastFit::useCovariates?"true":"false");
   fprintf(log,"\nTrait Options:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--makeResiduals [%s]\n",FastFit::makeResiduals?"true":"false");
   fprintf(log,"--inverseNormal [%s]\n",FastFit::inverseNormal?"true":"false");
   fprintf(log,"--traitName [%s]\n",FastFit::traitName.c_str());
 fprintf(log,"\nModel Options:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--recessive [%s]\n",PreMeta::recessive?"true":"false");
   fprintf(log,"--dominant [%s]\n",PreMeta::dominant?"true":"false");
   fprintf(log,"\nKinship Source:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--kinPedigree [%s]\n",FastTransform::pedKin?"true":"false");
   fprintf(log,"--kinGeno [%s]\n",FastTransform::empKin?"true":"false");
   fprintf(log,"--kinFile [%s]\n",FastTransform::readInEmp.c_str());
   fprintf(log,"--kinxFile [%s]\n",FastTransform::readInEmpX.c_str());
   fprintf(log,"--kinSave [%s]\n",OutputKin::outputKin?"true":"false");
   fprintf(log,"\nKinship Options:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--kinMaf [%g]\n",KinshipEmp::q);
   fprintf(log,"--kinMiss [%g]\n",KinshipEmp::miss);
   fprintf(log,"\nChromosome X:\n");
   fprintf(log,"============================\n");
   fprintf(log,"xLabel [%s]\n",PreMeta::xLabel.c_str());
   fprintf(log,"xStart [%d]\n",PreMeta::Xstart);
   fprintf(log,"xEnd [%d]\n\n",PreMeta::Xend);
   fprintf(log,"maleLabel [%d]\n",PreMeta::maleLabel);
   fprintf(log,"femaleLabel [%d]\n\n",PreMeta::femaleLabel);
   fprintf(log,"\nOthers:\n");
   fprintf(log,"============================\n");
 fprintf(log,"cpus [%d]\n\n",KinshipEmp::cpus);
}

#endif
