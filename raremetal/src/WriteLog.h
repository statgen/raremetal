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

void WriteLog(FILE * log)
{
   fprintf(log,"The following parameters are in effect:\n\n");
   fprintf(log,"List of Studies:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--summaryFiles [%s]\n",Meta::summaryFiles.c_str());
   fprintf(log,"--covFiles [%s]\n",Meta::covFiles.c_str());
   fprintf(log,"\nGrouping Methods:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--groupFile [%s]\n",GroupFromAnnotation::groupFile.c_str());
   fprintf(log,"--annotatedVcf [%s]\n",GroupFromAnnotation::vcfInput.c_str());
   fprintf(log,"--annotation [%s]\n",GroupFromAnnotation::function.c_str());
   fprintf(log,"--writeVcf [%s]\n",Meta::outvcf?"ON":"OFF");
   fprintf(log,"\nQC Options:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--hwe [%g]\n",Meta::HWE);
   fprintf(log,"--callRate [%g]\n",Meta::CALLRATE);
   fprintf(log,"\nAssociation Methods:\n");
   fprintf(log,"============================\n");
   fprintf(log,"--burden [%s]\n",Meta::Burden?"true":"false");
   fprintf(log,"--MB [%s]\n",Meta::MB?"true":"false");
   fprintf(log,"--SKAT [%s]\n",Meta::SKAT?"true":"false");
   fprintf(log,"--VT [%s]\n",Meta::VTa?"true":"false");
fprintf(log,"--condition [%s]\n",Meta::cond.c_str());
   //fprintf(log,"--permute [%s]\n",Meta::VTp?"true":"false");
   fprintf(log,"\nOther Options:\n");
   fprintf(log,"============================\n");
fprintf(log,"--tabix [%s]\n",Meta::tabix?"true":"false");
fprintf(log,"--correctGC [%s]\n",Meta::correctGC?"true":"false");
   fprintf(log,"--prefix [%s]\n",Meta::prefix.c_str());
   fprintf(log,"--maf [%g]\n",Meta::MAF_cutoff);
   fprintf(log,"--longOutput [%s]\n",Meta::fullResult?"true":"false");
   fprintf(log,"--tabulateHits [%s]\n",Meta::report?"true":"false");
   //fprintf(log,"--founderAF [%s]\n",Meta::founderAF?"true":"false");
   fprintf(log,"--dosage [%s]\n", Meta::dosage ? "ON" : "OFF");
   fprintf(log,"--hitsCutoff [%g]\n",Meta::report_pvalue_cutoff);
   fprintf(log,"--altMAF [%s]\n", Meta::altMAF ? "ON" : "OFF");
}

#endif
