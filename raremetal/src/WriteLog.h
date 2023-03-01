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
#include "Meta.h"
#include "GroupFromAnnotation.h"

void WriteLog(Meta &meta, GroupFromAnnotation &group, FILE *log)
{
    fprintf(log, "The following parameters are in effect:\n\n");
    fprintf(log, "List of Studies:\n");
    fprintf(log, "============================\n");
    fprintf(log, "--summaryFiles [%s]\n", meta.summaryFiles.c_str());
    fprintf(log, "--covFiles [%s]\n", meta.covFiles.c_str());
//   fprintf(log, "Implemented Meta-analysis Methods:\n");
//   fprintf(log,"============================\n");
//   fprintf(log, "--useExact [%s]\n", meta.useExactMetaMethod? "true":"false");
    fprintf(log, "\nGrouping Methods:\n");
    fprintf(log, "============================\n");
    fprintf(log, "--groupFile [%s]\n", group.groupFile.c_str());
    fprintf(log, "--annotatedVcf [%s]\n", group.vcfInput.c_str());
    fprintf(log, "--annotation [%s]\n", group.function.c_str());
    fprintf(log, "--writeVcf [%s]\n", meta.outvcf ? "ON" : "OFF");
    fprintf(log, "\nQC Options:\n");
    fprintf(log, "============================\n");
    fprintf(log, "--hwe [%g]\n", meta.HWE);
    fprintf(log, "--callRate [%g]\n", meta.CALLRATE);
    fprintf(log, "\nAssociation Methods:\n");
    fprintf(log, "============================\n");
    fprintf(log, "--burden [%s]\n", meta.Burden ? "true" : "false");
    fprintf(log, "--MB [%s]\n", meta.MB ? "true" : "false");
    fprintf(log, "--SKAT [%s]\n", meta.SKAT ? "true" : "false");
    fprintf(log, "--VT [%s]\n", meta.VTa ? "true" : "false");
    fprintf(log, "--condition [%s]\n", meta.cond.c_str());
    //fprintf(log,"--permute [%s]\n",meta.VTp?"true":"false");
    fprintf(log, "\nOther Options:\n");
    fprintf(log, "============================\n");
    fprintf(log, "--tabix [%s]\n", meta.tabix ? "true" : "false");
    fprintf(log, "--correctGC [%s]\n", meta.correctGC ? "true" : "false");
    fprintf(log, "--prefix [%s]\n", meta.prefix.c_str());
    fprintf(log, "--maf [%g]\n", meta.MAF_cutoff);
    fprintf(log, "--longOutput [%s]\n", meta.fullResult ? "true" : "false");
    fprintf(log, "--tabulateHits [%s]\n", meta.report ? "true" : "false");
    //fprintf(log,"--founderAF [%s]\n",meta.founderAF?"true":"false");
    fprintf(log, "--dosage [%s]\n", meta.dosage ? "ON" : "OFF");
    fprintf(log, "--hitsCutoff [%g]\n", meta.report_pvalue_cutoff);
    fprintf(log, "--altMAF [%s]\n", meta.altMAF ? "ON" : "OFF");
}

#endif
