//////////////////////////////////////////////////////////////////////
// Main.cpp
// (c) 2013-2017 Sai Chen, Shuang Feng, Dajiang Liu, Goncalo Abecasis
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile RareMetalWorker.
//
// All computer programs have bugs. Use this file at your own risk.
//
// February 13, 2013
//

#include "Parameters.h"
#include "math.h"
#include <iostream>
#include <stdio.h>
#include <time.h>
#include "Error.h"
#include "Meta.h"
#include "GroupFromAnnotation.h"
#include "WriteLog.h"
#include "unistd.h"

int main(int argc, char *argv[])
{
    printf("\nRAREMETAL %s -- A Tool for Rare Variants Meta-Analyses for Quantitative Traits\n"
                   "          (c) 2012-2017 Shuang Feng, Dajiang Liu, Sai Chen, Goncalo Abecasis\n\n", VERSION);
    printf("\nPlease go to \"http://genome.sph.umich.edu/wiki/RAREMETAL\" for the newest version.\n");
    PhoneHome::allThinning = 100;
    time_t initialTime = time(0);

    Meta meta;
    GroupFromAnnotation group;

    BEGIN_LONG_PARAMETERS(additional)
        LONG_PARAMETER_GROUP("List of Studies")
            LONG_STRINGPARAMETER("summaryFiles", &meta.summaryFiles)
            LONG_STRINGPARAMETER("covFiles", &meta.covFiles)
        LONG_PARAMETER_GROUP("Grouping Methods")
            LONG_STRINGPARAMETER("groupFile", &group.groupFile)
            LONG_STRINGPARAMETER("annotatedVcf", &group.vcfInput)
            LONG_STRINGPARAMETER("annotation", &group.function)
            LONG_PARAMETER("writeVcf", &meta.outvcf)
        LONG_PARAMETER_GROUP("QC Options")
            LONG_DOUBLEPARAMETER("hwe", &meta.HWE)
            LONG_DOUBLEPARAMETER("callRate", &meta.CALLRATE)
        LONG_PARAMETER_GROUP("Association Methods")
            LONG_PARAMETER("burden", &meta.Burden)
            LONG_PARAMETER("MB", &meta.MB)
            LONG_PARAMETER("MAB", &meta.MAB)
            LONG_PARAMETER("BBeta", &meta.BBeta)
            LONG_PARAMETER("SKAT", &meta.SKAT)
            //	LONG_PARAMETER("SKATO", &meta.SKATO)
            LONG_PARAMETER("VT", &meta.VTa)
            LONG_STRINGPARAMETER("condition", &meta.cond)
            //LONG_PARAMETER("permute", &meta.VTp)
        LONG_PARAMETER_GROUP("Other Options")
            //LONG_PARAMETER("tabix", &meta.tabix)
            LONG_PARAMETER("labelHits", &group.labelHits)
            LONG_STRINGPARAMETER("geneMap", &group.mapFile)
            LONG_PARAMETER("correctGC", &meta.correctGC)
            LONG_STRINGPARAMETER("prefix", &meta.prefix)
            //LONG_STRINGPARAMETER("mapFile", &group.mapFile)
            LONG_DOUBLEPARAMETER("maf", &meta.MAF_cutoff)
            LONG_PARAMETER("longOutput", &meta.fullResult)
            //LONG_PARAMETER("founderAF", &meta.founderAF)
            LONG_PARAMETER("tabulateHits", &meta.report)
            LONG_PARAMETER("dosage", &meta.dosage)
            LONG_DOUBLEPARAMETER("hitsCutoff", &meta.report_pvalue_cutoff)
            LONG_PARAMETER("altMAF", &meta.altMAF)
            LONG_STRINGPARAMETER("range", &meta.Region)
            LONG_PARAMETER("useExact", &meta.useExactMetaMethod)
            LONG_PARAMETER("normPop", &meta.normPop)
            LONG_STRINGPARAMETER("popFile", &meta.popfile_name)
            // related binary trait: Dajiang's method
            LONG_PARAMETER("relateBinary", &meta.relateBinary)
            LONG_PARAMETER("debug", &meta.debug)
            LONG_PARAMETER("matchOnly", &meta.matchOnly)
            LONG_PARAMETER("matchByAbs", &meta.matchByAbs)
            LONG_DOUBLEPARAMETER("matchDist", &meta.matchDist)
            LONG_DOUBLEPARAMETER("minMatchMAF", &meta.minMatchMAF)
            LONG_DOUBLEPARAMETER("maxMatchMAF", &meta.maxMatchMAF)
            //	LONG_PARAMETER("noAdjustUnmatch",&meta.noAdjustUnmatch)
            LONG_STRINGPARAMETER("dosageOptionFile", &meta.dosageOptionFile)
            LONG_PARAMETER("sumCaseAC", &meta.sumCaseAC)
            LONG_PARAMETER("heterogeneity", &meta.bHeterogeneity)
            LONG_PARAMETER("logPvalue", &meta.logP)
        LONG_PHONEHOME(VERSION)
    END_LONG_PARAMETERS();


    ParameterList pl;
    pl.Add(new LongParameters("Options:", additional));

    pl.Read(argc, argv);
    pl.Status();

    PhoneHome::checkVersion("raremetalworker", VERSION);

    if (meta.relateBinary) {
      // RPW 2020-03-23: code related to relateBinary has been commented out elsewhere in the program,
      // thus this option should not be allowed
      error("--relateBinary is no longer supported\n");
    }

    if (meta.normPop && !meta.useExactMetaMethod) {
      error("population correction only works for exact method. Please specify --useExact!\n");
    }

    if (meta.useExactMetaMethod) {
      printf("\n\nWARNING: This method only works for unrelated samples! Plus if you have covariates, please make sure --makeResiduals is specified when running Raremetalworker!\n\n");
    }

    if (meta.normPop && meta.popfile_name == "") {
      error("Must provide --popFile when --normPop is specified");
    }

    meta.setLogFile();
    WriteLog(meta, group, meta.log);

    time_t now;
    time(&now);
    printf("Analysis started at: %s\n", ctime(&now));
    fprintf(meta.log, "Analysis started at: %s\n", ctime(&now));

    if (meta.summaryFiles == "") {
      error("--summaryFiles can not be empty.\n");
    }

    if (group.groupFile != "" && group.vcfInput != "") {
      printf("Warning: you have entered both groupfile and annotated VCF file. Groups will be read from the group file only.\n");
      group.vcfInput = "";
    }

    try {
      String path;
      path = argv[0];
      meta.Prepare();
      group.Run(path, meta.log);
      meta.PoolSummaryStat(group);
      if (!meta.skipOutput) {
        meta.WriteSingleVariantResults(group);
      }
      meta.Run(group);

      time(&now);
      printf("\nAnalysis ends at: %s\n", ctime(&now));
      fprintf(meta.log, "\nAnalysis ends at: %s\n", ctime(&now));
      time_t endTime = time(0);
      int timeSec = (int) (endTime - initialTime);
      double timeHour = timeSec / 3600.0;
      printf("Analysis took %d seconds (~ %.1f hours).\n", timeSec, timeHour);
      fprintf(meta.log, "Analysis took %d seconds (~ %.1f hours).\n", timeSec, timeHour);

      fclose(meta.log);
    }
    catch (std::exception &e) {
      printf("Exiting, exception thrown: %s\n\n", e.what());
      fprintf(meta.log, "Exiting, exception thrown: %s\n\n", e.what());
      PhoneHome::completionStatus("Exception");
      return (-1);
    }
    PhoneHome::completionStatus("SUCCESS");
    return 0;
}
