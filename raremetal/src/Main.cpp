////////////////////////////////////////////////////////////////////// 
// Main.cpp 
// (c) 2013-2014 Shuang Feng, Dajiang Liu, Goncalo Abecasis
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

int main(int argc, char ** argv)
{
   printf("\nRAREMETAL %s -- A Tool for Rare Variants Meta-Analyses for Quantitative Traits\n"
	 "          (c) 2012-2013 Shuang Feng, Dajiang Liu, Goncalo Abecasis\n\n",VERSION);
printf("\nPlease go to \"http://genome.sph.umich.edu/wiki/RAREMETAL\" for the newest version.\n");
  PhoneHome::allThinning = 100;
    time_t initialTime = time(0);

   BEGIN_LONG_PARAMETERS(additional)
      LONG_PARAMETER_GROUP("List of Studies")
      LONG_STRINGPARAMETER("summaryFiles", &Meta::summaryFiles)
      LONG_STRINGPARAMETER("covFiles", &Meta::covFiles)
      LONG_PARAMETER_GROUP("Grouping Methods")
      LONG_STRINGPARAMETER("groupFile", &GroupFromAnnotation::groupFile)
      LONG_STRINGPARAMETER("annotatedVcf", &GroupFromAnnotation::vcfInput)
      LONG_STRINGPARAMETER("annotation", &GroupFromAnnotation::function)
      LONG_PARAMETER("writeVcf", &Meta::outvcf)
      LONG_PARAMETER_GROUP("QC Options")
      LONG_DOUBLEPARAMETER("hwe", &Meta::HWE)
      LONG_DOUBLEPARAMETER("callRate", &Meta::CALLRATE)
      LONG_PARAMETER_GROUP("Association Methods")
      LONG_PARAMETER("burden", &Meta::Burden)
      LONG_PARAMETER("MB", &Meta::MB)
      LONG_PARAMETER("SKAT", &Meta::SKAT)
      //LONG_PARAMETER("SKATO", &Meta::SKATO)
      LONG_PARAMETER("VT", &Meta::VTa)
      LONG_STRINGPARAMETER("condition", &Meta::cond)
      //LONG_PARAMETER("permute", &Meta::VTp)
      LONG_PARAMETER_GROUP("Other Options")
      //LONG_PARAMETER("tabix", &Meta::tabix)
      LONG_PARAMETER("labelHits", &GroupFromAnnotation::labelHits)
      LONG_STRINGPARAMETER("geneMap", &GroupFromAnnotation::mapFile)
      LONG_PARAMETER("correctGC", &Meta::correctGC)
      LONG_STRINGPARAMETER("prefix", &Meta::prefix)
//LONG_STRINGPARAMETER("mapFile", &GroupFromAnnotation::mapFile)
      LONG_DOUBLEPARAMETER("maf", &Meta::MAF_cutoff)
      LONG_PARAMETER("longOutput", &Meta::fullResult)
      //LONG_PARAMETER("founderAF", &Meta::founderAF)
      LONG_PARAMETER("tabulateHits", &Meta::report)
      LONG_PARAMETER("dosage", &Meta::dosage)
      LONG_DOUBLEPARAMETER("hitsCutoff", &Meta::report_pvalue_cutoff)

      LONG_PHONEHOME(VERSION)
      END_LONG_PARAMETERS();


   ParameterList pl;
   pl.Add(new LongParameters("Options:", additional));

   pl.Read(argc, argv);
   pl.Status();

   PhoneHome::checkVersion("raremetalworker",VERSION);

   FILE * logFile;
   String filename;
   if(Meta::prefix=="")
      filename = "raremetal.log";
   else if(Meta::prefix.Last()=='.' || Meta::prefix.Last()=='/')
      filename = Meta::prefix + "raremetal.log";
   else
      filename = Meta::prefix + ".raremetal.log";

   logFile = freopen(filename,"wt",stderr);
   WriteLog(logFile);

    time_t now;
    time(&now);
    printf("Analysis started at: %s\n", ctime(&now));
    fprintf(logFile,"Analysis started at: %s\n", ctime(&now));

   try {
      if(Meta::summaryFiles=="")
      {
          error("--summaryFiles can not be empty.\n");
      }
      GroupFromAnnotation group;

      if(GroupFromAnnotation::groupFile!="" 
		&& GroupFromAnnotation::vcfInput !="")
      {
	 printf("Warning: you have entered both groupfile and annotated VCF file. Groups will be read from the group file only.\n");
	 GroupFromAnnotation::vcfInput="";
      }

      String path;
      path = argv[0];
      Meta meta;
      meta.Prepare();
      group.Run(path,logFile);
      meta.PoolSummaryStat(group,logFile);
      meta.Run(group,logFile);

   time(&now);
       printf("\nAnalysis ends at: %s\n", ctime(&now));
       fprintf(logFile,"\nAnalysis ends at: %s\n", ctime(&now));
       time_t endTime = time(0);
       int timeSec = (int) (endTime - initialTime);
       double timeHour = timeSec/3600.0;
       printf("Analysis took %d seconds (~ %.1f hours).\n",timeSec, timeHour);
       fprintf(logFile,"Analysis took %d seconds (~ %.1f hours).\n",timeSec, timeHour);

      fclose(logFile);
   }
   catch(std::exception& e)
   {
      printf("Exiting, exception thrown: %s\n\n",e.what());
      fprintf(logFile,"Exiting, exception thrown: %s\n\n",e.what());
      PhoneHome::completionStatus("Exception");
      return(-1);
   }
   PhoneHome::completionStatus("SUCCESS");
   return 0;
}
