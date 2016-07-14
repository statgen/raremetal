////////////////////////////////////////////////////////////////////// 
// Main.cpp 
// (c) 2012-2016 Shuang Feng, Dajiang Liu, Goncalo Abecasis
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

#include "Parameters.h"
#include "Pedigree.h"
#include "math.h"
#include <iostream>
#include "FastFit.h"
#include "TraitTransformations.h"
#include "TransformResiduals.h"
#include <stdio.h>
#include <ctime>
#include "Error.h"
#include "AutoFit.h"
#include "PreMeta.h"
#include "OutputKin.h"
#include "WriteLog.h"
#include "BgzfFileType.h"
//#include <lapacke.h>
#include "DataQC.h"
#include "AutoFit2.h"
#include "InputFile.h"
#include "WritePDF.h"

int main(int argc, char ** argv)
{
	time_t initialTime = time(0);
	PhoneHome::allThinning = 100;
	bool noeof = false;

	printf("\nRAREMETALWORKER %s -- A Forerunner of RareMetal\n"
		"          (c) 2012-2016 Shuang Feng, Dajiang Liu, Sai Chen, Goncalo Abecasis\n\n", VERSION);

	printf("\nIf you feel the software is useful, please cite:\n  RAREMETAL: fast and powerful meta-analysis for rare variants; Shuang Feng, Dajiang Liu, Xiaowei Zhan, Mary Kate Wing, Goncalo R. Abecasis; Bioinformatics 2014\n");

	printf("\nPlease go to \"http://genome.sph.umich.edu/wiki/RAREMETAL_DOWNLOAD_%26_BUILD#Where_to_Download\" for the latest version.\n");

	String pedfile, datfile;
	double tol = 0.0001;

	BEGIN_LONG_PARAMETERS(additional)
	LONG_PARAMETER_GROUP("Input Files")
	LONG_STRINGPARAMETER("ped", &pedfile)
	LONG_STRINGPARAMETER("dat", &datfile)
	LONG_STRINGPARAMETER("vcf", &PreMeta::vcfInput)
	LONG_PARAMETER("dosage", &PreMeta::dosage)
	LONG_STRINGPARAMETER("flagDosage", &PreMeta::dosageFlag)
	LONG_PARAMETER("noeof", &noeof)
	LONG_PARAMETER_GROUP("Output Files")
	LONG_STRINGPARAMETER("prefix", &PreMeta::outputFile)
	LONG_INTPARAMETER("LDwindow", &PreMeta::window)
	LONG_PARAMETER("zip", &PreMeta::zipOutput)
	LONG_PARAMETER("thin", &WritePDF::thinPoints)
	LONG_PARAMETER("labelHits", &GroupFromAnnotation::labelHits)
	LONG_PARAMETER_GROUP("VC Options")
       //LONG_PARAMETER("vcShared", &AutoFit::fitSharedEnvironment)
	LONG_PARAMETER("vcX", &AutoFit::fitX)
	LONG_PARAMETER("separateX", &FastFit::separateX)
       //LONG_PARAMETER("useCovariates", &FastFit::useCovariates)
	LONG_PARAMETER_GROUP("Trait Options")
	LONG_PARAMETER("makeResiduals", &FastFit::makeResiduals)
	LONG_PARAMETER("inverseNormal", &FastFit::inverseNormal)
	LONG_STRINGPARAMETER("traitName", &FastFit::traitName)
	LONG_PARAMETER_GROUP("Model Options")
	LONG_PARAMETER("recessive", &PreMeta::recessive)
	LONG_PARAMETER("dominant", &PreMeta::dominant)
       //LONG_PARAMETER("additive", &PreMeta::additive)
	LONG_PARAMETER_GROUP("Kinship Source")
	LONG_PARAMETER("kinPedigree", &FastTransform::pedKin)
	LONG_PARAMETER("kinGeno", &FastTransform::empKin)
	LONG_STRINGPARAMETER("kinFile", &FastTransform::readInEmp)
	LONG_STRINGPARAMETER("kinxFile", &FastTransform::readInEmpX)
	LONG_PARAMETER("kinSave", &OutputKin::outputKin)
	LONG_PARAMETER_GROUP("Kinship Options")
	LONG_DOUBLEPARAMETER("kinMaf", &KinshipEmp::q)
	LONG_DOUBLEPARAMETER("kinMiss", &KinshipEmp::miss)
	LONG_PARAMETER_GROUP("Chromosome X")
	LONG_STRINGPARAMETER("xLabel", &PreMeta::xLabel)
	LONG_INTPARAMETER("xStart", &PreMeta::Xstart)
	LONG_INTPARAMETER("xEnd", &PreMeta::Xend)
	LONG_INTPARAMETER("maleLabel", &PreMeta::maleLabel)
	LONG_INTPARAMETER("femaleLabel", &PreMeta::femaleLabel)
	LONG_PARAMETER_GROUP("others")
	LONG_INTPARAMETER("cpu", &KinshipEmp::cpus)
	LONG_PARAMETER("kinOnly", &OutputKin::kinOnly)
	LONG_STRINGPARAMETER("geneMap",&GroupFromAnnotation::mapFile)
	LONG_PARAMETER("mergedVCFID",&FastTransform::mergedVCFID)
	LONG_PHONEHOME(VERSION)
	BEGIN_LEGACY_PARAMETERS()
	LONG_PARAMETER("useCovariates", &FastFit::useCovariates)
	LONG_STRINGPARAMETER("range", &PreMeta::Region)
    LONG_STRINGPARAMETER("variantList", &PreMeta::varListName) // list of variants to calculate score and cov

    END_LONG_PARAMETERS();

    ParameterList pl;
    pl.Add(new LongParameters("Options:", additional));
    pl.Read(argc, argv);

    PhoneHome::checkVersion("raremetalworker",VERSION);

    if(noeof)
    {
       // Set that the eof block is not required.
    	BgzfFileType::setRequireEofBlock(false);
    }


    String logFile;
    if(PreMeta::outputFile != "")
    {
    	if(PreMeta::outputFile.Last()=='.' || PreMeta::outputFile.Last()=='/')
    		logFile = PreMeta::outputFile + "singlevar.log";
    	else
    		logFile = PreMeta::outputFile + ".singlevar.log";
    }
    else
    	logFile = "singlevar.log";

    FILE * log;
    pl.Status();

    log = freopen (logFile,"wt",stderr);

    WriteLog(pedfile,datfile,noeof,log);
    time_t now;
    time(&now);
    printf("Analysis started at: %s\n", ctime(&now));
    fprintf(log,"Analysis started at: %s\n", ctime(&now));

    try {
    	Pedigree ped;
    	clock_t start;
    	start = clock();
    	printf("Loading PED/DAT files ...\n");
    	fprintf(log,"Loading input files ...\n");
    	printf("  Loading DAT files ...");
    	fprintf(log,"  Loading DAT files ...");
    	ped.Prepare(datfile);
    	printf("done.\n  Loading PED files ...");
    	fprintf(log,"done.\n  Loading PED files ...");
    	ped.Load(pedfile);
    	printf("done.\n\n");
    	fprintf(log,"done.\n\n");

    	printf("Loading data used %.1f seconds.\n\n",( std::clock() - start ) / (double)CLOCKS_PER_SEC);
    	SanityCheck checkData;
    	checkData.Check(ped,log);

    	FastFit engine;
    	FastTransform trans;
    	KinshipEmp kin_emp;

       //check if the sample IDs from PED and VCF files are consistent
    	trans.ScreenSampleID(ped,FastFit::useCovariates);


    	if(OutputKin::kinOnly)
    	{
    		OutputKin::outputKin =true;
    		start = clock();
    		kin_emp.SetupEmpKin(ped,trans.genotypedSamplePED,trans.genotypedSampleVCF,trans.samplePEDIDHash,checkData.skippedSNPs,checkData.chromosomeVCF,log);
    		if(AutoFit::fitX || FastFit::separateX)
    			kin_emp.SetupEmpKinX(ped,trans.genotypedSamplePED,trans.genotypedSampleVCF,trans.samplePEDIDHash,checkData.skippedSNPs,log);
    		return 0;
    	}

       //setup empirical kinship if applicable
    	if(FastTransform::empKin)
    	{
    		start = clock();
    		kin_emp.SetupEmpKin(ped,trans.genotypedSamplePED,trans.genotypedSampleVCF,trans.samplePEDIDHash,checkData.skippedSNPs,checkData.chromosomeVCF,log);
	  //printf("Handling empirical kinship matrix used %.1f minutes.\n\n",( std::clock() - start ) / (double)CLOCKS_PER_SEC/60.0);
	  //fprintf(log,"Handling empirical kinship matrix used %.1f minutes.\n\n",( std::clock() - start ) / (double)CLOCKS_PER_SEC/60.0);
    	}

    	StringArray traits;

    	if(!FastFit::traitName.IsEmpty())
    	{
    		FILE * trFile = fopen(FastFit::traitName,"r");
    		if(trFile!=NULL)
    		{
    			String buffer;
    			while (!feof(trFile))
    			{
    				buffer.ReadLine(trFile);
    				traits.AddTokens(buffer, SEPARATORS);
    			}
    			fclose(trFile);
    		}
    		else
    		{
    			traits.AddTokens(FastFit::traitName,"/");
    		}
    	}
    	else 
    	{
    		traits = ped.traitNames;
    	}

    	for(int i=0;i<traits.Length();i++)
    	{
    		if(i==traits.Length()-1 && !AutoFit::fitX)
    			FastFit::CleanKin=true;

    		FastFit::traitName = traits[i];
    		printf("Analyzing trait \"%s\" ...\n",traits[i].c_str());
    		fprintf(log,"Analyzing trait \"%s\" ...\n",traits[i].c_str());
    		fflush(stdout);

    		IFILE SCOREoutput=NULL;
    		IFILE SCOREoutput_rec=NULL;
    		IFILE SCOREoutput_dom=NULL;
    		IFILE SCOREcov=NULL;
    		IFILE SCOREcov_rec=NULL;
    		IFILE SCOREcov_dom=NULL;
    		String filename,filename_COV,filenamePlots;
    		String filename_rec,filename_COV_rec,filenamePlots_rec;
    		String filename_dom,filename_COV_dom,filenamePlots_dom;

    		if(PreMeta::outputFile != "")
    		{
    			if(PreMeta::outputFile.Last()=='.' || PreMeta::outputFile.Last()=='/')
    			{
    				filename = PreMeta::outputFile + traits[i] + ".singlevar.score";
    				filename_COV = PreMeta::outputFile + traits[i] + ".singlevar.cov";
    				filenamePlots = PreMeta::outputFile + traits[i] + ".plots.pdf";
    				if(PreMeta::recessive)
    				{
    					filename_rec = PreMeta::outputFile + traits[i] + ".recessive.singlevar.score";
    					filename_COV_rec = PreMeta::outputFile + traits[i] + ".recessive.singlevar.cov";
    					filenamePlots_rec = PreMeta::outputFile + traits[i] + ".recessive.plots.pdf";
    				}
    				if(PreMeta::dominant)
    				{
    					filename_dom = PreMeta::outputFile + traits[i] + ".dominant.singlevar.score";
    					filename_COV_dom = PreMeta::outputFile + traits[i] + ".dominant.singlevar.cov";
    					filenamePlots_dom = PreMeta::outputFile + traits[i] + ".dominant.plots.pdf";
    				}
    			}
    			else
    			{
    				filename = PreMeta::outputFile + "." + traits[i] + ".singlevar.score";
    				filename_COV = PreMeta::outputFile + "." + traits[i] + ".singlevar.cov";
    				filenamePlots= PreMeta::outputFile + "." + traits[i] + ".plots.pdf";
    				if(PreMeta::recessive)
    				{
    					filename_rec = PreMeta::outputFile + "." + traits[i] + ".recessive.singlevar.score";
    					filename_COV_rec = PreMeta::outputFile + "." + traits[i] + ".recessive.singlevar.cov";
    					filenamePlots_rec = PreMeta::outputFile + "." + traits[i] + ".recessive.plots.pdf";
    				}
    				if(PreMeta::dominant)
    				{
    					filename_dom = PreMeta::outputFile + "." + traits[i] + ".dominant.singlevar.score";
    					filename_COV_dom = PreMeta::outputFile + "." + traits[i] + ".dominant.singlevar.cov";
    					filenamePlots_dom = PreMeta::outputFile + "." + traits[i] + ".dominant.plots.pdf";
    				}
    			}
    		}
    		else
    		{
    			filename_COV = traits[i] + ".singlevar.cov";
    			filename = traits[i] + ".singlevar.score";
    			filenamePlots = traits[i] + ".plots.pdf";
    			if(PreMeta::recessive)
    			{
    				filename_COV_rec = traits[i] + ".recessive.singlevar.cov";
    				filename_rec = traits[i] + ".recessive.singlevar.score";
    				filenamePlots_rec = traits[i] + ".recessive.plots.pdf";
    			}
    			if(PreMeta::dominant)
    			{
    				filename_COV_dom = traits[i] + ".dominant.singlevar.cov";
    				filename_dom = traits[i] + ".dominant.singlevar.score";
    				filenamePlots_dom = traits[i] + ".dominant.plots.pdf";
    			}
    		}

    		if(PreMeta::zipOutput)
    		{
    			filename += ".txt.gz";
    			filename_COV += ".txt.gz";
    			filename_rec += ".txt.gz";
    			filename_COV_rec += ".txt.gz";
    			filename_dom += ".txt.gz";
    			filename_COV_dom += ".txt.gz";
    			SCOREoutput = ifopen(filename,"w",InputFile::BGZF);
    			SCOREcov = ifopen(filename_COV,"w",InputFile::BGZF);
                if (SCOREoutput==NULL || SCOREcov==NULL)
                    error("cannot open score output or cov file\n");
    			if(PreMeta::recessive)
    			{
    				SCOREoutput_rec = ifopen(filename_rec,"w",InputFile::BGZF);
    				SCOREcov_rec = ifopen(filename_COV_rec,"w",InputFile::BGZF);
                    if (SCOREoutput_rec==NULL || SCOREcov_rec==NULL)
                        error("cannot open recessive score output or cov file\n");
    			}
    			if(PreMeta::dominant)
    			{
    				SCOREoutput_dom = ifopen(filename_dom,"w",InputFile::BGZF);
    				SCOREcov_dom = ifopen(filename_COV_dom,"w",InputFile::BGZF);
                    if (SCOREoutput_dom==NULL || SCOREcov_dom==NULL)
                        error("cannot open dominant score output or cov file\n");
    			}
    		}
    		else
    		{
    			filename += ".txt";
    			filename_COV += ".txt";
    			filename_rec += ".txt";
    			filename_COV_rec += ".txt";
    			filename_dom += ".txt";
    			filename_COV_dom += ".txt";
    			SCOREoutput = ifopen(filename,"w",InputFile::UNCOMPRESSED);
    			SCOREcov = ifopen(filename_COV,"w",InputFile::UNCOMPRESSED);
                if (SCOREoutput==NULL || SCOREcov==NULL)
                    error("cannot open score output or cov file\n");
    			if(PreMeta::recessive)
    			{
    				SCOREoutput_rec = ifopen(filename_rec,"w",InputFile::UNCOMPRESSED);
    				SCOREcov_rec = ifopen(filename_COV_rec,"w",InputFile::UNCOMPRESSED);
                    if (SCOREoutput_rec==NULL || SCOREcov_rec==NULL)
                        error("cannot open recessive score output or cov file\n");
    			}
    			if(PreMeta::dominant)
    			{
    				SCOREoutput_dom = ifopen(filename_dom,"w",InputFile::UNCOMPRESSED);
    				SCOREcov_dom = ifopen(filename_COV_dom,"w",InputFile::UNCOMPRESSED);
                    if (SCOREoutput_dom==NULL || SCOREcov_dom==NULL)
                        error("cannot open dominant score output or cov file\n");
    			}
    		}

    		start = clock();
    		engine.FitModels(SCOREoutput,SCOREoutput_rec,SCOREoutput_dom,ped,tol,trans,kin_emp,log);
    		printf("  Fitting model used %.1f minutes.\n",( std::clock() - start ) / (double)CLOCKS_PER_SEC/60.0);
    		fprintf(log,"  Fitting model used %.1f minutes.\n",( std::clock() - start ) / (double)CLOCKS_PER_SEC/60.0);

	  //printf("\n  Generating association statistics and LD matrices.\n");
	  //fprintf(log,"\n  Generating association statistics and LD matrices.\n");
    		fflush(stdout);
    		String path;
    		path = argv[0];
    		GroupFromAnnotation group;
    		group.Run(path);
    		PreMeta preMeta;
    		start = clock();
    		preMeta.Run(SCOREoutput,SCOREoutput_rec,SCOREoutput_dom,SCOREcov,SCOREcov_rec,SCOREcov_dom,ped,trans,engine,group,log,checkData,kin_emp);

	  //printf("\n  Generating association statistics and LD matrices used %.1f minutes.\n",( std::clock() - start ) / (double)CLOCKS_PER_SEC/60.0);
    		fprintf(log,"  Generating association statistics and LD matrices used %.1f minutes.\n",( std::clock() - start ) / (double)CLOCKS_PER_SEC/60.0);

    		ifclose(SCOREoutput);
    		ifclose(SCOREcov);

    		if(PreMeta::recessive)
    		{
    			ifclose(SCOREoutput_rec);
    			ifclose(SCOREcov_rec);
    		}
    		if(PreMeta::dominant)
    		{
    			ifclose(SCOREoutput_dom);
    			ifclose(SCOREcov_dom);
    		}
    		String scoreFile,covFile;
    		printf("Single variant statistics stored in \n  %s\n",filename.c_str());
    		fprintf(log,"  Single variant statistics stored in \n  %s\n",filename.c_str());

    		if(PreMeta::recessive)
    		{
    			printf("  %s\n",filename_rec.c_str());
    			fprintf(log,"  %s\n",filename_rec.c_str());
    		}

    		if(PreMeta::dominant)
    		{
    			printf("  %s", filename_dom.c_str());
    			fprintf(log,"  %s", filename_dom.c_str());
    		}

    		printf("\n");
    		fprintf(log,"\n");

    		printf("Variance-covariance matrix stored in \n  %s\n",filename_COV.c_str());
    		fprintf(log,"Variance-covariance matrix stored in \n  %s\n",filename_COV.c_str());

    		if(PreMeta::recessive)
    		{
    			printf("  %s\n",filename_COV_rec.c_str());
    			fprintf(log,"  %s\n",filename_COV_rec.c_str());
    		}

    		if(PreMeta::dominant)
    		{
    			printf("  %s", filename_COV_dom.c_str());
    			fprintf(log,"  %s", filename_COV_dom.c_str());
    		}

    		printf("\n");
    		fprintf(log,"\n");

    		printf("QQ and Manhattan plots saved in \n  %s\n",filenamePlots.c_str());
    		fprintf(log,"QQ and Manhattan plots saved in \n  %s\n",filenamePlots.c_str());

    		if(PreMeta::recessive)
    		{
    			printf("  %s\n",filenamePlots_rec.c_str());
    			fprintf(log, "  %s\n",filenamePlots_rec.c_str());
    		}

    		if(PreMeta::dominant)
    		{
    			printf("  %s", filenamePlots_dom.c_str());
    			fprintf(log,"  %s", filenamePlots_dom.c_str());
    		}

	  if ( PreMeta::zipOutput ) { // do tabix at last
	  // tabix score
	  	printf( "\nTabixing output .gz files...\n" );
	  	fprintf(log, "\nTabixing output .gz files...\n" );
	  	String cmd = String("tabix -c \"#\" -s 1 -b 2 -e 2 -f ") + filename;
	  	int sys_status = system( cmd.c_str() );
	  	if ( sys_status == 0 ) {
	  		printf( "\nSingle variant stats %s has been tabixed\n",filename.c_str() );
	  		fprintf(log, "\nSingle variant stats %s has been tabixed\n",filename.c_str() );
	  	}
	  	else {
	  		printf("\nUnable to tabix %s\n", filename.c_str());
	  		fprintf(log, "\nUnable to tabix %s\n", filename.c_str());
	  	}
	  // tabix cov		
	  	cmd = String("tabix -c \"#\" -s 1 -b 2 -e 2 -f ") + filename_COV;
	  	sys_status = system( cmd.c_str() );
	  	if ( sys_status == 0 ) {
	  		printf( "\nVar-cov matrix %s has been tabixed\n",filename_COV.c_str() );
	  		fprintf(log, "\nVar-cov matrix %s has been tabixed\n",filename_COV.c_str() );
	  	}
	  	else {
	  		printf("\nUnable to tabix %s\n", filename_COV.c_str());
	  		fprintf(log, "\nUnable to tabix %s\n", filename_COV.c_str());
	  	}  	
	  }

	  printf("\n");
	  fprintf(log,"\n");
	}
	printf("\nLog file listing current options stored in %s\n",logFile.c_str());
	fprintf(log,"\nLog file listing current options stored in %s\n",logFile.c_str());

    //fprintf(log,"\nlog file has been saved in %s\n\n",logFile.c_str());
	time(&now);
	printf("\nAnalysis ends at: %s\n", ctime(&now));
	fprintf(log,"\nAnalysis ends at: %s\n", ctime(&now));
	time_t endTime = time(0);
	int timeSec = (int) (endTime - initialTime);
	double timeHour = timeSec/3600.0;
	printf("Analysis took %d seconds (~ %.1f hours).\n",timeSec, timeHour);
	fprintf(log,"Analysis took %d seconds (~ %.1f hours).\n",timeSec, timeHour);
	fclose(log);
       //fclose(stderr);
}
catch(std::exception& e)
{
  printf("Exiting, exception thrown: %s\n\n",e.what());
  fprintf(log,"Exiting, exception thrown: %s\n\n",e.what());
  PhoneHome::completionStatus("Exception");
  return(-1);
}
PhoneHome::completionStatus("SUCCESS");
return 0;
}
