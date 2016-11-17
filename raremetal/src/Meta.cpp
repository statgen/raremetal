#include "StringArray.h"
#include "Error.h"
#include <stdio.h>
#include "math.h"
#include <iostream>
#include "StringHash.h"
#include "MixChidist.h"
#include "MathSVD.h"
#include "Meta.h"
#include "MetaUtility.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "Calculate_mvt_pvalue.h"
#include "My_mvt.h"
#include "InputFile.h"
#include "SummaryFileReader.h"
#include "QuickIndex.h"
#include <Eigen/SVD>
#include <Eigen/Dense>
#include "QuadProg.h"
#include <iterator>

String Meta::summaryFiles="";
String Meta::covFiles = "";
String Meta::popfile_name = "";
String Meta::prefix = "";
bool Meta::correctGC = false;
bool Meta::Burden = false;
bool Meta::MB = false;
bool Meta::MAB = false;
bool Meta::BBeta = false;
bool Meta::SKAT = false;
//bool Meta::SKATO = false;
bool Meta::VTa = false;
bool Meta::outvcf =false;
bool Meta::fullResult =false;
bool Meta::report = false;
double  Meta::report_pvalue_cutoff = 1e-06;
bool Meta::founderAF = false;
double Meta::HWE = 1e-05;
double Meta::CALLRATE = 0.95;
double Meta::MAF_cutoff = 0.05;
int Meta::marker_col = 2;
int Meta::cov_col = 3;
bool Meta::altMAF = false;
bool Meta::tabix = false;
bool Meta::dosage = false;
bool Meta::RegionStatus = false;
String Meta::Region = "";
String Meta::cond = "";
String Meta::Chr = "";
int Meta::Start = -1;
int Meta::End = -1;
bool Meta::useExactMetaMethod = false; // use Jingjing's exact method
bool Meta::normPop = false; // correct for population stratification

Meta::Meta( FILE * plog )
{
	log = plog;
	Nsamples = -1;
}

Meta::~Meta()
{}

//This function read all study names from a file
//and save the names in the StringArray files
void Meta::Prepare()
{
//	if (useExactMetaMethod) {
//		printf("WARNING: this method is still under development and not recommended to use in this version!!!\n");
//	}
	setMetaStatics();
	openMetaFiles();

	if (normPop)
		FitPgamma();
	
	//if conditional analysis says yes
	if(cond!="")
		prepareConditionalAnalysis();
}

// fit pgamma based on af in each study
// need to read through all score files
void Meta::FitPgamma()
{
	load1kgPopMAF();
	pgamma.Dimension( scorefile.Length(),nPop );
	for(int study=0;study<scorefile.Length();study++) {
		fitpGammaForSingleStudy(study);
	}
}

void Meta::fitpGammaForSingleStudy( int study)
{
	Matrix X; // 1kg mafs
	Vector Y; // sample maf
	int Nmarkers = 2000;
	X.Dimension( Nmarkers, nPop );
	int current_index = 0;

	/** load mafs from sites **/
	String filename = scorefile[study];
	IFILE file;
	file = ifopen(filename,"r");
	if(file == NULL)
		error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.summary.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.summary.file.gz\n",scorefile[study].c_str());
	bool first_line = 1;
	bool adjust;
	while (!ifeof(file)) {	
		String buffer;
		buffer.ReadLine(file);
		if (first_line) {
			tellRvOrRmw( buffer, adjust, marker_col, cov_col );
			first_line = 0;
			continue;
		}
		if(buffer.FindChar('#') != -1)
			continue;
		StringArray tokens;
		tokens.AddTokens(buffer, SEPARATORS);
		if(tokens[0].Find("chr")!=-1)
			tokens[0] = tokens[0].SubStr(3);

		String chr_pos = tokens[0] + ":" + tokens[1];
		String markername = tokens[0] + ":" + tokens[1] + ":" + tokens[2] + ":" + tokens[3];
		bool fail = isDupMarker( tokens[0], chr_pos );
		if (fail)
			continue;
		int c1,c2,c3;
//printf("%s\n",markername.c_str());			
		if(dosage)
		{
			c3 = tokens[4].AsDouble()*tokens[6].AsDouble()*tokens[6].AsDouble();
			c2 = tokens[4].AsDouble()*2.0*tokens[6].AsDouble()*(1.0-tokens[6].AsDouble());
			c1 = tokens[4].AsDouble()*(1.0-tokens[6].AsDouble())*(1.0-tokens[6].AsDouble());
		}
		else
		{
			c1 = tokens[10-adjust].AsDouble();
			c2 = tokens[11-adjust].AsDouble();
			c3 = tokens[12-adjust].AsDouble();
		}	
		if ((tokens[2]=="." && tokens[3]==".") || (tokens[2]==tokens[3] && c1+c2!=0 && c2+c3!=0))
			continue;
		if(tokens[8-adjust].AsDouble()<CALLRATE || tokens[9-adjust].AsDouble()<HWE)
			continue;
		double current_AC;
		int  current_N;
		current_N = c1+c2+c3;
		current_AC = 2*c3+c2;	
//		if (current_AC==0)
//			continue;

		/** save to maf matrix **/
		Y.Push( current_AC/current_N/2 );
//		X[current_index][0] = 1;
		bool status = add1kgMafToX( X, markername, current_index );
		if (!status) { // not found in 1kg marker, use 1/4N
			for(int i=0; i<nPop; i++)
				X[current_index][i] = (double)0.25/current_N;
		}
		current_index++;
		if (current_index >= Nmarkers)
			break;
	}
	ifclose(file);
	if (current_index < Nmarkers)
		X.Dimension( current_index,nPop );
	hashToRemoveDuplicates.Clear();

	// perform regression & save to fk
	setGammaFromRegression( study,X,Y);
}

bool Meta::add1kgMafToX( Matrix & X, String & markername, int current_index )
{
//printf("markername=%s\n",markername.c_str());
	// match allele first
	std::map<String, std::vector<double> >::iterator ptr;
	ptr = af1KG.find(markername);
	if ( ptr == af1KG.end() ) { // check if flip
		StringArray tokens;
		tokens.AddTokens( markername, ":" );
		String flip_markername = tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" + tokens[2];
		ptr = af1KG.find(flip_markername);
		if (ptr==af1KG.end()) { // does not exist...
			return false;
		}
		else { // make a flip vector
			for(int i=0; i<nPop; i++)
				X[current_index][i] = 1 - ptr->second[i];
		}
	}
	else {
		for(int i=0; i<nPop; i++)
			X[current_index][i] = ptr->second[i];
	}
	return true;
}

// do regression, set coefficient
// fit linear model with beta >=0
void Meta::setGammaFromRegression( int study, Matrix & X, Vector & Y )
{
	if (X.rows != Y.Length())
		error("[Meta::setGammaFromRegression] X.row=%d, Y length=%d. Something is wrong...\n",X.rows,Y.Length());
	if (Y.Length()==0)
		error("No available non-zero maf in study #%d. Please check your input file!\n",study+1);

/*print x & y
printf("x:\n");
for(int i=0;i<X.rows;i++) {
	for(int j=0;j<X.cols;j++)
		printf("%g,",X[i][j]);
	printf("\n");
}
printf("\ny:\n");
for(int i=0;i<Y.Length();i++)
	printf("%g,",Y[i]);
printf("\n");
*/


	// (X'X)^-1
	Matrix transX;
	transX.Transpose(X);
	Matrix inv;
	inv.Product(transX,X);
	SVD svd;
	svd.InvertInPlace(inv);
	// X'Y
	Vector XY;
	XY.Dimension(nPop);
	for(int i=0;i<nPop; i++)
		XY[i] = Y.InnerProduct(transX[i]);
//	printf("study=%d",study);
	for(int i=0;i<nPop; i++) {
		pgamma[study][i] = inv[i].InnerProduct(XY);	
//printf(",%g",pgamma[study][i]);
	}	
//printf("\n\n");

/*
	Matrix transX;
	transX.Transpose(X);
	Matrix XfX;
	XfX.Product(transX,X);
	Vector XY;
	XY.Dimension(nPop);
	for(int i=0;i<nPop; i++)
		XY[i] = Y.InnerProduct(transX[i]);

	Matrix Ci;
	Ci.Dimension( X.cols,X.cols,0 );
	for(int i=0;i<X.cols;i++)
		Ci[i][i] = 1;
	Vector Ci0;
	Ci0.Dimension(X.cols,0);

	// note QuadProg requires -y & -ci0 (but ci0 is zero vector)
	Vector betas;
	betas.Dimension(X.cols,0);
	Matrix Ce;
	Ce.Dimension(X.cols,1,1);
	Vector Ce0;
	Ce0.Dimension(1,1);
	solve_quadprog(XfX, XY, Ce, Ce0,Ci,Ci0, betas);
printf("study=%d",study);	
	for(int i=0;i<nPop;i++) {
		pgamma[study][i] = 0 - betas[i];
printf(",%g",pgamma[study][i]);		
	}
printf("\n\n");
*/
}


//this function will read through summary statistics of each study and pool the information.
//At the end, single variant meta-analysis will be completed.
// information will be stored for gene-based analysis
void Meta::PoolSummaryStat(GroupFromAnnotation & group)
{	
	total_N=0;
	flipCount =0;
	skip_count=0;
	
	SampleSize.Dimension( scorefile.Length() );
	for(int s=0; s<scorefile.Length(); s++)
		SampleSize[s] = -1;

	if (useExactMetaMethod) {
		Ydelta.Dimension( scorefile.Length() );
		Ysigma2.Dimension( scorefile.Length() );
		for(int s=0; s<scorefile.Length(); s++) // no need to initialize ydelta
			Ysigma2[s] = -1;
		for(int s=0;s<scorefile.Length();s++) {
			bool y_status = updateYstat(s);
			if (!y_status)
				error("cannot update y stat at study #%d\n", s+1);
		}
		// caculate overall mean & adjust
		int n = 0;
		for(int s=0; s<scorefile.Length(); s++)
			n += SampleSize[s];
		double ymean = 0;
		for(int s=0; s<scorefile.Length(); s++)
			ymean += (double)SampleSize[s] / (double)n * Ydelta[s];
		for(int s=0; s<scorefile.Length(); s++)
			Ydelta[s] -= ymean;
//		residual_adj = 0;
//		for(int s=0; s<scorefile.Length(); s++) {
//			double new_r = (double)(SampleSize[s]-1) * Ysigma2[s] + (double)SampleSize[s]*Ydelta[s]*Ydelta[s];
//			residual_adj += new_r;
//		}
//for(int s=0; s<scorefile.Length(); s++)
//printf("s=%d,Ydelta=%g,Ysigma2=%g\n",s,Ydelta[s],Ysigma2[s]);		
//printf("ymean=%g,ra=%g,n=%d\n",ymean,residual_adj,n);				
//		residual_adj /= (double)(n-1);		
	}

	// pool stats by reading
	for(int study=0;study<scorefile.Length();study++)
	{
		int duplicateSNP=0;
		//Set up data structure to help remove duplicates
		StringIntHash hashToRemoveDuplicates; //hash for skip duplicate markers
		varCount = 0;
		dupcheck_chr = "";
		
		Vector chisq_study_i;
		
		printf("Pooling summary statistics from study %d ...\n",study+1);
		
		//read in summary statistics.
		//maf and summary stat are saved. SNPpos in the file are hashed.
		SummaryFileReader covReader;
		if(cond!="" && cond_status[study]) {
			String covFilename = covfile[study];
			bool cov_status = covReader.ReadTabix(covFilename);
			if (!cov_status)
				error("Cannot open cov file: %s\n", covFilename.c_str());
		}
		
		String filename = scorefile[study];
		IFILE file;
		file = ifopen(filename,"r");
		if(file == NULL)
			error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.summary.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.summary.file.gz\n",scorefile[study].c_str());
		
		String buffer;
		buffer.ReadLine(file);
		bool adjust;
		tellRvOrRmw( buffer, adjust, marker_col, cov_col );
//		if (!status)
//			error("File %s is neither RMW or rvtest score file!\n", filename.c_str());
		ifclose(file);
		FormatAdjust.Push(adjust);
		
		file = ifopen(filename,"r");
		bool pass_header = 0;
		while (!ifeof(file)) {
			buffer.ReadLine(file);
			if (!pass_header) {
				if (SampleSize[study]==-1) {
					if(buffer.Find("##AnalyzedSamples")!=-1) {
						StringArray tokens;
						tokens.AddTokens(buffer,"=");
						SampleSize[study] = tokens[1].AsInteger();
						continue;
					}
				}
				if(buffer.FindChar('#') ==-1 && buffer.Find("CHROM")==-1)
					pass_header = 1;				
			}
			if (!pass_header)
				continue;
			double current_chisq;
			bool status = poolSingleRecord( study, current_chisq, duplicateSNP, adjust, buffer, covReader);
			if (status)
				chisq_study_i.Push(current_chisq);
		}
		ifclose(file);
		
		total_N += SampleSize[study];
		//calculate GC
		chisq_study_i.Sort();
		if (chisq_study_i.size==0) {
			printf("\nWarning: no GC calculated in chisq_study_i!\n");
			GCbyStudy.Push(0);
		}
		else
			GCbyStudy.Push(chisq_study_i[0.5]/0.456);
		if(duplicateSNP>0) {
			printf("\nWarning: %d variants are skipped from analysis due to duplicate records in study %s. Please check log for details.\n\n",duplicateSNP,scorefile[study].c_str());
			fprintf(log,"\nWarning: %d variants are skipped from analysis due to duplicate records in study %s. Please check log for details.\n\n",duplicateSNP,scorefile[study].c_str());
		}
		printf("  done\n");
	}
	if(skip_count>1)
		printf("Warning: %d variants have been excluded from analysis due to un-matched alleles. Please check log file for details.\n",skip_count);
			
	//finalize direction array
	for(int i=0;i<directions.Length();i++)
		while(directions[i].Length()<scorefile.Length())
			directions[i]+='?';
		
	//calculate pooled allele frequencies
	setPooledAF();
		
	printf("\nPerforming Single variant meta analysis ...\n");
	
	//calculate final results here
	IFILE output;
	String filename;	
	IFILE vcfout;
	String vcf_filename;
	printSingleMetaHeader( filename, output );
	if (outvcf)
		printOutVcfHeader( vcf_filename, vcfout );
		
	//for annotation purpose
	target_chr="";
	target_pos = 0, target=0;
	target_pvalue=_NAN_;
	//Sort variants by chr and pos
	for(int i=0;i<SNPmaf_maf.Length();i++)
		printSingleMetaVariant( group, i, output, vcfout );
	
	plotSingleMetaGC( output, 1 );
	if ( cond != "" )
		plotSingleMetaGC( output, 0 );
	printf("\n  done.\n\n");
			
	ifclose(output);
	if(outvcf) {
		ifclose(vcfout);
		printf("\n  VCF file based on superset of variants from pooled studies has been saved \n    %s\n",vcf_filename.c_str());
	}
}
		
void Meta::Run(GroupFromAnnotation & group)
{
	if(outvcf)
		return;
	
	if(!SKAT && !Burden && !MB && !VTa) {
		printf("\nWarning: none of the gene-level tests was requested; only single variant meta-analysis was done.\n");
		return;
	}

	maf = new Vector [group.annoGroups.Length()];
	stats = new Vector [group.annoGroups.Length()];
	
	if(cond!="")
	{
		cond_stats = new Vector [group.annoGroups.Length()];
		cond_cov = new Matrix [group.annoGroups.Length()];
	}
	cov = new Matrix [group.annoGroups.Length()];
	singlePvalue = new Vector [group.annoGroups.Length()];
	singleEffSize = new Vector [group.annoGroups.Length()];
	
	if ( group.annoGroups.Length()<=0 )
		error("No valid groups. Unable to do gene-based test. Check if your group file is empty!\n");
	
	loadSingleStatsInGroup( group );
	
	printf("\nChecking if all groups are analyzed...\n");
	fprintf(log,"\nChecking if all groups are analyzed...\n");
	int geneCounter=0;
	for(int g=0;g<group.annoGroups.Length();g++)
	{
		if(maf[g].Length()==0)
		{
			geneCounter++;
			if(geneCounter==1)
			{
				printf("\t\tPlease check log file for groups that are not analyzed.\n");
				fprintf(log,"Warning: the following groups has no qualifed variants to group and are skipped:\n\t");
			}
			fprintf(log,"%s, ",group.annoGroups[g].c_str());
		}
	}
	printf("\n\tdone!\n");
	fprintf(log,"\n\tdone!\n");
	
	//loop through cov matrices of all studies and update cov
	loadSingleCovInGroup( group );
	
	// update stats in group
	for(int g=0;g<group.annoGroups.Length();g++)
	{
		for(int i=0;i<cov[g].rows;i++)
			for(int j=i+1;j<cov[g].cols;j++)
				cov[g][j][i] = cov[g][i][j];
		if(cond!="")
		{
			for(int i=0;i<cond_cov[g].rows;i++)
				for(int j=i+1;j<cond_cov[g].cols;j++)
					cond_cov[g][j][i] = cond_cov[g][i][j];
		}
	}
	String method = "";
	if(Burden)
	{
		method = "burden";
		BurdenAssoc(method,group,maf,stats,cond_stats,cov,cond_cov,singleEffSize,singlePvalue);
	}
	if(MB)
	{
		method = "MB";
		BurdenAssoc(method,group,maf,stats,cond_stats,cov,cond_cov,singleEffSize,singlePvalue);
	}
	if(MAB) {
		method = "MAB";
		BurdenAssoc(method,group,maf,stats,cond_stats,cov,cond_cov,singleEffSize,singlePvalue);
	}
	if(BBeta) {
		method = "BBeta";
		BurdenAssoc(method,group,maf,stats,cond_stats,cov,cond_cov,singleEffSize,singlePvalue);
	}
	
	if(VTa)
		VTassoc( group );
	
	if(SKAT)
		SKATassoc( group );
	
	pdf.CloseFile();
	
	//housekeeping
	if(maf) delete [] maf;
	if(stats) delete [] stats;
	if(cov) delete [] cov;
	if(singlePvalue) delete [] singlePvalue;
	if(singleEffSize) delete [] singleEffSize;
	printf("\nQQ plots and manhattan polts have been saved in %s.\n",pdf_filename.c_str());
}

//calculate the genotype covariance between G and X
void Meta::CalculateGenotypeCov(SummaryFileReader & covReader,String chr, int pos,int study,Vector & result)
{
	result.Dimension(commonVar_study[study].Length(),0.0);
	for(int cv=0;cv<commonVar_study[study].Length();cv++)
	{
		if (chr != common_chr[commonVar_study[study][cv]])
			continue;
		if (abs(common_pos[commonVar_study[study][cv]] - pos) > 1000000)
			continue;
		if(common_pos[commonVar_study[study][cv]] > pos)
		{
			if(!covReader.ReadRecord(chr,pos))
			continue;
			StringArray tmp;
			tmp.AddTokens(covReader.marker_nearby,",");
			String pos_str;
			pos_str = common_pos[commonVar_study[study][cv]];
			int loc = tmp.Find(pos_str);
			if(loc==-1)
			continue;
			tmp.Clear();
			tmp.AddTokens(covReader.marker_cov,",");
			result[cv] =  tmp[loc].AsDouble()*SampleSize[study];
		}
		else
		{
			int loc = commonVar_markers_in_window[study][cv].Find(pos);
			if(loc==-1)
			continue;
			result[cv] = commonVar_marker_cov[study][cv][loc]*SampleSize[study];
			/*
			if(pos==79760530)
			printf("chr %s pos %d GXcov %g\n",chr.c_str(),pos,commonVar_marker_cov[study][cv][loc]);
		*/    }
	}
}

double Meta::GrabGenotypeCov(SummaryFileReader & covReader,int study,String chr1,String pos1,String chr2,String pos2,String & SNP1, String & SNP2)
{
	double result = 0.0;
	if(chr1!=chr2)
	return result;
	
	int skip = SNPexclude.Integer(study+":"+chr1+":"+pos1);
	if(skip!=-1)
	return result;
	skip = SNPexclude.Integer(study+":"+chr2+":"+pos2);
	if(skip!=-1)
	return result;
	
	StringArray flip,tmp;
	int marker_idx;
	
	if(pos1.AsInteger()<pos2.AsInteger())
	{
		//Check in current hash table first
		marker_idx = covReader.markerPosHash.Integer(chr1 + ":"+pos1);
		//if this record has been hashed
		if(marker_idx!=-1)
		{
			tmp.AddTokens(covReader.markerNearby[marker_idx],",");
			int loc = tmp.Find(pos2);
			if(loc==-1)
				return result;
			tmp.Clear();
			tmp.AddTokens(covReader.markerNearbyCov[marker_idx],",");
			result = tmp[loc].AsDouble()*SampleSize[study];
			return result;
		}
		if(!covReader.ReadRecord(chr1,pos1.AsInteger()))
			return result;
		tmp.Clear();
		tmp.AddTokens(covReader.marker_nearby,",");
		int loc = tmp.Find(pos2);
		if(loc==-1)
			return result;
		tmp.Clear();
		tmp.AddTokens(covReader.marker_cov,",");
		result = tmp[loc].AsDouble()*SampleSize[study];
	}
	else
	{
		//Check in current hash table first
		marker_idx = covReader.markerPosHash.Integer(chr2 + ":"+pos2);
		//if this record has been hashed
		if(marker_idx!=-1)
		{
			StringArray tmp;
			tmp.AddTokens(covReader.markerNearby[marker_idx],",");
			int loc = tmp.Find(pos1);
			if(loc==-1)
			return result;
			tmp.Clear();
			tmp.AddTokens(covReader.markerNearbyCov[marker_idx],",");
			result = tmp[loc].AsDouble()*SampleSize[study];
			return result;
		}
		if(!covReader.ReadRecord(chr2,pos2.AsInteger()))
		return result;
		tmp.Clear();
		tmp.AddTokens(covReader.marker_nearby,",");
		int loc = tmp.Find(pos1);
		if(loc==-1)
		return result;
		tmp.Clear();
		tmp.AddTokens(covReader.marker_cov,",");
		result = tmp[loc].AsDouble()*SampleSize[study];
	}
	
	//if this marker is flipped then markers from the entire row
	//should have cov multiply by -1.0.
	double factor1=1.0,factor2=1.0;
	if(flipSNP.Integer(study+":"+chr1+":"+pos1)!=-1)
	factor1=-1.0;
	if(flipSNP.Integer(study+":"+chr2+":"+pos2)!=-1)
	factor2=-1.0;
	result *= (factor1*factor2);
	
	return result;
}

void Meta::CalculateXXCov(int study,Matrix & result)
{
	for(int i=0;i<result.rows;i++)
	{
		result[i][i] = commonVar_marker_cov[study][i][0]*SampleSize[study];
		for(int j=i+1;j<result.cols;j++)
		{
			if(common_chr[commonVar_study[study][j]] != common_chr[commonVar_study[study][i]])
			continue;
			if(common_pos[commonVar_study[study][i]]<common_pos[commonVar_study[study][j]])
			{
				int loc = commonVar_markers_in_window[study][i].Find(common_pos[commonVar_study[study][j]]);
				if(loc==-1)
				continue;
				result[i][j] = result[j][i] = commonVar_marker_cov[study][i][loc]*SampleSize[study];
			}
			else
			{
				int loc = commonVar_markers_in_window[study][j].Find(common_pos[commonVar_study[study][i]]);
				if(loc==-1)
				continue;
				result[i][j] = result[j][i] = commonVar_marker_cov[study][j][loc]*SampleSize[study];
			}
		}
	}
	SVD svd;
	svd.InvertInPlace(result);
}

/**** for Meta::Prepare() ******/

// set static stat in meta class
void Meta::setMetaStatics()
{
	// region
	if ( Region != "" ) {
		RegionStatus = true;
		printf("Restrict analysis to region %s!\n", Region.c_str());
		StringArray tf;
		tf.AddTokens(Region, ":-");
		Chr = tf[0];
		Start = tf[1].AsInteger();
		End = tf[2].AsInteger();
	}

	// pdf file name
	if(prefix.Last()=='.' || prefix.Last()=='/')
		pdf_filename = prefix + "meta.plots.pdf";
	else
		pdf_filename = prefix + ".meta.plots.pdf";
}


// open pdf
// load list of summary files
// load list of cov files if needed
void Meta::openMetaFiles()
{
	pdf.OpenFile(pdf_filename);

	// summary file
	if(summaryFiles!="") {
		IFILE inFile = ifopen(summaryFiles,"r");
		if(inFile==NULL)
			error("FATAL ERROR! Please check file name for --summaryFiles  %s\n",summaryFiles.c_str());
		String buffer;
		while (!ifeof(inFile)) {
			buffer.ReadLine(inFile);
			if(buffer.FindChar('#')!=-1)
				continue;
			scorefile.AddTokens(buffer, "\n");
		}
		ifclose(inFile);
	}
	else
		error("FATAL ERROR! --studyName can not be empty! Usage: --summaryFiles your.list.of.summary.files\n");
	
	// cov file
	if(covFiles!="") {
		IFILE inFile = ifopen(covFiles,"r");
		if(inFile==NULL)
			error("Cannot open file %s\n",covFiles.c_str());
		String buffer;
		while (!ifeof(inFile)) {
			buffer.ReadLine(inFile);
			if(buffer.FindChar('#')!=-1)
				continue;
			covfile.AddTokens(buffer, "\n");
		}
		ifclose(inFile);
		//check if summary files and cov files have the same length
		if(scorefile.Length()!=covfile.Length())
			error("There are %d summary files and %d covariance files. Please check to make sure the same number of files are included in the list.\n");
	}
	else if(!Burden || !MB || !VTa || !SKAT)
		error("Covariance files are essential to do gene-level tests. Pleaes ues --covFiles your.list.of.cov.files option.\n");

}


// prepare for cond analysis
void Meta::prepareConditionalAnalysis()
{
	// read list of cond variants
	IFILE condFile = ifopen(cond,"r");
	if(condFile==NULL)
		error("Can not open file %s.\n",cond.c_str());

	String buffer;
	while (!ifeof(condFile)) {
		buffer.ReadLine(condFile);
		if(buffer.FindChar('#')!=-1 || buffer.FindChar(':')==-1)
			continue;
		StringArray tmpMarker;
		tmpMarker.AddTokens(buffer, " \t");
		for(int i=0;i<tmpMarker.Length();i++) {
			commonVar.Push(tmpMarker[i]);
			StringArray tmp;
			tmp.AddTokens(tmpMarker[i],":");
			common_chr.Push(tmp[0]);
			common_pos.Push(tmp[1]);
			common_ref.Push(tmp[2]);
			common_alt.Push(tmp[3]);
			conditionVar.SetInteger(tmp[0]+":"+ tmp[1],commonVar.Length()-1);
		}
	}
	ifclose(condFile);
	if(common_pos.Length()==0)
		error("No variant to be conditioned on. Condition file %s might be in the wrong format.\n",cond.c_str());

	// set variable
	int numStudy = scorefile.Length();
	commonVar_study = new IntArray [numStudy];
	commonVar_effSize = new Vector [numStudy];
	commonVar_U = new Vector [numStudy];
	commonVar_V = new Vector [numStudy];
	commonVar_markers_in_window = new IntArray * [numStudy];
	commonVar_marker_cov = new Vector * [numStudy];
	XX_inv = new Matrix [numStudy];
	cond_status = new bool [numStudy];
	commonVar_betaHat = new Vector [numStudy];

	// set marker stats
	bool status = setCondMarkers();
	if (!status) {
		cond = "";
		printf("\nWarning: None of the variants to be conditioned on are found in any of the studies. Conditional analysis options will be ignored.\n\n");
	}
}


// read cond marker info from score & cov file
// if no marker info was load, return false
bool Meta::setCondMarkers()
{
	for(int s=0;s<scorefile.Length();s++) {
		SummaryFileReader statReader,covReader;
		String filename = scorefile[s];
		String covFilename = covfile[s];
			
		if(!statReader.ReadTabix(filename))
			error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.summary.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.summary.file.gz\n",scorefile[s].c_str());
		
		bool adjust;
		String buffer;
		IFILE dup = ifopen(filename,"r");
		buffer.ReadLine(dup);
		tellRvOrRmw( buffer, adjust, marker_col, cov_col );
		ifclose(dup);
		
		if (SampleSize[s]==-1) {
			dup=ifopen(filename,"r");
			while (!ifeof(dup)) {
				buffer.ReadLine(dup);
				//get sample size
				if(buffer.Find("##AnalyzedSamples")!=-1) {
					StringArray tokens;
					tokens.AddTokens(buffer,"=");
					SampleSize[s] = tokens[1].AsInteger();
					break;
				}
			}
			ifclose(dup);
		}
			
		if(!covReader.ReadTabix(covFilename))
			error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.covariance.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.covariance.file.gz.\n",covfile[s].c_str());
			
		for(int i=0;i<commonVar.Length();i++) {
			//if this variant is genotyped in this study
			bool cov_status = covReader.ReadRecord(common_chr[i],common_pos[i]);
			if(!cov_status)
				continue;
			if(!statReader.ReadRecord(common_chr[i],common_pos[i]))
				continue;
			StringArray record;
			record.AddTokens(statReader.buffer,"\t");
				
			if((record[2]==common_ref[i] && record[3]==common_alt[i]) || (record[3]==common_ref[i] && record[2]==common_alt[i])) {
				double v = record[14-adjust].AsDouble();
				if(v>0) {
					commonVar_study[s].Push(i);
					commonVar_V[s].Push(v*v);
					commonVar_effSize[s].Push(record[15-adjust].AsDouble());
					commonVar_U[s].Push(record[13-adjust].AsDouble());
				}
			}
		}
		
		int dim = commonVar_study[s].Length();
		if(dim!=0)
			cond_status[s] =true;
		else {
			cond_status[s] = false;
			continue;
		}
		commonVar_markers_in_window[s] = new IntArray [dim];
		commonVar_marker_cov[s] = new Vector [dim];
			
		bool cov_status = covReader.ReadTabix(covFilename);
		if ( !cov_status )
			error("Cannot open cov file: %s\n", covFilename.c_str());

		for(int i=0;i<dim;i++) {
			int idx = commonVar_study[s][i];
			if(!covReader.ReadRecord(common_chr[idx],common_pos[idx]))
				continue;
			StringArray tmp_markers;
			tmp_markers.AddTokens(covReader.marker_nearby,",");
			for(int j=0;j<tmp_markers.Length();j++)
				commonVar_markers_in_window[s][i].Push(tmp_markers[j].AsInteger());
			tmp_markers.Clear();
			tmp_markers.AddTokens(covReader.marker_cov,",");
			for(int j=0;j<tmp_markers.Length();j++)
				commonVar_marker_cov[s][i].Push(tmp_markers[j].AsDouble());
		}
		XX_inv[s].Dimension(commonVar_study[s].Length(),commonVar_study[s].Length(),0.0);
		CalculateXXCov(s,XX_inv[s]);
		//calculate beta hat for common variants
		commonVar_betaHat[s].Dimension(dim);
		for(int i=0;i<commonVar_betaHat[s].dim;i++)
			commonVar_betaHat[s][i] = XX_inv[s][i].InnerProduct(commonVar_U[s]);
	}

	// check status of cond markers
	bool status = false;
	for(int i=0;i<scorefile.Length();i++) {
		if(cond_status[i]) {
			status=true;
			break;
		}
	}
	return status;	
}

// read RMW output header that records y
// update y and y-mean
bool Meta::updateYstat( int study )
{
	bool found_y = 0;
	int status = 0;
	String fname = scorefile[study];
	IFILE file = ifopen(fname, "r");
	if (file==NULL)
		error("Cannot open file: %s!\n",fname.c_str());
	while( !ifeof(file) ) {
		String buffer;
		buffer.ReadLine(file);
		if (SampleSize[study]==-1) {
			if (buffer.Find("##AnalyzedSamples")!=-1) {
				StringArray tokens;
				tokens.AddTokens(buffer, "=");
				SampleSize[study] = tokens[1].AsInteger();
			}
		}
/*		if (buffer.Find("#AnalyzedTrait") != -1) {
			StringArray tokens;
			tokens.AddTokens(buffer, "\t ");
			Ysigma2[study] = tokens[7].AsDouble();
			status++;
			if (status==2)
				break;
		}*/
		if (buffer.Find("#")==-1)
			break;
		if (found_y) { // name min 25% median 75th max mean(7th) variance(8th)
			StringArray tokens;
			tokens.AddTokens(buffer, "\t ");
			Ydelta[study] = tokens[6].AsDouble();
			Ysigma2[study] = tokens[7].AsDouble();
			status++;
			found_y = 0;
			if (status==1)
				break;
		}
		if (buffer.Find("##TraitSummaries") != -1) {
			found_y = 1;
			continue;
		}
	}
	if (status ==1)
		return 1;
	else
		return 0;
}

void Meta::UpdateDirection(int & direction_idx,int study,char marker_direction,String & chr_pos,bool exclude)
{
	//if this chr_pos is not shown before, but study>0, then assume studies before this study are all monomorphic at this site.
	if(direction_idx==-1)
	{
		if(study>0)
		{
			String tmp = "";
			for(int i=0;i<study;i++)
				tmp+='?';
			tmp += marker_direction;
			directions.Push(tmp);
			directionByChrPos.SetInteger(chr_pos,directions.Length()-1);
		}
		else if(study==0)
		{
			String tmp;
			tmp=marker_direction;
			directions.Push(tmp);
			directionByChrPos.SetInteger(chr_pos,directions.Length()-1);
		}
		if(exclude)
			refalt.Push(".:.");
	}
	else
	{
		int l = directions[direction_idx].Length();
		if(l==study)
			directions[direction_idx] += marker_direction;
		if(l<study)
		{
			for(int i=0;i<study-l;i++)
				directions[direction_idx] += '?';
			directions[direction_idx]+=marker_direction;
		}
		if(l>study)
		{
			printf("WARNING:there is duplicate marker in some of the studies at %s.\n",chr_pos.c_str());
			fprintf(log,"WARNING:there is duplicate marker in some of the studies at %s.\n",chr_pos.c_str());
			directions[direction_idx]+=marker_direction;
		}
	}
}


void Meta::UpdateExcludedMarker( int & study, String & chr_pos,int filter,String markername)
{
	skip_count++;
	String skip_SNP;
	skip_SNP = study;
	skip_SNP += ":";
	skip_SNP +=chr_pos;
	
	SNPexclude.SetInteger(skip_SNP,skip_count);
	if(filter==0)
	{
		fprintf(log,"Warning: variant %s from study %d failed to pass HWE or CALLRATE filter and it is excluded from meta analysis.\n",markername.c_str(),study);
		return;
	}
	if(filter==1)
	{
		fprintf(log,"Warning: variant %s from study %d has at least one allele missing but is polymorphic; and is excluded from meta-analysis.\n",markername.c_str(), study);
		return;
	}
	if(filter==2)
	{
		int idx = directionByChrPos.Integer(chr_pos);
		fprintf(log,"Warning: variant %s from study %d is excluded due to unmatched alleles (expecting to see %s).\n",markername.c_str(),study,refalt[idx].c_str());
		return;
	}
}

void Meta::UpdateStrIntHash(String & chr_pos, int val, StringIntHash & sihash)
{
	int idx = sihash.Find(chr_pos);
	if(idx==-1)
		sihash.SetInteger(chr_pos,val);
	else
	{
		int old_N = sihash.GetCount(idx);
		old_N += val;
		sihash.SetInteger(idx,old_N);
	}	
}

void Meta::UpdateStrDoubleHash(String & chr_pos, double val, StringDoubleHash & sdhash)
{
	int idx = sdhash.Find(chr_pos);
	if(idx==-1)
		sdhash.SetDouble(chr_pos,val);
	else
	{
		double old_N = sdhash.Double(chr_pos);
		old_N += val;
		sdhash.SetDouble(chr_pos,old_N);
	}	
}

void Meta::UpdateACInfo(String & chr_pos,double AC)
{
	int idx = usefulAC.Find(chr_pos);
	if(idx==-1)
		usefulAC.SetDouble(chr_pos,AC);
	else
	{
		int old_N = usefulAC.Double(idx);
		old_N += AC;
		usefulAC.SetDouble(idx,old_N);
	}
}


// simply add up u and v
// for exact method, v passed here is vk/sigmak. V will be further adjusted in print meta record
void Meta::UpdateStats(int study, String & markerName,double stat,double vstat,bool flip)
{
	//update SNPstat
	double flip_factor = 1.0;
	if(flip)
		flip_factor = -1.0;
	vstat = vstat*vstat;
	stat *= flip_factor;
	if (useExactMetaMethod) { // differnt from u & v RMW uses
		stat *= Ysigma2[study];
		vstat *= Ysigma2[study]; // note we're recording Vexact/sigma2 here
	}
	int stat_idx = SNPstat.Find(markerName);
	if(stat_idx<0)
		SNPstat.SetDouble(markerName,stat);
	else
	{
		double prev =  SNPstat.Double(stat_idx);
		prev += stat;
		SNPstat.SetDouble(markerName,prev);
	}
	
	//update SNP_Vstat
	stat_idx = SNP_Vstat.Find(markerName);
	if(stat_idx==-1)
		SNP_Vstat.SetDouble(markerName,vstat);
	else
	{
		double prev =  SNP_Vstat.Double(stat_idx);
		prev += vstat;
		SNP_Vstat.SetDouble(markerName,prev);
	}
}

char Meta::GetDirection(String & chr_pos,double effsize,bool flip)
{
	char direction = '+';
	if(flip)
	{
		if(effsize>0)
		direction = '-';
	}
	else
	{
		if(effsize<0)
		direction = '-';
	}
	return direction;
}

int Meta::MatchTwoAlleles(String refalt_current,int & idx,String & chr_pos )
{
	//if matched an allele pair
	String refalt_hashed = refalt[idx];
	if(refalt_hashed ==".:.")
	{
		refalt[idx]=refalt_current;
		return 0;
	}
	StringArray tmp,tmp_current;
	tmp.AddTokens(refalt_hashed,":");
	tmp_current.AddTokens(refalt_current,":");
	if(tmp[1]==".")
	{
		if(tmp_current[0]!=tmp[0] && tmp_current[1]!=tmp[0])
			return 2;
		else
		{
			refalt[idx] = refalt_current;
			return 0;
		}
	}
	if(refalt_current == refalt_hashed)
		return 0;
	if(tmp[1]+":"+tmp[0] == refalt_current)
		return 1;
	//if previously saved alleles are from monomorphic sites
	if(tmp[1]==tmp[0])
	{
		StringArray tmp2;
		tmp2.AddTokens(refalt_current,":");
		int idx = directionByChrPos.Integer(chr_pos);
		if(tmp2[0]==tmp[0])
		{
			refalt[idx] = refalt_current;
			return 0;
		}
		if(tmp2[1]==tmp[0])
		{
			refalt[idx] = tmp2[1]+":"+tmp2[0];
			return 1;
		}
	}
	return 2;
}

int Meta::MatchOneAllele(String ref_current, int & idx)
{
	//if match ref allele return 0
	// if match alt allele return 1
	// if match nothing, return 2
	if(refalt[idx]==".:.")
	{
		String name;
		name = ref_current +":" + ref_current;
		refalt[idx] = name;
		return 0;
	}
	StringArray tmp;
	tmp.AddTokens(refalt[idx],":");
	if(tmp[0]=="." || tmp[1]==".")
	{
		if(ref_current == tmp[0])
			return 0;
		if(ref_current == tmp[1])
			return 0;
		else
			return 2;
	}
	
	if(tmp[1] == ref_current)
		return 0;
	if(tmp[0] == ref_current)
		return 1;
	return 2;
}


/**************************************************/


/*** pool summary stat related ****/

// pool single record from score file
// if u2/v2 is generated, return true
bool Meta::poolSingleRecord( int study, double & current_chisq, int & duplicateSNP, bool adjust, String & buffer, SummaryFileReader & covReader )
{
	StringArray tokens;
	tokens.AddTokens(buffer, SEPARATORS);
	if (tokens[0].Find("#")!=-1)
		return 0;

	if(tokens[0].Find("chr")!=-1)
		tokens[0] = tokens[0].SubStr(3);
	
	String chr_pos = tokens[0] + ":" + tokens[1];
	int direction_idx = directionByChrPos.Integer(chr_pos);
	
	//CHECK duplicate markers
	bool is_dup = isDupMarker( tokens[0], chr_pos );
	if (is_dup) {
		duplicateSNP++;
		fprintf(log,"Warning: variant %s from study %s is skipped because of duplicate records in the same study.\n",chr_pos.c_str(),scorefile[study].c_str());
		return 0;
	}
	
	//POOLING STEP1: if fail HWE or CALLRATE then skip this record
	//if a variant has a missing allele but not mornomorphic then exclude this variant without updating the total sample size
	//if(((tokens[2]=="." || tokens[3]==".") && (c1+c2!=0 && c2+c3!=0)) || (tokens[2]=="." && tokens[3]==".") || (tokens[2]==tokens[3] && c1+c2!=0 && c2+c3!=0))

	if (tokens[2]=="0")
		tokens[2] = ".";
	if (tokens[3]=="0")
		tokens[3] = ".";

	//check allele counts to see if the site is monomorphic
	int c1,c2,c3;
	if(dosage)
	{
		c3 = tokens[4].AsDouble()*tokens[6].AsDouble()*tokens[6].AsDouble();
		c2 = tokens[4].AsDouble()*2.0*tokens[6].AsDouble()*(1.0-tokens[6].AsDouble());
		c1 = tokens[4].AsDouble()*(1.0-tokens[6].AsDouble())*(1.0-tokens[6].AsDouble());
	}
	else
	{
		c1 = tokens[10-adjust].AsDouble();
		c2 = tokens[11-adjust].AsDouble();
		c3 = tokens[12-adjust].AsDouble();
	}

	bool is_fail = 0;
	int filter_type = 1;
	if ((tokens[2]=="." && tokens[3]==".") || (tokens[2]==tokens[3] && c1+c2!=0 && c2+c3!=0))
		is_fail = 1;
	if(tokens[8-adjust].AsDouble()<CALLRATE || tokens[9-adjust].AsDouble()<HWE) {
	// hwe & call rate filter
		is_fail = 1;
		filter_type = 0;
	}
	
	if (is_fail) {
		char direct = '!';
		UpdateDirection(direction_idx,study,direct,chr_pos,true);
		UpdateStrIntHash( chr_pos, SampleSize[study], usefulSize );
		UpdateExcludedMarker(study,chr_pos,filter_type,tokens[0]+":"+tokens[1]+":"+tokens[2]+":"+tokens[3]);
		return 0;
	}

	double current_AC;
	int  current_N;
	current_N = c1+c2+c3;
	current_AC = 2*c3+c2;
	
	//STEP2: check if this position has been hashed. If yes, match alleles; if not, hash position and ref alt alleles.
	int marker_idx = directionByChrPos.Integer(chr_pos);
	String refalt_current = tokens[2]+':'+ tokens[3];
	String markername = chr_pos + ":" + refalt_current;
	bool flip=false;
	double u = tokens[13-adjust].AsDouble();
	double v = tokens[14-adjust].AsDouble();
	double raw_v = v;
	double raw_u = u;
	if ( useExactMetaMethod ) {
		String snp_no_allele = tokens[0]+":"+tokens[1];
		double raw_af = (double)current_AC / current_N / 2;
		if (normPop) { // adjust for population stratification by regression
				double fk = getAFtilda( markername, raw_af,study );
				addToMapStrVec( variant_fk, study, chr_pos, fk, scorefile.Length()+1);
				addToMapStrVec( variant_nk, study, chr_pos, (double)current_N, scorefile.Length()+1);
		}
		else {
			addToMapStrVec( variant_fk, study, chr_pos, raw_af, scorefile.Length()+1);
			addToMapStrVec( variant_nk, study, chr_pos, (double)current_N, scorefile.Length()+1);
		}
	}
	if(marker_idx == -1) //if this position is never hashed before
	{
		UpdateStrIntHash(chr_pos, current_N, recSize); // add size to --altMAF when site is observed
		if(c2+c3==0||c1+c2==0) 
		{ //if monomorphic, then hash in ref=alt allele for this position and update direction with '?' usefulAC as count of ref allele
			int count=0;
			setRefAltHashKey( refalt_current, tokens, c1, c2 );
			refalt.Push(refalt_current);
			char direct = '?';
			UpdateDirection(direction_idx,study,direct,chr_pos,false);
			//should only hash the count of alternative allele
			UpdateACInfo(chr_pos,count);
			if(SampleSize[study]!=current_N)
				UpdateStrIntHash(chr_pos,SampleSize[study]-current_N, usefulSize);
		}
		else
		{ //if not monomorphic, then hash in ref:alt allele and update direction with alt direction, update usefulSize if count is different from sampleSize, update usefulAC with alt allele count.
			refalt.Push(refalt_current);
			char direct = GetDirection(chr_pos,tokens[15-adjust].AsDouble(),false);
			UpdateDirection( direction_idx,study,direct,chr_pos,false);
			UpdateACInfo(chr_pos,current_AC);
			if(SampleSize[study]!=current_N)
				UpdateStrIntHash(chr_pos,SampleSize[study]-current_N, usefulSize);
			UpdateStats( study, chr_pos,u,v,flip);
		}
	}
	else //if this position has been seen from previous studies
	{
		int match;
		//if the previous records are all monomorphic
		if(c2+c3==0 || c1+c2==0)
		{
			//if monomorphic, then check the following: ref?=alt,if yes, check if ref of this study is the same as the one of the alleles previously hashed for this position (if not the same, then exclude marker;if they are the same update counts and statistics).
			int count=0;
			setPolyMatch( match, chr_pos, refalt_current,tokens, marker_idx );
			// set count from match
			if(tokens[2]=="." || tokens[3]==".") {
				if(match==1)
					count = 0;
				else if(match==0)
					count =0;
			}
			else {
				if(match==0)
					count = 2*c3+c2;
				else if(match==1)
					count = 2*c1+c2;
			}
			
			if(match==2)
			{ //if allels do not match, then exclude this variant and continue
				UpdateExcludedMarker(study,chr_pos,2,tokens[0]+":"+tokens[1]+":"+tokens[2]+":"+tokens[3]);
				UpdateStrIntHash(chr_pos,SampleSize[study],usefulSize);
				char direct = '!';
				UpdateDirection(direction_idx,study,direct,chr_pos,true);
				return 0;
			}
			if(match==1)
				flip=true;
			
			UpdateStrIntHash(chr_pos, current_N, recSize); // add --altMAF
			UpdateACInfo(chr_pos,count);
			if(SampleSize[study]!=current_N)
				UpdateStrIntHash(chr_pos,SampleSize[study]-current_N,usefulSize);
			char direct = '?';
			UpdateDirection( direction_idx,study,direct,chr_pos,false);
		}
		else //if site is not monomorphic,then need to match alleles and upate statistics
		{
			match = MatchTwoAlleles(refalt_current,marker_idx,chr_pos);
			if(match==2)
			{
				UpdateExcludedMarker(study,chr_pos,2,tokens[0]+":"+tokens[1]+":"+tokens[2]+":"+tokens[3]);
				UpdateStrIntHash(chr_pos,SampleSize[study],usefulSize);
				char direct = '!';
				UpdateDirection( direction_idx,study,direct,chr_pos,true);
				return 0;
			}
			UpdateStrIntHash(chr_pos, current_N,recSize);
			if(match==1)
			{
				flip=true;
				if(SampleSize[study]!=current_N)
					UpdateStrIntHash(chr_pos,SampleSize[study]-current_N,usefulSize);
				UpdateACInfo(chr_pos,2*c1+c2);
				char direct = GetDirection(chr_pos,tokens[15-adjust].AsDouble(),true);
				UpdateDirection( direction_idx,study,direct,chr_pos,false);
				UpdateStats( study,chr_pos,u,v,flip);
			}
			if(match==0)
			{
				if(SampleSize[study]!=current_N)
					UpdateStrIntHash(chr_pos,SampleSize[study]-current_N,usefulSize);
				UpdateACInfo(chr_pos,current_AC);
				char direct = GetDirection(chr_pos,tokens[15-adjust].AsDouble(),false);
				UpdateDirection( direction_idx,study,direct,chr_pos,false);
				UpdateStats( study,chr_pos,u,v,flip);
			}
		}
	}
	if (useExactMetaMethod)  { // get it back for calculating gc
		v = raw_v;
		u = raw_u;
	}
	
	if(flip) {
		flipCount++;
		flipSNP.SetInteger(study+":"+chr_pos,flipCount);
	}
	
	//if a variant is monomorphic then update the count of sample size and generate warning
	if(c1+c2==0 || c2+c3==0)
		return 0;
	
	//push the chisq statistics for GC calculation
	if(tokens[14-adjust].AsDouble()>0)
		current_chisq = u*u/(v*v);
	
	//update SNP_cond_stat and SNP_cond_V
	//if(cond!="" && cond_status[study] && tokens[14-adjust].AsDouble()>0.0)
	if(cond!="" && tokens[14-adjust].AsDouble()>0.0) {
		//if this variant is not the one to be conditioned upon
		if(conditionVar.Integer(tokens[0]+":"+tokens[1])==-1)
			updateSNPcond( study, flip, adjust, chr_pos, tokens, covReader );
	}
	return 1;
}

// if dup marker, return true
bool Meta::isDupMarker( String & chr_str, String & chr_pos )
{
	bool is_dup = 0;
	if(chr_str != dupcheck_chr) {
		hashToRemoveDuplicates.Clear();
		dupcheck_chr = chr_str;
	}
	if(hashToRemoveDuplicates.Integer(chr_pos) ==-1) {
		varCount++;
		hashToRemoveDuplicates.SetInteger(chr_pos, varCount);
	}
	else
		is_dup = 1;
	return is_dup;
}

void Meta::setRefAltHashKey( String & refalt_current, StringArray & tokens, int c1, int c2 )
{
	if(tokens[2]=="." || tokens[3]==".") {
		if(tokens[2]==".")
			refalt_current = tokens[3] + ":.";
		else
			refalt_current =tokens[2] + ":.";
	}
	else if(tokens[2]==tokens[3])
		refalt_current  = tokens[2] + ":.";
	else {
		refalt_current = tokens[2] + ':' + tokens[3];
		if(c1+c2==0)
			refalt_current = tokens[3] + ':' + tokens[2];
	}
}

void Meta::setPolyMatch( int & match, String & chr_pos, String & refalt_current, StringArray & tokens, int marker_idx )
{
	if(tokens[2]!=tokens[3]) {
		if(tokens[2]==".")
			match = MatchOneAllele(tokens[3],marker_idx);
		else if (tokens[3]==".")
			match = MatchOneAllele(tokens[2],marker_idx);
		else
			match = MatchTwoAlleles(refalt_current,marker_idx,chr_pos);
	}
	else
		match = MatchOneAllele(tokens[3],marker_idx);
}

//update SNP_cond_stat and SNP_cond_V
void Meta::updateSNPcond( int study, bool flip, int adjust, String & chr_pos, StringArray & tokens, SummaryFileReader & covReader )
{
	Vector GX;
	CalculateGenotypeCov(covReader,tokens[0],tokens[1].AsInteger(),study,GX);
	double cond_u = tokens[13-adjust].AsDouble() - GX.InnerProduct(commonVar_betaHat[study]);
	int idx;
	idx = SNPstat_cond.Find(chr_pos);
	double stat = 0.0;
	if(idx!=-1)
	{
		stat = SNPstat_cond.Double(idx);
		if(!flip)
			stat += cond_u;
		else
			stat += -1.0*cond_u;
		SNPstat_cond.SetDouble(idx,stat);
	}
	else
	{
		if(!flip)
			stat += cond_u;
		else
			stat += -1.0*cond_u;
		SNPstat_cond.SetDouble(chr_pos,stat);
	}
	
	Vector tmp;
	for(int i=0;i<GX.dim;i++)
	tmp.Push(GX.InnerProduct(XX_inv[study][i]));
	double v = tokens[14-adjust].AsDouble();
	double cond_v_part = tmp.InnerProduct(GX);
	double cond_V = v*v-cond_v_part;
	idx = SNP_Vstat_cond.Find(chr_pos);
	if(idx!=-1)
	{
		
		double cond_Vstat = SNP_Vstat_cond.Double(idx);
		cond_Vstat += cond_V;
		SNP_Vstat_cond.SetDouble(idx,cond_Vstat);
	}
	else
		SNP_Vstat_cond.SetDouble(chr_pos,cond_V);	
}


// after reading & pooling stat from score file, set overall af of each marker
// results stored in SNPmaf
void Meta::setPooledAF()
{
	StringArray chr_AC,unique_chr,SNPname_AC;
	IntArray pos_AC;		
	//get the unique chromosomes
	for(int i=0;i<directionByChrPos.Capacity();i++) {
		if(!directionByChrPos.SlotInUse(i))
			continue;
		String SNPname = directionByChrPos[i];
		int idx = directionByChrPos.Integer(SNPname);
		StringArray tmp;
		tmp.AddTokens(SNPname,":");
		chr_AC.Push(tmp[0]);
		pos_AC.Push(tmp[1].AsInteger());
		if(directions[idx].Length()<scorefile.Length()) {
			for(int l=directions[idx].Length();l<scorefile.Length();l++)
				directions[idx] += '?';
		}
		SNPname += ':';
		SNPname += refalt[idx];
		SNPname_AC.Push(SNPname);
		if(unique_chr.Find(tmp[0])==-1)
			unique_chr.Push(tmp[0]);
	}
	QuickIndex chr_AC_idx(chr_AC);
	unique_chr.Sort();
	StringArray chr_cp,character_chr;

	// pool
	for(int i=0;i<unique_chr.Length();i++) {
		if(unique_chr[i].AsInteger()<=22 && unique_chr[i].AsInteger()>=1)
			chr_cp.Push(unique_chr[i]);
		else
			character_chr.Push(unique_chr[i]);
	}
	for(int i=0;i<character_chr.Length();i++)
		chr_cp.Push(character_chr[i]);
	unique_chr = chr_cp; //now unique_chr are sorted as 1,2,...,22,X,Y,M...
	chr_cp.Clear();
	character_chr.Clear();
	for(int i=0;i<unique_chr.Length();i++)
	{
		IntArray pos_i;
		StringArray SNPname_i;
		for(int j=0;j<chr_AC.Length();j++) {
//printf("chr_AC=%d,chr_AC_idx=%d\n",chr_AC.Length(), chr_AC_idx.Length());
			if(chr_AC[chr_AC_idx[j]] == unique_chr[i]) {
				pos_i.Push(pos_AC[chr_AC_idx[j]]);
				SNPname_i.Push(SNPname_AC[chr_AC_idx[j]]);
			}
		}
		QuickIndex pos_i_idx(pos_i);
		for(int j=0;j<pos_i.Length();j++) {
			StringArray tmp;
			tmp.AddTokens(SNPname_i[pos_i_idx[j]],":");
			double AC = usefulAC.Double(tmp[0]+":"+tmp[1]);
				
			// usefulSize has the # of samples to be excluded
			// recSize has # of samples truly added from vcf/ped. Use this when altMAF is toggled
			int N;
			if ( this->altMAF ) {
				N = recSize.Integer( tmp[0]+":"+tmp[1] );
				if ( N==-1 )
					N = 0;
			}
			else { // use default
				N = usefulSize.Integer(tmp[0]+":"+tmp[1]);
				if(N!=-1)
					N = total_N-N;
				else
					N = total_N;
			}
				
			double maf;
			if(founderAF)
				maf = AC/(2.0*N);
			else
				maf = AC/(2.0*N);
				
			int idx = directionByChrPos.Integer(tmp[0]+":"+tmp[1]);
			if(directions[idx].FindChar('+')==-1 && directions[idx].FindChar('-')==-1)
				maf = 0.0;
				
			SNPmaf_maf.Push(maf);
			SNPmaf_name.Push(SNPname_i[pos_i_idx[j]]);
			SNP_effect_N.Push(N);
			SNPmaf.SetInteger(SNPname_i[pos_i_idx[j]],SNPmaf_maf.Length()-1);
		}
	}
}


void Meta::printSingleMetaHeader( String & filename, IFILE & output )
{
	if(prefix =="")
		filename = "meta.singlevar.results";
	else if(prefix.Last()=='.' || prefix.Last()=='/')
		filename = prefix +  "meta.singlevar.results";
	else
		filename = prefix + ".meta.singlevar.results";

	output=ifopen(filename,"w",InputFile::UNCOMPRESSED);	
	ifprintf(output,"##Method=SinglevarScore\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",total_N);
	if(cond=="")
		ifprintf(output,"#CHROM\tPOS\tREF\tALT\tN\tPOOLED_ALT_AF\tDIRECTION_BY_STUDY\tEFFECT_SIZE\tEFFECT_SIZE_SD\tH2\tPVALUE\n");
	else
		ifprintf(output,"#CHROM\tPOS\tREF\tALT\tN\tPOOLED_ALT_AF\tDIRECTION_BY_STUDY\tEFFECT_SIZE\tEFFECT_SIZE_SD\tH2\tPVALUE\tCOND_EFFSIZE\tCOND_EFFSIZE_SD\tCOND_H2\tCOND_PVALUE\n");
}

void Meta::printOutVcfHeader( String & vcf_filename, IFILE & vcfout )
{
	if(prefix=="")
		vcf_filename = "pooled.variants.vcf";
	else if(prefix.Last()=='.' || prefix.Last()=='/')
		vcf_filename = prefix + "pooled.variants.vcf";
	else
		vcf_filename = prefix + ".pooled.variants.vcf";
		
	vcfout = ifopen(vcf_filename,"w",InputFile::UNCOMPRESSED);
	ifprintf(vcfout,"#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}


// print ith record  of single var result
void Meta::printSingleMetaVariant( GroupFromAnnotation & group, int i, IFILE & output, IFILE & vcfout )
{
	if (Nsamples==-1) {
		for (int i=0;i<SNP_effect_N.Length();i++) {
			if (SNP_effect_N[i]>Nsamples)
				Nsamples = SNP_effect_N[i];
		}
		if (Nsamples==-1)
			error("No available SNPs to analysis!");
	}

	String SNPname = SNPmaf_name[i];
	StringArray tmp;
	tmp.AddTokens(SNPname,":");
	String SNPname_noallele = tmp[0]+":"+tmp[1];
	int N = SNP_effect_N[i];
	double U = SNPstat.Double(SNPname_noallele);
	double V = SNP_Vstat.Double(SNPname_noallele);
	double maf = SNPmaf_maf[i];
	if (useExactMetaMethod) {
		std::map< String, std::vector<double> >::iterator pfk = variant_fk.find(SNPname_noallele);
		std::map< String, std::vector<double> >::iterator pnk = variant_nk.find(SNPname_noallele);
		if (pfk==variant_fk.end())
			error("[Meta::printSingleMetaVariant] Cannot find %s!\n",SNPname_noallele.c_str());
		double exact_maf;
		if (normPop) {
		// fill in variant_f
			double f = 0;
			for(int i=0;i<scorefile.Length();i++)
				f += pfk->second[i] / N * pnk->second[i];
			exact_maf = f;
		}
		else
			exact_maf = maf;	
		// here V = sum(vk / sigmak)
		double nkdeltak = 0;
		for(int i=0;i<scorefile.Length();i++)
			nkdeltak += pnk->second[i] * Ydelta[i];
		double nkdeltakfk = 0;
		for(int i=0;i<scorefile.Length();i++)
			nkdeltakfk += pnk->second[i] * pfk->second[i] * Ydelta[i];
		U = U - 2*exact_maf * nkdeltak + 2*nkdeltakfk;
		double v2 = 0;
		for(int i=0;i<scorefile.Length();i++)
			v2 += pnk->second[i] * pfk->second[i]*pfk->second[i];
		double new_r = 0;
		for(int i=0;i<scorefile.Length();i++)
			new_r += (pnk->second[i]-1) * Ysigma2[i] + pnk->second[i] * Ydelta[i]*Ydelta[i];
		new_r /= (N-1);
		V = new_r * (V + 4*v2 - 4 * N * exact_maf * exact_maf);
		// update U and V for gene based test
		pfk->second[scorefile.Length()] = exact_maf;
		pnk->second[scorefile.Length()] = new_r;
		SNPstat.SetDouble( SNPname_noallele,U );
		SNP_Vstat.SetDouble( SNPname_noallele,V );
	}
	int direction_idx = directionByChrPos.Integer(SNPname_noallele);
		
	String direction = directions[direction_idx];
	IntArray pvalue1_idx,pvalue5_idx;


	if(maf==0.0 || maf==1.0) {
		if(cond!="")
			ifprintf(output,"%s\t%s\t%s\t%s\t%d\t%g\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n",tmp[0].c_str(),tmp[1].c_str(),tmp[2].c_str(),tmp[3].c_str(),N,maf,direction.c_str());
		else
			ifprintf(output,"%s\t%s\t%s\t%s\t%d\t%g\t%s\tNA\tNA\tNA\tNA\n",tmp[0].c_str(),tmp[1].c_str(),tmp[2].c_str(),tmp[3].c_str(),N,maf,direction.c_str());
	}
	
	if(maf<=0.0 || maf>=1.0 || V==0.0)
		return;

	double chisq = U*U/V;
//printf("U=%g,V=%g,chisq=%g\n",U,V,chisq);	
	//chisq_before_GC.Push(chisq);
	double pvalue = pchisq(chisq,1,0,0);
	double effSize = U/V;
	double cond_pvalue = _NAN_,cond_effSize = _NAN_;
	bool disect=false;
	double h2 = V*effSize*effSize/N;
	double effSize_se = 1.0/sqrt(V);
	while(pvalue==0.0) {
		disect=true;
		chisq *= 0.999;
		pvalue = pchisq(chisq,1,0,0);
	}
		
	singleVarPvalue.SetDouble(SNPname,pvalue);
	singleVarEff.SetDouble(SNPname,effSize);
		
	pvalueAll.Push(pvalue);
	if(maf<0.01) {
		pvalue1.Push(pvalue);
		pvalue1_idx.Push(pvalueAll.Length()-1);
	}
	if(maf<0.05) {
		pvalue5.Push(pvalue);
		pvalue5_idx.Push(pvalueAll.Length()-1);
	}
	chr_plot.Push(tmp[0]);
	pos_plot.Push(tmp[1].AsInteger());
			
	if(cond!="") { // print conditional analysis results
		bool cond_disect=false;
		double cond_U, cond_V, chisq;
		double cond_h2=_NAN_;
		double cond_effSize_se = _NAN_;
		if(conditionVar.Integer(SNPname_noallele)==-1)
		{
			cond_U = SNPstat_cond.Double(SNPname_noallele);
			cond_V = SNP_Vstat_cond.Double(SNPname_noallele);
			chisq = cond_U*cond_U/cond_V;
			cond_pvalue = pchisq(chisq,1,0,0);
			cond_effSize = cond_U/cond_V;
			while(cond_pvalue==0.0)
			{
				disect=true;
				chisq *= 0.999;
				cond_pvalue = pchisq(chisq,1,0,0);
			}
			cond_h2 = cond_V*cond_effSize*cond_effSize/N;
			cond_effSize_se = 1.0/sqrt(cond_V);
		}
		else
		{
			cond_effSize=0.0;
			cond_pvalue = 1.0;
			cond_h2=0.0;
			cond_effSize_se = 0.0;
		}
		ifprintf(output,"%s\t%s\t%s\t%s\t%d\t%g\t%s\t%g\t%g\t%g\t%s%g\t%g\t%g\t%g\t%s%g\n",tmp[0].c_str(),tmp[1].c_str(),tmp[2].c_str(),tmp[3].c_str(),N,maf,direction.c_str(),effSize,effSize_se,h2,disect?"<":"",pvalue,cond_effSize,cond_effSize_se,cond_h2,cond_disect?"<":"",cond_pvalue);
			
		pvalueAll_cond.Push(cond_pvalue);
		if(maf<0.01)
			pvalue1_cond.Push(cond_pvalue);
		if(maf<0.05)
			pvalue5_cond.Push(cond_pvalue);		
	}
	else
		ifprintf(output,"%s\t%s\t%s\t%s\t%d\t%g\t%s\t%g\t%g\t%g\t%s%g\n",tmp[0].c_str(),tmp[1].c_str(),tmp[2].c_str(),tmp[3].c_str(),N,maf,direction.c_str(),effSize,effSize_se,h2,disect?"<":"",pvalue);

	if(outvcf)
		ifprintf(vcfout,"%s\t%s\t%s\t%s\t%s\t.\t.\tALT_AF=;\n",tmp[0].c_str(),tmp[1].c_str(),tmp[1].c_str(),tmp[2].c_str(),tmp[3].c_str());

	// annotation if needed
	annotateSingleVariantToGene( group, pvalue, cond_pvalue, tmp );
}

//Annotate single variants to gene
void Meta::annotateSingleVariantToGene( GroupFromAnnotation & group, double pvalue, double cond_pvalue, StringArray & tmp )
{
	// initialize
	if(geneLabel.Length()==0) {
		target_chr="";
		target_pos=0;
		target=0;
		target_pvalue = _NAN_;
	}	

	bool skip = 0;
	if(pvalue > 0.05/SNPmaf_maf.Length()) {
		if (cond != "") {
			if ( cond_pvalue > 0.05/SNPmaf_maf.Length() )
				skip = 1;
		}
		else
			skip = 1;
	}
	if (skip) {
		geneLabel.Push("");
		return;
	}

	// only annotate those significant signals
	if(target_chr==tmp[0] && target_pos>tmp[1].AsInteger()-1000000) { //annotated already
		if(pvalue<target_pvalue)
		{ //if this is the higher peak then update the annotation of this locus
			String current_anno = group.AnnotateSingleVar(tmp[0],tmp[1].AsInteger());
			int distance = 1000;
			while(current_anno=="" && distance<=1000000)
			{
				current_anno = group.AnnotateSingleVar(tmp[0],tmp[1].AsInteger()+distance);
				if(current_anno=="")
					current_anno = group.AnnotateSingleVar(tmp[0],tmp[1].AsInteger()-distance);
				distance += 1000;
			}
			if(geneLabel[target]!="" && current_anno!="" && geneLabel[target].Find(current_anno)==-1)
				current_anno = geneLabel[target] + "/" + current_anno;
			if(geneLabel[target]!="" && current_anno=="")
				current_anno = geneLabel[target];				
			geneLabel.Push(current_anno);
			geneLabel[target]="";
			target = geneLabel.Length()-1;
			target_pvalue = pvalue;
			target_pos = tmp[1].AsInteger();
			target_chr = tmp[0];
		}
		else //otherwise, leave the original annotation of this locus
			geneLabel.Push("");
	}
	else
	{
		geneLabel.Push(group.AnnotateSingleVar(tmp[0],tmp[1].AsInteger()));
		target_chr = tmp[0];
		target_pos = tmp[1].AsInteger();
		target_pvalue = pvalue;
		target = geneLabel.Length()-1;
	}
}

// single variant GC plot by MAF
void Meta::plotSingleMetaGC( IFILE & output, bool calc_gc )
{
	String title,demo1,demo2,demo3;
	title = "single variant analysis";
	if (pvalueAll.size <= 0)
		error("No available p values. Something goes wrong?\n");
	double GC1 = GetGenomicControlFromPvalue(pvalueAll);
	demo1="GC=";
	demo1 += GC1;
	double GC2 = 0;
	double GC3 = 0;
	if (pvalue1.size <= 0)
		printf("\nWarning: no MAF<0.01 variants. This group GC = 0!\n");
	else
		GC2 = GetGenomicControlFromPvalue(pvalue1);
	demo2="GC=";
	demo2 += GC2;
	if (pvalue5.size <= 0)
		printf("\nWarning: no MAF<0.05 variants. This group GC = 0!\n");
	else
		GC3 = GetGenomicControlFromPvalue(pvalue5);
	demo3="GC=";
	demo3 += GC3;
	writepdf.Draw(pdf,geneLabel,pvalueAll,pvalue1,pvalue5,chr_plot,pos_plot,title,demo1,demo2,demo3,true);
		
	//Calculate genomic control
	if (calc_gc) {
		ifprintf(output,"#Genomic Control for pooled sample is: %g\n",GC1);
			printf("  Genomic Control for all studies are listed below:\n");
		for(int s=0;s<scorefile.Length();s++) {
			printf("  %g\t",GCbyStudy[s]);
			ifprintf(output,"#Genomic Control for study %d is: %g\n",s,GCbyStudy[s]);
		}
	}
}

/****************************************************************************************/


/*********** for Meta::Run **************/
void Meta::loadSingleStatsInGroup( GroupFromAnnotation & group )
{
	for(int g=0;g<group.annoGroups.Length();g++)
	{
		IntArray del; //del has the SNPs to be deleted from SNPlist
		maf[g].Dimension(0);
		stats[g].Dimension(0);
		int count = group.SNPlist[g].Length();
		singlePvalue[g].Dimension(0);
		singleEffSize[g].Dimension(0);
		
		for(int m=0;m<count;m++)
		{
			String newSNP;
			bool flipStatus = false;
			double af=0.0;
			double singleP = singleVarPvalue.Double(group.SNPlist[g][m]);
			
			if(singleP!=_NAN_)
				af = SNPmaf_maf[SNPmaf.Integer(group.SNPlist[g][m])];
			
			double singleEff = singleVarEff.Double(group.SNPlist[g][m]);
			if(singleP == _NAN_)
			{
				flipStatus=true;
				StringArray SNPflip;
				SNPflip.AddTokens(group.SNPlist[g][m],":");
				newSNP = SNPflip[0]+":"+SNPflip[1]+":"+SNPflip[3]+":"+SNPflip[2];
				singleP = singleVarPvalue.Double(newSNP);
				
				if(singleP!=_NAN_)
					af = SNPmaf_maf[SNPmaf.Integer(newSNP)];
				
				singleEff = singleVarEff.Double(newSNP);
			}
			
			if(singleP==_NAN_)
			{
				del.Push(m);
				fprintf(log,"Warning: variant %s is excluded from group %s for the following reasons:monomorphic or failed QC.\n",group.SNPlist[g][m].c_str(),group.annoGroups[g].c_str());
					continue;
			}
			if(af==0.0 || af==1.0)
			{
				del.Push(m);
				fprintf(log,"Warning: variant %s is excluded from group %s for the following reasons: monomorphic.\n",group.SNPlist[g][m].c_str(),group.annoGroups[g].c_str());
					continue;
			}
			
			if(af>0.5 && 1.0-af>MAF_cutoff)
			{
				del.Push(m);
				fprintf(log,"Warning: variant %s is excluded from group %s for the following reason: maf>cutoff.\n",group.SNPlist[g][m].c_str(),group.annoGroups[g].c_str());
					continue;
			}
			
			if(af<=0.5 && af>MAF_cutoff)
			{
				del.Push(m);
				fprintf(log,"Warning: variant %s is excluded from group %s for the following reason: maf>cutoff.\n",group.SNPlist[g][m].c_str(),group.annoGroups[g].c_str());
					continue;
			}
			
			double tmp;
			StringArray tmp_name;
			tmp_name.AddTokens(group.SNPlist[g][m],":");
			if(flipStatus)
			{
				group.SNPlist[g][m] = newSNP;
				StringArray tmp_name;
				tmp_name.AddTokens(newSNP,":");
				tmp = SNPstat.Double(tmp_name[0]+":"+tmp_name[1]);
			}
			else
				tmp = SNPstat.Double(tmp_name[0]+":"+tmp_name[1]);
			if(af>0.5)
			{
				tmp *= -1.0;
				af = 1.0-af;
				singleEff *= -1.0;
				StringArray tokens;
				tokens.AddTokens(group.SNPlist[g][m],":_/");
				group.SNPlist[g][m] = tokens[0]+":"+tokens[1]+":"+tokens[3]+":"+tokens[2];
				//update flipSNP to get the covaraince right
				for(int s=0;s<scorefile.Length();s++)
				{
					String markername = s+":"+tokens[0]+":"+tokens[1];
					if(flipSNP.Integer(markername)!=-1)
						flipSNP.Delete(markername);
					else {
						flipCount++;
						flipSNP.SetInteger(markername,flipCount);
					}
				}
			}
			stats[g].Push(tmp);
			maf[g].Push(af);
			
			singlePvalue[g].Push(singleP);
			singleEffSize[g].Push(singleEff);
		}
		//delete the SNPs that are not genotyped in any of the studies,
		//or not polymorphic, or SNPs with maf>cutoff.
		for(int i=del.Length()-1;i>=0;i--)
		{
			group.SNPlist[g].Delete(del[i]);
			group.SNPNoAllele[g].Delete(del[i]);
		}
		count = group.SNPlist[g].Length();
		cov[g].Dimension(count,count,0.0);
		if(cond!="")
		{
			cond_stats[g].Dimension(count,0.0);
			cond_cov[g].Dimension(count,count,0.0);
		}
	}
}

//loop through cov matrices of all studies and update cov
void Meta::loadSingleCovInGroup( GroupFromAnnotation & group )
{
	for(int study=0;study<covfile.Length();study++)
	{
		SummaryFileReader covReader;
		StringIntHash markerPosHash;
		StringArray markersInWindow,markersCov;
		String covFilename = covfile[study];
		setFromRvOrRmwAdjust( FormatAdjust[study], marker_col, cov_col );
/*		if(tabix) {
			if(!covReader.ReadTabix(covFilename))
				error("Cannot open file %s.\n",covfile[study].c_str());
		}
		else {*/
		printf("Reading cov matrix from study %d ...\n",study+1);
		String filename = covfile[study];
		IFILE covfile_;
		covfile_  = ifopen(filename,"r");
		if(covfile_ == NULL)
			error("ERROR! Cannot open file: %s! Input cov file has to be bgzipped and tabix indexed using the following command:\n bgzip yourfile.singlevar.cov.txt; tabix -c \"#\" -s 1 -b 2 -e 2 yourprefix.singlevar.cov.txt.gz\n",filename.c_str());
		Tabix covtabix;
		String tabix_name = filename + ".tbi";
		StatGenStatus::Status libstatus = covtabix.readIndex( tabix_name.c_str() );
		if ( RegionStatus ) {
			if ( libstatus != StatGenStatus::SUCCESS )
				error("Cannot open tabix file %s!\n", tabix_name.c_str());
			bool status = SetIfilePosition( covfile_, covtabix, Chr, Start );
			if ( !status )
				error( "Cannot find position %s:%d-%d in cov file %s!\n", Chr.c_str(), Start, End, filename.c_str() );
		}
		String buffer;
		StringArray tokens;
		int m=0;
		while (!ifeof(covfile_)) {
			buffer.ReadLine(covfile_);
			if(buffer.FindChar('#')!=-1 || buffer.Find("CHROM")!=-1)
				continue;
			tokens.Clear();
			tokens.AddTokens(buffer, "\t ");
			if ( RegionStatus ) {
				if ( tokens[1].AsInteger() > End || tokens[0] != Chr ) // out of this region or into another chromosome
					break;
			}
			if(tokens[0].Find("chr")!=-1)
				tokens[0] = tokens[0].SubStr(3);
			String SNP = tokens[0] + ":" + tokens[1];
			//exclude variants that have different REF and ALT
			if(SNPexclude.Integer(study+":"+SNP)!=-1)
				continue;
			m++;
			markerPosHash.SetInteger(SNP,m);
			markersInWindow.Push(tokens[marker_col]);
			markersCov.Push(tokens[cov_col]);
		}
		ifclose(covfile_);
		printf("done\n");

		//   printf("Updating group stats ...\n");
		//update group statistics
		for(int g=0;g<group.annoGroups.Length();g++) {
			//printf("doing group %d\n",g);
			int count = group.SNPlist[g].Length();
			StringArray chr,pos;
			for(int i=0;i<count;i++) {
				StringArray tmp;
				tmp.AddTokens(group.SNPNoAllele[g][i],":");
				chr.Push(tmp[0]);
				pos.Push(tmp[1]);
//				ref_allele.Push(tmp[2]);
//				alt_allele.Push(tmp[2]);
			} //now pos has all the positions of markers in group g.
			
			Matrix cov_i,GX;
			cov_i.Dimension(count,count,0.0);
			if(cond!="") {
				GX.Dimension(count,XX_inv[study].cols,0.0);
				//GX_beta.Dimension(count);
			}
/*			
			if(tabix)
			for(int m=0;m<count;m++) {
				int skip = SNPexclude.Integer(study+":"+group.SNPNoAllele[g][m]);
				if(skip!=-1)
				//printf("SNP %s in stuy %d is skipped %d.\n",group.SNPNoAllele[g][m].c_str(),study,skip);
					continue;
				if(cond!="") {
					StringArray name;
					name.AddTokens(group.SNPlist[g][m],":");
					CalculateGenotypeCov(covReader,chr[m],pos[m].AsInteger(),study,GX[m]);
					//GX_beta[m] = GX[m].InnerProduct(commonVar_betaHat[study]);
					//cond_stats[g][m] = stats[g][m]-GX_beta[m];
					cond_stats[g][m] = SNPstat_cond.Double(name[0]+":"+name[1]);
				}
				for(int s=m;s<count;s++) {
					cov_i[m][s] = GrabGenotypeCov(covReader,study,chr[m],pos[m],chr[s],pos[s],group.SNPlist[g][m],group.SNPlist[g][s]);
				}
			}
			else
*/			
			for(int m=0;m<count;m++) {
				int loc = markerPosHash.Integer(group.SNPNoAllele[g][m]);
				//If this SNP is not genotpyed in this study, then skip
				if(loc==-1)
					continue;
				
				//If this SNP from this study should be excluded due to non consistent ref/alt allele
				//then skip
				int skip = SNPexclude.Integer(study+":"+group.SNPNoAllele[g][m]);
				if(skip!=-1)
				continue;				
				//if this marker is flipped then markers from the entire row
				//should have cov multiply by -1.0.
				int flip = flipSNP.Integer(study+":"+group.SNPNoAllele[g][m]);
				double multiplyFactor=1.0;
				if(flip!=-1)
					multiplyFactor=-1.0;
				//read through markersInWindow and find the selected markers
				StringArray markers,markerscov;
				markers.AddTokens(markersInWindow[loc-1],",");
				markerscov.AddTokens(markersCov[loc-1],",");
				
				//fill in the mth row of cov
				for(int s=m;s<count;s++) {
					int p = markers.Find(pos[s]);
					if(p==-1)
						continue;
					else {
						String markerName = study+":"+chr[s]+":"+pos[s];
						//String markername = markerName + ":"+ref_allele[s] + ":"+alt_allele[s];
						//if the marker in window is supposed to be excluded due to non-consistent ref/alt allele
						//then skip
						int skip = SNPexclude.Integer(markerName);
						if(skip!=-1)
							continue;
						int flip = flipSNP.Integer(markerName);
						double factor=1.0;
						if(flip!=-1)
							factor=-1.0;
						cov_i[m][s]= multiplyFactor*factor*markerscov[p].AsDouble()*SampleSize[study];
						if (useExactMetaMethod) {
							String markername = chr[m] + ":" + pos[m];
							String markername2 = chr[s] + ":" + pos[s];
							std::map< String, std::vector<double> >::iterator pnk = variant_nk.find(markername);
							/*if (p == variant_nk.end()) {
								markername = markerName + ":" + alt_allele[s] + ":"+ref_allele[s];
								p = variant_nk.find(markername);
							}*/
							if( pnk==variant_nk.end())
								error("error2137\n");
							std::map< String, std::vector<double> >::iterator pfk = variant_fk.find(markername);
							std::map< String, std::vector<double> >::iterator pnk2;
							std::map< String, std::vector<double> >::iterator pfk2;
							if (m==s) {
								pnk2 = pnk;
								pfk2 = pfk;
							}
							else {
								pnk2 = variant_nk.find(markername2);
								pfk2 = variant_fk.find(markername2);
							}
							int nk = pnk->second[study];
							double nkfk2 = nk * pfk->second[study] * pfk2->second[study];
							double new_r = pnk->second[scorefile.Length()];
							double ff = pfk->second[scorefile.Length()]*pfk2->second[scorefile.Length()];
							// Vexact=Vrmw*sigma4
							cov_i[m][s] = new_r*(cov_i[m][s]*Ysigma2[study] - 4*nk*ff + 4*nkfk2);
						}
					}
				}
				//fill in GX
				if(cond!="")
				{
					StringArray name;
					name.AddTokens(group.SNPlist[g][m],":");
					cond_stats[g][m] = SNPstat_cond.Double(name[0]+":"+name[1]);
					String pos_str;
					for(int s=0;s<commonVar_study[study].Length();s++)
					{
						if(pos[m].AsInteger()<common_pos[commonVar_study[study][s]]) { // direct fill in
							pos_str = common_pos[commonVar_study[study][s]];
							int p = markers.Find(pos_str);
							if(p==-1)
								continue;
							GX[m][s]= markerscov[p].AsDouble()*SampleSize[study];
						}
						else { // need to put the smaller one first
							int loc = markerPosHash.Integer(common_chr[commonVar_study[study][s]]+":"+common_pos[commonVar_study[study][s]]);
							//If this SNP is not genotpyed in this study, then skip
							if(loc==-1)
								continue;
							//read through markersInWindow and find the selected markers
							StringArray new_markers, new_markerscov;
							new_markers.AddTokens(markersInWindow[loc-1],",");
							new_markerscov.AddTokens(markersCov[loc-1],",");
							pos_str = pos[m];
							int p = new_markers.Find(pos_str);
							if (p==-1)
								continue;
							GX[m][s]= new_markerscov[p].AsDouble()*SampleSize[study];
						}
					}
				}
			}
			cov[g].Add(cov_i);
			if(cond!="")
			{
				cond_cov[g].Add(cov_i);
				Matrix GX_trans,tmp,extra_cov_i;
				GX_trans.Transpose(GX);
				tmp.Product(GX,XX_inv[study]);
				extra_cov_i.Product(tmp,GX_trans);
				extra_cov_i.Multiply(-1.0);
				cond_cov[g].Add(extra_cov_i);
			}
		}
	}	
}

void Meta::BurdenAssoc(String method,GroupFromAnnotation & group,Vector *& maf,Vector *& stats,Vector *& cond_stats,Matrix *& cov,Matrix *& cond_cov,Vector *& singleEffSize,Vector *& singlePvalue)
{
	printf("Performing %s tests ...\n",method.c_str());
	//calculate final results here
	
	Vector pvalue_burden,pvalue_burden_cond;
	
	IFILE output;
	String filename;
	openMetaResultFile( prefix, filename, output, method );
	
	String method_out = method;
	method_out +="_";
	method_out += MAF_cutoff;
	
	IFILE reportOutput;
	String reportFile;
	if(report)
	{
		if(prefix =="")
			reportFile = "meta.tophits."+method+".tbl";
		else if(prefix.Last()=='.' || prefix.Last()=='/')
			reportFile = prefix +  "meta.tophits."+method+".tbl";
		else
			reportFile = prefix + ".meta.tophits."+method+".tbl";
		reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
		if(cond!="")
			ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\tCOND_EFFSIZES\tCOND_PVALUES\n");
		else
			ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\n");
	}
	ifprintf(output,"##Method=Burden\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",total_N);
	if(cond!="")
		if(fullResult)
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
		else
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
	else
		if(fullResult)
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
		else
			ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
	double numerator=_NAN_,denominator=_NAN_,chisq=0.0,pvalue=_NAN_,cond_num=_NAN_,cond_denom=_NAN_;
	StringArray chr_plot,geneLabels;
	Vector pos_plot;
	for(int g=0;g<group.annoGroups.Length();g++)
	{
		if(maf[g].Length()==0)
		continue;
		
		double average_af = maf[g].Average();
		double min_af = maf[g].Min();
		double max_af = maf[g].Max();
		
		String var;
		for(int i=0;i<maf[g].Length()-1;i++)
			var += group.SNPlist[g][i] + ";";
		var += group.SNPlist[g][maf[g].Length()-1];
		
		Vector weight;
		weight.Dimension(maf[g].Length());
		SetWeight(method,weight,maf[g]);
		// for burden test, need to 1/w
		for(int w=0;w<weight.Length();w++)
			weight[w] = 1/weight[w];
		
		numerator  = weight.InnerProduct(stats[g]);
		Vector tmp;
		tmp.Dimension(group.SNPlist[g].Length());
		
		for(int i=0;i<tmp.Length();i++)
			tmp[i] = weight.InnerProduct(cov[g][i]);
		denominator = tmp.InnerProduct(weight);
			
		if(cond!="") {
			cond_num = weight.InnerProduct(cond_stats[g]);
			for(int i=0;i<tmp.Length();i++)
				tmp[i] = weight.InnerProduct(cond_cov[g][i]);
			cond_denom = tmp.InnerProduct(weight);
		}
		
		if(denominator==0.0)
			continue;
		
		chisq = numerator*numerator/denominator;
		pvalue = pchisq(chisq,1,0,0);
		double effSize = numerator/denominator;
		double cond_chisq=_NAN_,cond_effSize=_NAN_,cond_pvalue=_NAN_;

		bool disect=false;
		while(pvalue==0.0)
		{
			disect=true;
			chisq *= 0.999;
			pvalue = pchisq(chisq,1,0,0);
		}
		bool cond_disect =false;
		if(cond!="")
		{
			if(cond_denom==0)
			{
				cond_effSize =0.0;
				cond_pvalue =1.0;
			}
			else
			{
				cond_chisq = cond_num*cond_num/cond_denom;
				cond_effSize = cond_num/cond_denom;
				cond_pvalue = pchisq(cond_chisq,1,0,0);
				while(cond_pvalue==0.0)
				{
					cond_disect=true;
					cond_chisq *= 0.999;
					cond_pvalue = pchisq(cond_chisq,1,0,0);
				}
			}
			pvalue_burden_cond.Push(cond_pvalue);
		}
		
		if(fullResult)
		{
			ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());
			
			for(int i=0;i<maf[g].Length()-1;i++)
				ifprintf(output,"%g,",maf[g][i]);
			ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);
			
			for(int i=0;i<singleEffSize[g].Length()-1;i++)
				ifprintf(output,"%g,",singleEffSize[g][i]);
			ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);
			
			for(int i=0;i<singlePvalue[g].Length()-1;i++)
				ifprintf(output,"%g,",singlePvalue[g][i]);
			ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
			
			if(cond!="")
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n",average_af,min_af,max_af,effSize,disect?"<":"",pvalue,cond_effSize,cond_disect?"<":"",cond_pvalue);
			else
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\n",average_af,min_af,max_af,effSize,disect?"<":"",pvalue);
		}
		else
		{
			if(cond!="")
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,effSize,disect?"<":"",pvalue,cond_effSize,cond_disect?"<":"",cond_pvalue);
			else
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,effSize,disect?"<":"",pvalue);
		}
		if(pvalue <report_pvalue_cutoff && report)
		{
			StringArray variants;
			variants.AddTokens(var,";");
			if(cond!="")
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect?"<":"",pvalue,cond_disect?"<":"",cond_pvalue,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
			else
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect?"<":"",pvalue,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
		}
		pvalue_burden.Push(pvalue);
		geneLabels.Push(group.annoGroups[g]);
		StringArray tmp_SNPname;
		tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
		chr_plot.Push(tmp_SNPname[0]);
		pos_plot.Push(tmp_SNPname[1].AsDouble());
	}
	String name = method;
	name += " (maf<";
	name += MAF_cutoff;
	name +=  ")";
	String extraname = "";
	String demo;
	if(pvalue_burden.Length()>0)
	{
		//Calculate genomic control
		double GC = GetGenomicControlFromPvalue(pvalue_burden);
		demo="GC = ";
		demo += GC;
		writepdf.Draw(pdf,geneLabels,pvalue_burden,chr_plot,pos_plot,name,extraname,demo,true);
		if(cond!="")
		{
			name += "conditional analysis";
			double GC = GetGenomicControlFromPvalue(pvalue_burden_cond);
			demo="GC = ";
			demo += GC;
			writepdf.Draw(pdf,geneLabels,pvalue_burden_cond,chr_plot,pos_plot,name,extraname,demo,true);
		}
	}
	ifclose(output);
	if(report)
		ifclose(reportOutput);
	printf("  done.\n\n");
}

void Meta::VTassoc( GroupFromAnnotation & group )
{
	printf("Performing Variable Threshold tests ...\n");
	//calculate final results here
	Vector pvalue_VT,pos_plot,cond_pvalue_VT;
	StringArray chr_plot,geneLabels;
	
	IFILE output;
	String filename;
	String method = "VT_";
	openMetaResultFile( prefix, filename, output, method );
	
	method += MAF_cutoff;
	IFILE reportOutput;
	if(report)
	{
		String reportFile;
		if(prefix =="")
		reportFile = "meta.tophits.VT.tbl";
		else if(prefix.Last()=='.' || prefix.Last()=='/')
		reportFile = prefix +  "meta.tophits.VT.tbl";
		else
			reportFile = prefix + ".meta.tophits.VT.tbl";
		reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
		ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\n");
	}
	
	ifprintf(output,"##Method=VT\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",total_N);
	if(fullResult)
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\t");
	else
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\t");
	
	if(cond!="")
		ifprintf(output,"EFFECT_SIZE\tMAF_CUTOFF\tCOND_PVALUE\n");
	else
		ifprintf(output,"\n");
	
	for(int g=0;g<group.annoGroups.Length();g++)
	{
		if(g>1 && g%1000==1)
		printf("Finished analyzing %d genes.\n",g-1);
		
		if(maf[g].Length()==0)
		continue;
		
		if(maf[g].Length()==1)
		{
			if(fullResult)
				ifprintf(output,"%s\t1\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t",group.annoGroups[g].c_str(),group.SNPlist[g][0].c_str(),maf[g][0],singleEffSize[g][0],singlePvalue[g][0],maf[g][0],maf[g][0],maf[g][0],singleEffSize[g][0],maf[g][0],singlePvalue[g][0]);
			else
				ifprintf(output,"%s\t1\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t",group.annoGroups[g].c_str(),group.SNPlist[g][0].c_str(),maf[g][0],maf[g][0],maf[g][0],singleEffSize[g][0],maf[g][0],singlePvalue[g][0]);
			pvalue_VT.Push(singlePvalue[g][0]);
			
			if(cond!="")
			{
				String SNPname_noallele;
				StringArray tmp;
				tmp.AddTokens(group.SNPlist[g][0],":");
				SNPname_noallele=tmp[0]+":"+tmp[1];
				double cond_pvalue_=_NAN_,cond_U=_NAN_,cond_V=_NAN_,chisq=_NAN_,cond_effSize_=_NAN_;
				bool disect=false;
				if(conditionVar.Integer(SNPname_noallele)==-1)
				{
					cond_U = SNPstat_cond.Double(SNPname_noallele);
					cond_V = SNP_Vstat_cond.Double(SNPname_noallele);
					chisq = cond_U*cond_U/cond_V;
					cond_pvalue_ = pchisq(chisq,1,0,0);
					cond_effSize_ = cond_U/cond_V;
					while(cond_pvalue_==0.0)
					{
						disect=true;
						chisq *= 0.999;
						cond_pvalue_ = pchisq(chisq,1,0,0);
					}
				}
				else {
					cond_effSize_=0.0;
					cond_pvalue_ = 1.0;
				}
				cond_pvalue_VT.Push(cond_pvalue_);
				ifprintf(output,"%g\t%g\t%s%g",cond_effSize_,maf[g][0],disect?"<":"",cond_pvalue_);
			}
			ifprintf(output,"\n");
			StringArray tmp_SNPname;
			tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
			chr_plot.Push(tmp_SNPname[0]);
			pos_plot.Push(tmp_SNPname[1].AsDouble());
			geneLabels.Push(group.annoGroups[g]);
			continue;
		}
		
		//STEP1: sort maf[g] to find # of cutoffs
		Vector cp_maf;
		cp_maf.Copy(maf[g]);
		cp_maf.Sort();
		Vector maf_cutoff;
		maf_cutoff.Push(cp_maf[0]);
		
		for(int i=1;i<cp_maf.Length();i++)
		{
			if(cp_maf[i]>maf_cutoff[maf_cutoff.Length()-1])
				maf_cutoff.Push(cp_maf[i]);
		} //now unique maf cutoffs are saved in maf_cutoff.
		double pvalue = _NAN_,cond_pvalue = _NAN_;
		bool condition_status = false;
		pvalue = VTassocSingle(group,maf_cutoff,reportOutput,output,g,condition_status,method);
		if(cond!="") {
			condition_status=true;
			cond_pvalue = VTassocSingle(group,maf_cutoff,reportOutput,output,g,condition_status,method);
		}
		pvalue_VT.Push(pvalue);
		if(cond!="")
			cond_pvalue_VT.Push(cond_pvalue);
		
		StringArray tmp_SNPname;
		tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
		chr_plot.Push(tmp_SNPname[0]);
		pos_plot.Push(tmp_SNPname[1].AsDouble());
		geneLabels.Push(group.annoGroups[g]);
	}
	
	String name = "VT (maf<";
	name +=  MAF_cutoff;
	name +=  ")";
	String extraname = "";
	String demo="";
	double GC = GetGenomicControlFromPvalue(pvalue_VT);
	demo="GC=";
	demo += GC;
	writepdf.Draw(pdf,geneLabels,pvalue_VT,chr_plot,pos_plot,name,extraname,demo,true);
	if(cond!="")
	{
		name += " Conditional Analysis";
		double GC = GetGenomicControlFromPvalue(cond_pvalue_VT);
		demo="GC=";
		demo += GC;
		writepdf.Draw(pdf,geneLabels,cond_pvalue_VT,chr_plot,pos_plot,name,extraname,demo,true);
	}
	
	ifclose(output);
	if(report)
		ifclose(reportOutput);
	printf("Done.\n\n");
}

double Meta::VTassocSingle(GroupFromAnnotation & group, Vector & maf_cutoff, IFILE reportOutput, IFILE output, int & g,bool condition,String & method)
{
	double pvalue=_NAN_,chosen_cutoff=_NAN_,chosen_effSize=_NAN_;
	double numerator=0.0,denominator=0.0,t_max=_NAN_;
	Vector weight,tmp,chosen_weight,score;
	Matrix cov_weight;
	weight.Dimension(maf[g].Length());
	tmp.Dimension(group.SNPlist[g].Length());
	cov_weight.Dimension(maf_cutoff.Length(),maf[g].Length());
	score.Dimension(maf_cutoff.Length());
	for(int i=0;i<maf_cutoff.Length();i++)
	{
		for(int w=0;w<weight.Length();w++)
		{
			if(maf[g][w]<=maf_cutoff[i])
				weight[w]=1.0;
			else
				weight[w]=0.0;
			cov_weight[i][w]=weight[w];
		}
		if(condition)
			numerator = weight.InnerProduct(cond_stats[g]);
		else
			numerator = weight.InnerProduct(stats[g]);
		for(int d=0;d<tmp.Length();d++)
		{
			if(condition)
				tmp[d] = weight.InnerProduct(cond_cov[g][d]);
			else
				tmp[d] = weight.InnerProduct(cov[g][d]);
		}
		denominator = tmp.InnerProduct(weight);
		
		if(denominator != 0.0)
		{
			double t_stat = fabs(numerator/sqrt(denominator));
			score[i]=t_stat;
			if(t_max==_NAN_)
			{
				t_max = t_stat;
				chosen_cutoff = maf_cutoff[i];
				chosen_weight.Copy(weight);
				chosen_effSize = numerator/denominator;
			}
			else
			{
				if(t_stat>t_max)
				{
					t_max = t_stat;
					chosen_cutoff = maf_cutoff[i];
					chosen_weight.Copy(weight);
					chosen_effSize = numerator/denominator;
				}
			}
		}
		else
			score[i]=0.0;
	}
	if(score.Max()==0.0)
	{
		printf("Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
		fprintf(log,"Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
		return pvalue;
	}
	Vector tmp_maf,tmp_eff,tmp_pvalue;
	for(int i=0;i<maf[g].Length();i++)
	{
		if(chosen_weight[i]==1.0)
			tmp_maf.Push(maf[g][i]);
	}
	for(int i=0;i<maf[g].Length();i++)
	{
		if(chosen_weight[i]==1.0)
		{
			tmp_eff.Push(singleEffSize[g][i]);
			tmp_pvalue.Push(singlePvalue[g][i]);
		}
	}
	
	double average_af = tmp_maf.Average();
	double min_af = tmp_maf.Min();
	double max_af = tmp_maf.Max();
	String var;
	for(int i=0;i<maf[g].Length()-1;i++)
	{
		if(chosen_weight[i]==1.0)
			var += group.SNPlist[g][i] + ";";
	}
	if(chosen_weight[maf[g].Length()-1]==1.0)
		var += group.SNPlist[g][maf[g].Length()-1];
	
	//STEP3: calculate covariance matrix for (U_1 ... U_#cutoff)
	Matrix cov_U,cov_U_tmp;
	if(condition)
		cov_U_tmp.Product(cov_weight,cond_cov[g]);
	else
		cov_U_tmp.Product(cov_weight,cov[g]);
	Matrix cov_weight_trans;
	cov_weight_trans.Transpose(cov_weight);
	cov_U.Product(cov_U_tmp,cov_weight_trans); //now, cov(U) is saved in cov_U
	//Calculate covariance matrix for (T_1 ... T_#cutoff)
	Matrix cov_T;
	cov_T.Dimension(cov_U.rows,cov_U.cols);
	cov2cor(cov_U,cov_T);
	//STEP4: calculate VT pvalue and report.
	int cutoff = maf_cutoff.Length();
	double * lower = new double [cutoff];
	double * upper = new double [cutoff];
	double * mean = new double [cutoff];
	
	for(int i=0;i<cutoff;i++)
	{
		mean[i] = 0.0;
		lower[i] = -t_max;
		upper[i] = t_max;
	}
	
	//Use pmvnorm to calculate the asymptotic p-value
	Vector result;
	pmvnorm(lower,upper,mean,cov_T,false,result);
	
	if(result[0]==-1.0)
	{
		if(!condition)
		{
			if(cond!="")
			{
				if(fullResult)
				{
					ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str());
					
					for(int i=0;i<tmp_maf.Length()-1;i++)
						ifprintf(output,"%g,",tmp_maf[i]);
					ifprintf(output,"%g\t",tmp_maf[tmp_maf.Length()-1]);
					
					for(int i=0;i<tmp_eff.Length()-1;i++)
						ifprintf(output,"%g,",tmp_eff[i]);
					ifprintf(output,"%g\t",tmp_eff[tmp_eff.Length()-1]);
					
					for(int i=0;i<tmp_pvalue.Length()-1;i++)
						ifprintf(output,"%g,",tmp_pvalue[i]);
					ifprintf(output,"%g\t",tmp_pvalue[tmp_pvalue.Length()-1]);
					
					ifprintf(output,"%g\t%g\t%g\t%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\t",average_af,min_af,max_af,chosen_effSize,chosen_cutoff);
				}
				else
					ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str(),average_af,min_af,max_af,chosen_effSize,chosen_cutoff);
			}
			else
			{
				ifprintf(output,"\n");
			}
		}
		else
		{
			ifprintf(output,"%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\n",chosen_effSize,chosen_cutoff);
		}
	}
	else
	{
		if(1.0-result[0]==0.0)
		{
			//           printf("gene %s has result %g\n",group.annoGroups[g].c_str(),1.0-result[0]);
			printf("Using Shuang's algorithm to calculate MVN pvalue for gene %s ... ",group.annoGroups[g].c_str());
				if(maf_cutoff.Length()>20)
			{
				while(maf_cutoff.Length()>20)
					maf_cutoff.Delete(0);
				
				double numerator,denominator,t_max=_NAN_;
				Vector weight,tmp,chosen_weight,score;
				Matrix cov_weight;
				weight.Dimension(maf[g].Length());
				tmp.Dimension(group.SNPlist[g].Length());
				cov_weight.Dimension(maf_cutoff.Length(),maf[g].Length());
				for(int i=0;i<maf_cutoff.Length();i++)
				{
					for(int w=0;w<weight.Length();w++)
					{
						if(maf[g][w]<=maf_cutoff[i])
							weight[w]=1.0;
						else
							weight[w]=0.0;
						cov_weight[i][w]=weight[w];
					}
					if(condition)
						numerator = weight.InnerProduct(cond_stats[g]);
					else
						numerator = weight.InnerProduct(stats[g]);
					
					for(int d=0;d<tmp.Length();d++)
					{
						if(condition)
							tmp[d] = weight.InnerProduct(cond_cov[g][d]);
						else
							tmp[d] = weight.InnerProduct(cov[g][d]);
					}
					denominator = tmp.InnerProduct(weight);
					if(denominator != 0)
					{
						double t_stat = fabs(numerator/sqrt(denominator));
						score.Push(t_stat);
						if(t_max==_NAN_)
						{
							t_max = t_stat;
							chosen_cutoff = maf_cutoff[i];
							chosen_weight.Copy(weight);
							chosen_effSize = numerator/denominator;
						}
						else
						{
							if(t_stat>t_max)
							{
								t_max = t_stat;
								chosen_cutoff = maf_cutoff[i];
								chosen_weight.Copy(weight);
								chosen_effSize = numerator/denominator;
							}
						}
					}
					else
						score.Push(0.0);
				}
				if(score.Max()==0.0)
				{
					printf("Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
					fprintf(log,"Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
					return pvalue;
					printf("completed!\n");
				}
				Vector tmp_maf,tmp_eff,tmp_pvalue;
				for(int i=0;i<maf[g].Length();i++)
				{
					if(chosen_weight[i]==1.0)
						tmp_maf.Push(maf[g][i]);
				}
				
				for(int i=0;i<maf[g].Length();i++)
				{
					if(chosen_weight[i]==1.0)
					{
						tmp_eff.Push(singleEffSize[g][i]);
						tmp_pvalue.Push(singlePvalue[g][i]);
					}
				}
				average_af = tmp_maf.Average();
				min_af = tmp_maf.Min();
				max_af = tmp_maf.Max();
				
				String var;
				for(int i=0;i<maf[g].Length()-1;i++)
				{
					if(chosen_weight[i]==1.0)
						var += group.SNPlist[g][i] + ";";
				}
				if(chosen_weight[maf[g].Length()-1]==1.0)
					var += group.SNPlist[g][maf[g].Length()-1];
				//STEP3: calculate covariance matrix for (U_1 ... U_#cutoff)
				Matrix cov_U,cov_U_tmp;
				if(condition)
					cov_U_tmp.Product(cov_weight,cond_cov[g]);
				else
					cov_U_tmp.Product(cov_weight,cov[g]);
				Matrix cov_weight_trans;
				cov_weight_trans.Transpose(cov_weight);
				cov_U.Product(cov_U_tmp,cov_weight_trans); //now, cov(U) is saved in cov_U
				//Calculate covariance matrix for (T_1 ... T_#cutoff)
				Matrix cov_T;
				cov_T.Dimension(cov_U.rows,cov_U.cols);
				cov2cor(cov_U,cov_T);
				
				pvalue = CalculateMVTPvalue(score,cov_T,t_max);
				printf("completed!\n");
			}
			else
			{
				pvalue = CalculateMVTPvalue(score,cov_T,t_max);
				printf("completed!\n");
			}
		}
		else
		{
			pvalue = 1.0 - result[0];
		}
		
		if((condition && cond!="") || cond=="")
		{
			if(pvalue <report_pvalue_cutoff && report)
			{
				StringArray variants;
				variants.AddTokens(var,";");
				for(int v=0;v<tmp_maf.Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),pvalue,MAF_cutoff,chosen_cutoff,variants[v].c_str(),tmp_maf[v],tmp_eff[v],tmp_pvalue[v]);
			}
		}
		
		if(cond=="" || (!condition && cond!=""))
		{
			if(fullResult)
			{
				ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str());
				
				for(int i=0;i<tmp_maf.Length()-1;i++)
					ifprintf(output,"%g,",tmp_maf[i]);
				ifprintf(output,"%g\t",tmp_maf[tmp_maf.Length()-1]);
				
				for(int i=0;i<tmp_eff.Length()-1;i++)
					ifprintf(output,"%g,",tmp_eff[i]);
				ifprintf(output,"%g\t",tmp_eff[tmp_eff.Length()-1]);
				
				for(int i=0;i<tmp_pvalue.Length()-1;i++)
					ifprintf(output,"%g,",tmp_pvalue[i]);
				ifprintf(output,"%g\t",tmp_pvalue[tmp_pvalue.Length()-1]);
				
				ifprintf(output,"%g\t%g\t%g\t%g\t%g\t%g\t",average_af,min_af,max_af,chosen_effSize,chosen_cutoff,pvalue);
			}
			else
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str(),average_af,min_af,max_af,chosen_effSize,chosen_cutoff,pvalue);
			if(cond=="")
				ifprintf(output,"\n");
		}
		
		if(cond!="" && condition)
			ifprintf(output,"%g\t%g\t%g\n",chosen_effSize,chosen_cutoff,pvalue);
		
		if(pvalue>1.0)
		pvalue = 1.0;
	}
	if(lower) delete [] lower;
	if(upper) delete [] upper;
	if(mean) delete [] mean;
	return pvalue;
}

void Meta::SKATassoc( GroupFromAnnotation & group )
{
	printf("Performing SKAT ...\n");
	//calculate Q statistics here
	Vector pvalue_SKAT,pos_plot,cond_pvalue_SKAT;
	StringArray chr_plot,geneLabels;
	IFILE output;
	String filename;
	String method = "SKAT_";
	openMetaResultFile( prefix, filename, output, method );
	
	method += MAF_cutoff;
	IFILE reportOutput;
	if(report)
	{
		String reportFile;
		if(prefix =="")
			reportFile = "meta.tophits.SKAT.tbl";
		else if(prefix.Last()=='.' || prefix.Last()=='/')
			reportFile = prefix +  "meta.tophits.SKAT.tbl";
		else
			reportFile = prefix + ".meta.tophits.SKAT.tbl";
		reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
		ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE_DAVIES\tGENE_PVALUE_LIU\tMAF_CUTOFF\tACTUAL_CUTOFF\tVAR\tMAF\tEFFSIZE\tPVALUE\n");
	}
	
	ifprintf(output,"##Method=SKAT\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",total_N);
	if(fullResult)
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tSTATISTICS\tPVALUE_DAVIES\tPVALUE_LIU\n");
	else
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tSTATISTICS\tPVALUE_DAVIES\tPVALUE_LIU\n");
	
	double Qstat,pvalue,pvalue_liu;
	for(int g=0;g<group.annoGroups.Length();g++)
	{
		//      printf("Now working on group %s\n",group.annoGroups[g].c_str());
		if(g%1000==1 && g>1000)
			printf("Finished analyzing %d genes.\n",((int)g/1000)*1000);
		int n = maf[g].Length();
		Vector weight;
		weight.Dimension(n,0);
		if(maf[g].Length()==0)
			continue;
		
		double average_af = maf[g].Average();
		double min_af = maf[g].Min();
		double max_af = maf[g].Max();
		
		String var;
		for(int i=0;i<maf[g].Length()-1;i++)
			var += group.SNPlist[g][i] + ";";
		var += group.SNPlist[g][maf[g].Length()-1];
		
		//get weight based on maf
//		double alpha = 1.0;
//		double beta=25.0;
		String skat_method = "BBeta";
		SetWeight( skat_method, weight,maf[g] );
		Vector tmp,cond_tmp;
		tmp.Dimension(n);
		if(cond!="")
			cond_tmp.Dimension(n);
		for(int i=0;i<n;i++)
		{
			tmp[i] = weight[i]*stats[g][i];
			if(cond!="")
				cond_tmp[i] = weight[i]*cond_stats[g][i];
		}
		Qstat = tmp.InnerProduct(stats[g]);

		double cond_Qstat = _NAN_;
		if(cond!="")
			cond_Qstat = cond_tmp.InnerProduct(cond_stats[g]);
		double * lambda = new double [n];
		CalculateLambda(cov[g],weight,lambda);
		// check lambda for dead loop
		double lambda_sum=0;
		for( int i=0; i<n; i++ )
			lambda_sum += lambda[i];
		if ( lambda_sum < 0.0000000001 ) {
			fprintf(log,"Gene %s lambda sum is zero. Skipped!\n",group.annoGroups[g].c_str());
			continue;
		}

		double Qstat_dav = Qstat;
		double Qstat_liu = Qstat;
		double cond_pvalue=_NAN_,cond_pvalue_liu=_NAN_;
		
		pvalue = MixChidist(lambda, n, Qstat,"Davies");
		
		bool disect_davies=false,disect_liu = false,cond_disect_davies=false,cond_disect_liu=false;
		int disect_itr=0;
		while( (pvalue<=0.0 ||pvalue==2.0 || std::isnan(pvalue)) && disect_itr<10000)
		{
			disect_davies=true;
			Qstat_dav*=0.9999;
			pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
			disect_itr++;
		}
		while((pvalue<=0.0 ||pvalue==2.0 || std::isnan(pvalue)))
		{
			Qstat_dav*=0.99;
			pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
		}
		pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
		disect_itr=0;
		while( (pvalue_liu<=0.0 ||pvalue_liu==2.0 || std::isnan(pvalue_liu)) && disect_itr<10000)
		{
			disect_liu=true;
			Qstat_liu*=0.9999;
			pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
			disect_itr++;
		}
		while((pvalue_liu<=0.0 ||pvalue_liu==2.0 || std::isnan(pvalue_liu)))
		{
			Qstat_liu*=0.99;
			pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
		}
		
		if(cond!="")
		{
			double * lambda = new double [n];
			CalculateLambda(cond_cov[g],weight,lambda);
			Qstat_dav = cond_Qstat;
			Qstat_liu = cond_Qstat;
			
			cond_pvalue = MixChidist(lambda, n, cond_Qstat,"Davies");
			
			int disect_itr=0;
			while( (cond_pvalue<=0.0 ||cond_pvalue==2.0 || std::isnan(cond_pvalue)) && disect_itr<10000)
			{
				cond_disect_davies=true;
				Qstat_dav*=0.9999;
				pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
				disect_itr++;
			}
			while((cond_pvalue<=0.0 ||cond_pvalue==2.0 || std::isnan(cond_pvalue)))
			{
				Qstat_dav*=0.99;
				cond_pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
			}
			cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
			disect_itr=0;
			while( (cond_pvalue_liu<=0.0 ||cond_pvalue_liu==2.0 || std::isnan(cond_pvalue_liu)) && disect_itr<10000)
			{
				cond_disect_liu=true;
				Qstat_liu*=0.9999;
				cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
				disect_itr++;
			}
			while((cond_pvalue_liu<=0.0 ||cond_pvalue_liu==2.0 || std::isnan(cond_pvalue_liu)))
			{
				Qstat_liu*=0.99;
				cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
			}
			if(lambda) delete [] lambda;
		}
		if(std::isnan(pvalue_liu) || std::isnan(pvalue))
		{
			if(fullResult)
			{
				ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());
				
				for(int i=0;i<maf[g].Length()-1;i++)
					ifprintf(output,"%g,",maf[g][i]);
				ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);
				
				for(int i=0;i<singleEffSize[g].Length()-1;i++)
					ifprintf(output,"%g,",singleEffSize[g][i]);
				ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);
				
				for(int i=0;i<singlePvalue[g].Length()-1;i++)
					ifprintf(output,"%g,",singlePvalue[g][i]);
				ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
				if(cond=="")
					ifprintf(output,"%g\t%g\t%g\t%g\t-\t-\n",average_af,min_af,max_af,Qstat);
				else
					ifprintf(output,"%g\t%g\t%g\t%g\t-\t-\t%g\t-\t-\n",average_af,min_af,max_af,Qstat,cond_Qstat);
			}
			else {
				if(cond=="")
					ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t-\t-\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat);
				else
					ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t-\t-\t%g\t-\t-\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat,cond_Qstat);
			}
			continue;
		}
		//tabulate top results
		if(cond!="")
		{
			if((cond_pvalue <report_pvalue_cutoff && cond_pvalue_liu<report_pvalue_cutoff)  && report)
			{
				StringArray variants;
				variants.AddTokens(var,";");
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),cond_disect_davies?"<":"",cond_pvalue,cond_disect_liu?"<":"",cond_pvalue_liu,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
			}
			else
				if((pvalue <report_pvalue_cutoff || pvalue_liu<report_pvalue_cutoff)  && report)
			{
				StringArray variants;
				variants.AddTokens(var,";");
				for(int v=0;v<maf[g].Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
			}
		}
		if(fullResult)
		{
			ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());
			
			for(int i=0;i<maf[g].Length()-1;i++)
				ifprintf(output,"%g,",maf[g][i]);
			ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);
			
			for(int i=0;i<singleEffSize[g].Length()-1;i++)
				ifprintf(output,"%g,",singleEffSize[g][i]);
			ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);
			
			for(int i=0;i<singlePvalue[g].Length()-1;i++)
				ifprintf(output,"%g,",singlePvalue[g][i]);
			ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
			if(cond=="")
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%s%g\n",average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu);
			else
				ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%s%g\t%g\t%s%g\t%s%g\n",average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu,cond_Qstat,cond_disect_davies?"<":"",cond_pvalue,cond_disect_liu?"<":"",cond_pvalue_liu);
		}
		else
		{
			if(cond=="")
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu);
			else
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%s%g\t%g\t%s%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu,cond_Qstat,cond_disect_davies?"<":"",cond_pvalue,cond_disect_liu?"<":"",cond_pvalue_liu);
		}
		
		if(lambda) delete [] lambda;
		
		pvalue_SKAT.Push(pvalue_liu);
		if(cond!="")
			cond_pvalue_SKAT.Push(cond_pvalue_liu);
		StringArray tmp_SNPname;
		tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
		chr_plot.Push(tmp_SNPname[0]);
		pos_plot.Push(tmp_SNPname[1].AsDouble());
		geneLabels.Push(group.annoGroups[g]);
	}
	
	String name = "SKAT (maf<";
	name +=  MAF_cutoff;
	name +=  ")";
	String extraname = "";
	String demo = "";
	double GC = GetGenomicControlFromPvalue(pvalue_SKAT);
	demo="GC=";
	demo += GC;
	writepdf.Draw(pdf,geneLabels,pvalue_SKAT,chr_plot,pos_plot,name,extraname,demo,true);
	if(cond!="")
	{
		name += " conditional analysis";
		double GC = GetGenomicControlFromPvalue(cond_pvalue_SKAT);
		demo="GC=";
		demo += GC;
		writepdf.Draw(pdf,geneLabels,cond_pvalue_SKAT,chr_plot,pos_plot,name,extraname,demo,true);
	}
	ifclose(output);
	if(report)
		ifclose(reportOutput);
	printf("Done.\n\n");
}

/*
void Meta::SKAToptimized( GroupFromAnnotation & group)
{
	double rho;
}
*/

void Meta::CalculateLambda(Matrix & cov,Vector& weight, double * lambda)
{
	int n = cov.rows;
	//calculat sqrt(V)
	Eigen::MatrixXd cov_eigen(n,n);
	for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
		cov_eigen(i,j) = cov[i][j];
	
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_cov(cov_eigen, Eigen::ComputeThinU);
	Eigen::MatrixXd final_eigen(n,n);
	Eigen::MatrixXd final_eigen_rhs(n,n);
	Eigen::MatrixXd final_eigen_lhs(n,n);
	Eigen::VectorXd tmp(n);
	for(int i=0;i<n;i++)
		tmp[i]=sqrt(svd_cov.singularValues()[i]);
	Eigen::VectorXd weight_eigen(n);
	for(int i=0;i<n;i++)
		weight_eigen[i] = weight[i];
	final_eigen_rhs = svd_cov.matrixU()*tmp.asDiagonal()*svd_cov.matrixU().transpose();
	final_eigen_lhs = final_eigen_rhs*weight_eigen.asDiagonal();
	final_eigen = final_eigen_lhs*final_eigen_rhs;
	
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(final_eigen, Eigen::ComputeThinU);
	const Eigen::VectorXd S = svd.singularValues();
	for(int i=0;i<n;i++)
		lambda[i] = fabs(S[i]);
}

// print error to both log and std err
void Meta::ErrorToLog( const char* msg )
{
	fprintf(log, "Error [Meta.cpp]: ");
	fprintf(log, msg);
	fprintf(log, "\n");
	error(msg);
}


/*** pop stratification ***/

// load 1KG MAF
void Meta::load1kgPopMAF()
{
	IFILE popfile = ifopen(popfile_name, "r");
	if (popfile==NULL)
		error("Cannot open 1000G population MAF file: %s\n", popfile_name.c_str());
	bool header_line = 1;
	while(!ifeof(popfile)) {
		String buffer;
		buffer.ReadLine(popfile);
		StringArray tokens;
		tokens.AddTokens(buffer, '\t');
		if (header_line) {
			nPop = tokens.Length() - 8;
			header_line = 0;
			if (nPop<=0)
				error("No population in %s after column 8. Wrong format?\n",popfile_name.c_str());
			continue;
		}
		String markername = tokens[0]+":"+tokens[1]+":"+tokens[3]+":"+tokens[4];
		af1KG[markername].resize(nPop);
		for(int i=0; i<nPop; i++)
			af1KG[markername][i] = tokens[8+i].AsDouble();
	}
	ifclose(popfile);
}

/* this must be divided by N finally
void Meta::updateRegressedTotalAF( String & markername, double total )
{
	int marker_idx = regressedTotalAF.Find( markername );
	if (marker_idx != -1) {
		double prev = regressedTotalAF.Double(markername);
		prev += total;
		regressedTotalAF.SetDouble( markername, prev );
	}
	else
		regressedTotalAF.SetDouble( markername, total );	
}*/


// find the marker in 1KG again
// calculate fk~tilda
double Meta::getAFtilda( String & markername, double raw_af, int study )
{
//	double fk = raw_af - pgamma[study][0];
	double fk = raw_af;
	// match allele first & get fk
	std::map<String, std::vector<double> >::iterator ptr;
	ptr = af1KG.find(markername);
	if ( ptr == af1KG.end() ) { // check if flip
		StringArray tokens;
		tokens.AddTokens( markername, ":" );
		String flip_markername = tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" + tokens[2];
		ptr = af1KG.find(flip_markername);
		if (ptr==af1KG.end()) { // does not exist, use 1/4N
			for(int i=0;i<nPop;i++)
				fk -= pgamma[study][i] * 0.25;
		}
		else { // make a flip vector
			for(int i=0; i<nPop; i++)
				fk -= pgamma[study][i] * (1-(ptr->second[i]));
		}
	}
	else {
		for(int i=0; i<nPop; i++)
			fk -= pgamma[study][i] * (ptr->second[i]);
	}
	return fk;
}


void Meta::addToMapStrVec( std::map<String, std::vector<double> >& variant, int study, String & markername, double fk, int vsize)
{
	std::map<String, std::vector<double> >::iterator p = variant.find(markername);
	if (p==variant.end()) {
		variant[markername].resize( vsize );
		variant[markername][study] = fk;
	}
	else
		variant[markername][study] = fk; 
}

void Meta::SetWeight( String & method, Vector & weight, Vector& maf )
{
	if(method=="burden")// equal weight
		for(int w=0;w<weight.Length();w++)
				weight[w] = 1.0;
	else if(method=="MB") // weight by 1/sqrt( maf*(1-maf) )
		for(int w=0;w<weight.Length();w++)
			weight[w] = sqrt(maf[w]*(1.0-maf[w]));
	else if (method="MAB") {
		for(int w=0;w<weight.Length();w++)
			weight[w] = maf[w];
	}
	else if (method="BBeta") { // truncated beta
		double alpha = 0.5;
		double beta = 0.5;
		double f0 = 2 / Nsamples; // truncate at 4 alleles
		for(int w=0;w<weight.Length();w++) {
			double xmaf = maf[w];
			if (xmaf>0 && xmaf<f0)
				xmaf = f0;
			if (xmaf>0 && xmaf>(1-f0))
				xmaf = 1-f0;
			double beta_density = GetBetaDensity(alpha,beta,xmaf);
			weight[w] = (beta_density*beta_density);
		}
	}
	else
		error("Invalid weight %s!\n",method.c_str());
}

