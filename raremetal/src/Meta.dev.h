#ifndef __INITIAL_H__
#define __INITIAL_H__

#include "VcfRecord.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include "GroupFromAnnotation.h"
#include "MathMatrix.h"
#include "StringHash.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "SummaryFileReader.h"
#include "WritePDF.h"

class Meta
{
    public:
        Meta();
        ~Meta();

        //Input/Output options  
	static String studyName;
	static String prefix;
	static String cond;
	static bool correctGC;
	static bool outvcf;
	static bool tabix;
	static bool Burden;
	static bool MB;
	static bool SKAT;
	static bool SKATO;
	static bool VTa;
	static bool VTp;
	static bool report;
	static bool fullResult;
	static bool founderAF;
	static double CALLRATE;
	static double  report_pvalue_cutoff;
	static double MAF_cutoff;
	static double HWE;
	static int marker_col;
	static int cov_col;

	//saved single variant information from pooled studies
	StringArray studies;
   String pdf_filename;
	//flipSNP has all SNPs fliped saved. Key is study:chr:pos, value is an integer
	//SNPexclude has all SNPs with non-consistent ref/alt alleles saved. key is study:chr:pos; value is an integer
	StringIntHash flipSNP,SNPexclude,SNPmaf;
	StringDoubleHash singleVarPvalue,singleVarEff;
	StringDoubleHash SNPstat,SNP_Vstat,SNPstat_cond,SNP_Vstat_cond; //SNPstat has the single variants stats for all SNPs 
Vector SNPmaf_maf;
StringArray SNPmaf_name;

	//these are for conditional analysis
	StringArray commonVar,common_chr,common_ref,common_alt;
StringIntHash conditionVar;
	IntArray common_pos,FormatAdjust;
	Matrix  * XX_inv;
	IntArray * commonVar_study;
	Vector * commonVar_effSize;
	Vector * commonVar_V;
	Vector * commonVar_U;
	Vector * commonVar_betaHat;
	IntArray ** commonVar_markers_in_window;
	Vector ** commonVar_marker_cov;
	bool * cond_status;

	IntArray SampleSize;
	int total_N;

	PDF pdf;
	WritePDF writepdf;

	void Run(GroupFromAnnotation & group);
	void PoolSummaryStat(GroupFromAnnotation & group);
	void Prepare();
	void BurdenAssoc(String method,GroupFromAnnotation & group,Vector *& maf,Vector *& stats,Vector *& con_stats,Matrix *& cov,Matrix *& cond_cov,Vector *& singleEffSize,Vector *& singlePvalue);
	void CalculateXXCov(int study,Matrix & result);
	void CalculateGenotypeCov(SummaryFileReader & covReader,String chr, int pos,int study,Vector & result);
	double GrabGenotypeCov(SummaryFileReader & covReader,int study,String chr1,String pos1,String chr2,String pos2,String & SNP1, String & SNP2);
	double GetBetaDensity(double a, double b, double x);
	double CalculateCorrCoef(Vector & a,Vector & b);
	void RevertAllele(String SNP, String & newSNP);
void CalculateLambda(Matrix & cov, Vector & weight, double * lambda,bool optimal, double rho);
void CalculateMinPvalue(Matrix &cov, Vector & weight, double *lambda,double *rho,double & Qskat, double & Qburden)
};

#endif
