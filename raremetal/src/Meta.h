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
        Meta( FILE * plog );
        ~Meta();

    //Input/Output options  
	static String summaryFiles;
	static String covFiles;
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
	static bool report;
	static bool fullResult;
	static bool founderAF;
	static bool dosage;
	static double CALLRATE;
	static double report_pvalue_cutoff;
	static double MAF_cutoff;
	static double HWE;
	static int marker_col;
	static int cov_col;
	static bool altMAF; // if TRUE, exclude size of studies that do not contain that variant
	static bool GeneOnly; // only perform gene based test
	static bool RegionStatus; // if TRUE, restrict gene-based test to the specified region
	static String Region; // raw region option
	static String Chr;
	static int Start;
	static int End; // 3 variables to define region
	FILE * log;

	static bool useExactMetaMethod;

	//saved single variant information from pooled studies
	StringArray scorefile;
	StringArray covfile;
 	String pdf_filename;
 	StringIntHash hashToRemoveDuplicates;
	//flipSNP has all SNPs fliped saved. Key is study:chr:pos, value is an integer
	//SNPexclude has all SNPs with non-consistent ref/alt alleles saved. key is study:chr:pos; value is an integer
	StringIntHash flipSNP;
	StringIntHash SNPexclude;
	StringIntHash SNPmaf;
	//SNPstat has the single variants stats for all SNPs 
	StringDoubleHash singleVarPvalue;
	StringDoubleHash singleVarEff;
	StringDoubleHash SNPstat;
	StringDoubleHash SNP_Vstat;
	StringDoubleHash SNPstat_cond;
	StringDoubleHash SNP_Vstat_cond;
	Vector SNPmaf_maf;
	StringArray SNPmaf_name;
	IntArray SNP_effect_N;
	int flipCount;
	int skip_count;

	// for pooling summary stats
	StringIntHash usefulSize; // pooled N info
	StringDoubleHash usefulAC; // pooled AC info
	StringIntHash recSize; // #samples to include (--altMAF option only)
	StringArray refalt; // study:chr:pos as key and ref:alt as value
	StringIntHash directionByChrPos; // chr:pos as key and save the positions in directions array.
	StringArray directions; // stores the actual direction strings
	String dupcheck_chr;
	int varCount;
	Vector pvalueAll,pvalue1,pvalue5,pvalueAll_cond,pvalue1_cond,pvalue5_cond; // store some p value for plot
	// for annotation
	String target_chr;
	int target_pos;
	int target;
	double target_pvalue;

	// for exact method
	Vector Ydelta; // yk - ymean
	Vector Ysigma; // sigmak(variance), divided by nk-1
	StringDoubleHash V2; // sum(4*nk*fk^2)
	StringDoubleHash V3; // sum(2*nk*fk)
	StringDoubleHash residual_adj; // sigma outside

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

	Vector * maf;
	Vector * stats;
	Vector * cond_stats;
	Vector * singlePvalue;
	Vector * singleEffSize;
	Matrix * cov;
	Matrix * cond_cov;

	IntArray SampleSize;
	int total_N;
	Vector GCbyStudy;

	StringArray chr_plot,geneLabel;
	Vector pos_plot,chisq_before_GC;

	PDF pdf;
	WritePDF writepdf;

	void Prepare();
	void PoolSummaryStat(GroupFromAnnotation & group);
	void Run(GroupFromAnnotation & group);
	void CalculateGenotypeCov(SummaryFileReader & covReader,String chr, int pos,int study,Vector & result);
	double GrabGenotypeCov(SummaryFileReader & covReader,int study,String chr1,String pos1,String chr2,String pos2,String & SNP1, String & SNP2);
	void CalculateXXCov(int study,Matrix & result);

	void setMetaStatics();
	void openMetaFiles();
	void prepareConditionalAnalysis();
	bool setCondMarkers();

	bool updateYstat( int study );
	void UpdateDirection(int & direction_idx,int study,char marker_direction,String & chr_pos,bool exclude);
	void UpdateExcludedMarker(int & study, String & chr_pos, int filter,String markername);
	void UpdateStrIntHash(String & chr_pos, int val, StringIntHash & sihash);
	void UpdateACInfo(String & chr_pos,double AC);
	void UpdateStats(String & markerName,double stat,double vstat,bool flip);
	void updateExactStats( String & markerName, int study, int current_AC, int current_N );
	char GetDirection(String & chr_pos,double effsize,bool flip);
	int MatchTwoAlleles(String refalt_current,int & marker_idx,String & chr_pos);
	int MatchOneAllele(String refalt_current, int & marker_idx);
	
	bool poolSingleRecord( int study, double & current_chisq, int & duplicateSNP, bool adjust, String & buffer, SummaryFileReader & covReader );

	bool isDupMarker( String & chr_str, String & chr_pos );
	void setRefAltHashKey( String & refalt_current, StringArray & tokens, int c1, int c2 );
	void setPolyMatch( int & match, String & chr_pos, String & refalt_current, StringArray & tokens, int marker_idx );
	void updateSNPcond( int study, bool flip, int adjust, String & chr_pos, StringArray & tokens, SummaryFileReader & covReader );
	void setPooledAF();
	void printSingleMetaHeader( String & filename, IFILE & output );
	void printOutVcfHeader( String & vcf_filename, IFILE & vcfout );
	void printSingleMetaVariant( GroupFromAnnotation & group, int i, IFILE & output, IFILE & vcfout );
	void annotateSingleVariantToGene( GroupFromAnnotation & group, double pvalue, double cond_pvalue, StringArray & tmp );
	void plotSingleMetaGC( IFILE & output, bool calc_gc );

	void loadSingleStatsInGroup( GroupFromAnnotation & group );
	void loadSingleCovInGroup( GroupFromAnnotation & group );
	void BurdenAssoc(String method,GroupFromAnnotation & group,Vector *& maf,Vector *& stats,Vector *& con_stats,Matrix *& cov,Matrix *& cond_cov,Vector *& singleEffSize,Vector *& singlePvalue);
	void VTassoc( GroupFromAnnotation & group );
	double VTassocSingle(GroupFromAnnotation & group, Vector & maf_cutoff, IFILE reportOutput, IFILE output, int & g,bool condition,String & method);
	void SKATassoc( GroupFromAnnotation & group );

	void CalculateLambda(Matrix & cov, double * weight, double * lambda);
	void ErrorToLog(const char * msg);
};

#endif
