#ifndef __INITIAL_H__
#define __INITIAL_H__


#include "GroupFromAnnotation.h"
#include "MathMatrix.h"
#include "StringHash.h"

#define MATHLIB_STANDALONE

#include <Rmath.h>
#include "SummaryFileReader.h"
#include "WritePDF.h"

#include <map>

const double LN_10 = log(10);

class Meta
{
public:
    Meta();
    ~Meta();

    //Input/Output options  
    String summaryFiles;
    String covFiles;
    String prefix;
    String cond;
    bool correctGC;
    bool outvcf;
    bool tabix;
    bool Burden;
    bool MB;
    bool MAB;
    bool BBeta;
    bool SKAT;
    bool VTa;
    bool report;
    bool fullResult;
    bool founderAF;
    bool dosage;
    double CALLRATE;
    double report_pvalue_cutoff;
    double MAF_cutoff;
    double HWE;
    int marker_col;
    int cov_col;
    bool altMAF; // if TRUE, exclude size of studies that do not contain that variant
    bool RegionStatus; // if TRUE, restrict gene-based test to the specified region
    String Region; // raw region option
    String Chr;
    int Start;
    int End; // 3 variables to define region
    FILE *log;
    bool skipOutput;

    bool useExactMetaMethod;
    bool normPop;
    String popfile_name;
    bool simplifyCovLoad;
    bool relateBinary;
    bool debug;
    bool matchOnly;
    bool matchByAbs;
    double matchDist;
    double minMatchMAF; // MAF threshold for adjustment
    double maxMatchMAF;
    String dosageOptionFile;
    bool sumCaseAC;
    bool bHeterogeneity;
    bool logP; // write p-values in -log10 scale

    //saved single variant information from pooled studies
    StringArray scorefile;
    StringArray covfile;
    std::vector<bool> dosageOptions;
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
    StringDoubleHash SNP_heterog_stat;
    StringIntHash SNP_heterog_df;
    StringDoubleHash SNP_heterog_cond_stat;
    StringIntHash SNP_heterog_cond_df;
    StringIntHash groupAnchor;
    Vector SNPmaf_maf;
    StringArray SNPmaf_name;
    IntArray SNP_effect_N;
    int Nsamples;
    int flipCount;
    int skip_count;
    StringIntHash caseAC;
    StringIntHash controlAC;

    // for pooling summary stats
    StringIntHash usefulSize; // pooled N info
    StringDoubleHash usefulAC; // pooled AC info
    StringIntHash recSize; // #samples to include (--altMAF option only)
    StringArray refalt; // study:chr:pos as key and ref:alt as value
    StringIntHash directionByChrPos; // chr:pos as key and save the positions in directions array.
    StringArray directions; // stores the actual direction strings
    String dupcheck_chr;
    int varCount;
    Vector pvalueAll, pvalue1, pvalue5, pvalueAll_cond, pvalue1_cond, pvalue5_cond; // store some p value for plot
    // for annotation
    String target_chr;
    int target_pos;
    int target;
    double target_pvalue;

// for cov load
    StringIntHash markerPosHash;
    Vector markersExp; // for new format
    StringArray markersInWindow; // for old format
    StringArray markersCov;

// /net/fantasia/home/yjingj/METAL/1KG/MAF_1KG.txt
    // for exact method
    Vector Ydelta; // yk - ymean
    Vector Ysigma2; // sigmak(variance), divided by nk-1
    std::map<String, std::vector<double> > variant_fk; //fk of each variant for each study. will turn into fk-tilda if --normPop. last one is total maf
    std::map<String, std::vector<double> > variant_nk; // nk of each variant for each study. last one is new
//	StringDoubleHash V2; // sum(4*nk*fk^2)
//	StringDoubleHash residual_adj; // sigma tilda outside	
//	StringDoubleHash NkDeltak; // nk * deltak
//	StringDoubleHash regressedTotalAF; // for pop stratification, update f-tilda
    std::map<String, std::vector<double> > af1KG; // 1000g pop af for all markers...
    Matrix pgamma; // coefficients for population in each study
    int nPop; // number of populations

    Vector Const_binary; // Dajiang's C for binary
    Vector Ysigma2g;
    Vector Ysigma2e;
    Vector Ymean;

    //these are for conditional analysis
    StringArray commonVar, common_chr, common_ref, common_alt;
    StringIntHash conditionVar;
    IntArray common_pos, FormatAdjust;
    Matrix *XX_inv;
    IntArray *commonVar_study;
    Vector *commonVar_effSize;
    Vector *commonVar_V;
    Vector *commonVar_U;
    Vector *commonVar_betaHat;
    IntArray **commonVar_markers_in_window;
    Vector **commonVar_marker_cov;
    bool *cond_status;

    Vector *maf;
    Vector *stats;
    Vector *cond_stats;
    Vector *singlePvalue;
    Vector *singleEffSize;
    Matrix *cov;
    Matrix *cond_cov;

    IntArray SampleSize;
    int total_N;
    Vector GCbyStudy;

    StringArray chr_plot, geneLabel;
    Vector pos_plot, chisq_before_GC;

    PDF pdf;
    WritePDF writepdf;

    void Prepare();

    void PoolSummaryStat(GroupFromAnnotation &group);

    void Run(GroupFromAnnotation &group);

    void CalculateGenotypeCov(SummaryFileReader &covReader, String chr, int pos, int study, Vector &result);

    double GrabGenotypeCov(SummaryFileReader &covReader, int study, String chr1, String pos1, String chr2, String pos2,
                           String &SNP1, String &SNP2);

    void CalculateXXCov(int study, Matrix &result);

    void parseScoreFiles();
    void parseCovFiles();

    void setLogFile(FILE *plog = nullptr);

    void prepareConditionalAnalysis();

    bool setCondMarkers();

    bool updateYstat(int study);

    void UpdateDirection(int &direction_idx, int study, char marker_direction, String &chr_pos, bool exclude);

    void UpdateExcludedMarker(int &study, String &chr_pos, int filter, String markername);

    void UpdateStrIntHash(String &chr_pos, int val, StringIntHash &sihash);

    void UpdateStrDoubleHash(String &chr_pos, double val, StringDoubleHash &sdhash);

    void UpdateACInfo(String &chr_pos, double AC);

    void UpdateStats(int study, String &markerName, double stat, double vstat, bool flip);
    void UpdateHetStats(int study, String &markerName, double stat, double vstat, bool flip);

    char GetDirection(String &chr_pos, double effsize, bool flip);

    int MatchTwoAlleles(String refalt_current, int &marker_idx, String &chr_pos);

    int MatchOneAllele(String refalt_current, int &marker_idx);

    bool poolSingleRecord(int study, double &current_chisq, int &duplicateSNP, bool adjust, String &buffer,
                          SummaryFileReader &covReader);

    void poolHeterogeneity(int study, bool adjust, String &buffer, SummaryFileReader &covReader);

    bool isDupMarker(String &chr_str, String &chr_pos);

    void setRefAltHashKey(String &refalt_current, StringArray &tokens, int c1, int c2);

    void setPolyMatch(int &match, String &chr_pos, String &refalt_current, StringArray &tokens, int marker_idx);

    void updateSNPcond(int study, bool flip, int adjust, String &chr_pos, StringArray &tokens, SummaryFileReader &covReader);
    void UpdateHetCondStats(int study, bool flip, int adjust, String &chr_pos, StringArray &tokens, SummaryFileReader &covReader);

    void setPooledAF();

    void printSingleMetaHeader(String &filename, IFILE &output);

    void printOutVcfHeader(String &vcf_filename, IFILE &vcfout);

    void printSingleMetaVariant(GroupFromAnnotation &group, int i, IFILE &output, IFILE &vcfout);
    void WriteSingleVariantResults(GroupFromAnnotation &group);

    void annotateSingleVariantToGene(GroupFromAnnotation &group, double pvalue, double cond_pvalue, StringArray &tmp);

    void plotSingleMetaGC(IFILE &output, bool calc_gc);

    void loadSingleStatsInGroup(GroupFromAnnotation &group);

    void loadSingleCovInGroup(GroupFromAnnotation &group);

    void BurdenAssoc(String method, GroupFromAnnotation &group, Vector *&maf, Vector *&stats, Vector *&con_stats,
                     Matrix *&cov, Matrix *&cond_cov, Vector *&singleEffSize, Vector *&singlePvalue);

    void VTassoc(GroupFromAnnotation &group);

    double VTassocSingle(GroupFromAnnotation &group, Vector &maf_cutoff, IFILE reportOutput, IFILE output, int &g,
                         bool condition, String &method);

    void SKATassoc(GroupFromAnnotation &group);

    void CalculateLambda(Matrix &cov, Vector &weight, double *lambda);

    // pop correction
    void
    addToMapStrVec(std::map<String, std::vector<double> > &variant, int study, String &markername, double fk, int vsize,
                   double initial_number);

    void SetWeight(String &method, Vector &weight, Vector &maf);

    void FitPgamma();

    void fitpGammaForSingleStudy(int study);

    bool add1kgMafToX(Matrix &X, String &markername, int current_index);

    void setGammaFromRegression(int study, Matrix &X, Vector &Y);

    void load1kgPopMAF();

    void updateRegressedTotalAF(String &markername, double total);

    double getAFtilda(String &markername, double raw_af, int study);

    void addToMapStrVec(std::map<String, Vector> &variant, int study, String &markername, double fk);

    void ErrorToLog(const char *msg);

    // old & new cov format
    void readCovOldFormatLine(int study, StringArray &tokens, int &m);

    void readCovNewFormatLine(int study, StringArray &tokens, int &m);

    void addNewFormatCov(int mexp, String &cov_str, Vector &covs);

    void updateExactCov(int study, int m, int s, StringArray &chr, StringArray &pos, Matrix &cov_i);

    // cov matrix
    void updateGroupStats(GroupFromAnnotation &group, int study, int g, bool newFormat);

    void updateSingleVariantGroupStats(GroupFromAnnotation &group, int study, int g, Matrix &cov_i, Matrix &GX,
                                       StringArray &chr, StringArray &pos, int m, int gvar_count, bool newFormat);

    void updateSingleVariantGroupStatsOldFormat(GroupFromAnnotation &group, int study, int g, Matrix &cov_i, Matrix &GX,
                                                StringArray &chr, StringArray &pos, int loc, int m, int gvar_count,
                                                double multiplyFactor);

    void updateSingleVariantGroupStatsNewFormat(GroupFromAnnotation &group, int study, int g, Matrix &cov_i, Matrix &GX,
                                                StringArray &chr, StringArray &pos, int loc, int m, int gvar_count,
                                                double multiplyFactor);
};

struct SingleVariantResult {
  double u, v, n;
  double chisq;
  double pvalue;
  double log_pvalue; // -log10(pvalue)
  double effSize;
  bool disect;
  double h2;
  double effSize_se;

  SingleVariantResult() {}

  SingleVariantResult(double u, double v, double n) {
    setStats(u, v, n);
    calculate();
  }

  void setStats(double u, double v, double n) {
    this->u = u;
    this->v = v;
    this->n = n;
  }

  void calculate() {
    chisq = u * u / v;
    pvalue = pchisq(chisq, 1, 0, 0);
    log_pvalue = -pchisq(chisq, 1, 0, 1) / LN_10;
    effSize = u / v;
    h2 = v * effSize * effSize / n;
    effSize_se = 1.0 / sqrt(v);

    while (pvalue == 0.0) {
      disect = true;
      chisq *= 0.999;
      pvalue = pchisq(chisq, 1, 0, 0);
    }
  }
};

#endif
