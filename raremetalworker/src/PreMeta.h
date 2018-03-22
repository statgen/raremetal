////////////////////////////////////////////////////////////////////// 
// PreMeta.cpp 
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

#ifndef __PREMETA_H__
#define __PREMETA_H__

#include "FastFit.h"
#include "InputFile.h"
#include "GroupFromAnnotation.h"
#include "DataQC.h"
#include <savvy/reader.hpp>
#include <vector>
#include <map>
#include <string>

class PreMeta
{
public:
    PreMeta(FILE *logfile, IFILE &SCOREoutput1, IFILE &SCOREoutput_rec1, IFILE &SCOREoutput_dom1, IFILE &SCOREcov1,
            IFILE &SCOREcov_rec1, IFILE &SCOREcov_dom1);

    ~PreMeta();

    FILE *log;

    static bool vcfAnnotated;
    /**
     * Whether to automatically compress the output files within this range
     */
    static bool zipOutput;
    static bool checkRef;

    //Input/Output options
    static String vcfInput;
    /**
     * A filename prefix to be prepended to all output file names. See:
     * See https://genome.sph.umich.edu/wiki/RAREMETALWORKER#OUTPUT_FILE_FORMATS
     */
    static String outputFile;
    /**
     * --xStart takes an integer that described the start position of nonPAR region on chromosome X.
     * The default is 2699520 based on Human Genome build 19.
     *
     * --xStart takes an integer that described the end position of nonPAR region on chromosome X.
     * The default is 154931044 based on Human Genome build 19.
     */
    static int Xstart, Xend;
    /**
     * The string used to label chromosome X in input files. Defaults to "X"
     */
    static String xLabel;
    /**
     * Whether to read dosage information from the VCF file (in the DS field)
     */
    static bool dosage;
    static bool genoFromPed;
    static bool genoFromVCF;
    static bool FounderFreq;
    static bool correctGC;
    /**
     * RAREMETALWORKER generates LD matrices between a current marker that it is working on and all markers within this
     *  range. Default is 1e6 bases.
     */
    static int window;
    /**
     * The label used for males in PED files. Defaults to "1"
     */
    static int maleLabel;
    /**
     * The label used for females in PED files. Defaults to "2"
     */
    static int femaleLabel;
    /**
     * If --recessive is used, then RAREMETALWORKER generates recessive results in addition to the additive results.
     * A separate pdf file with QQ and Manhattan plots based on recessive results is generated with name
     *   yourprefix.traitname.recessive.plots.pdf.
     */
    static bool recessive;
    /**
     * TODO: should this say Dominant? (pasted from wiki)
     * If --dominant is used, then RAREMETALWORKER generates recessive results in addition to the additive results.
     *   A separate pdf file with QQ and Manhattan plots based on recessive results is generated with
     *   name yourprefix.traitname.dominant.plots.pdf.
     */
    static bool dominant;
    static bool additive;
    static bool calculateOR;
    static bool Simplify;
    static String Region;
    static String varListName;
    static bool splitUV;
    static bool newFormat;
    static bool printCaseAC;
    static bool debug; // if run debug info

    int warnings;
    int numFounders;
    int malehwewarning;

    savvy::indexed_reader reader;
    savvy::variant<std::vector<float>> record;
    // indicate a sample is case or control. used in printCaseAC
    std::map<int, bool> sampleCaseIndicator; // sample id in vcf -> if case, true
    int caseAC = 0;
    int controlAC = 0;

    //this is the position of each marker in a gene for output
    int pos, hom1, het, hom2;
    double founderMaf, markerMaf, callRate;
    String chr, refAllele, rareAllele;
    double mac, hwe_pvalue;
    Eigen::VectorXf genotype, genotype_rec, genotype_dom;
    Eigen::VectorXf transGeno, transGeno_rec, transGeno_dom;
    Vector residuals;
    Vector chisq_before_GCcorrect_rec, chisq_before_GCcorrect_dom, chisq_before_GCcorrect;
    Vector pvalueAll, pvalueAll_rec, pvalueAll_dom;
    Vector pos_plot, pos_plot_rec, pos_plot_dom;
    StringArray chr_plot, chr_plot_rec, chr_plot_dom;
    Vector maf_GC, maf_GC_rec, maf_GC_dom;
    IFILE SCOREoutput, SCOREoutput_rec, SCOREoutput_dom;
    IFILE SCOREcov, SCOREcov_rec, SCOREcov_dom;
    Matrix genotypeAll, genotypeAll_rec, genotypeAll_dom; //this matrix has the all marker genotype within a window
    IntArray markerInLD;

    Vector chisqred;
    double lambda;
    Matrix projectionMat, inv;
    std::map<std::string, std::vector<int> > varList; // variant list if needed

    bool SetCaseControlMap(String &pedname);

    void CalculateProjectionMatrix(FastTransform &trans, FastFit &engine, Vector &sigma2);

    void Run(String &pedname, Pedigree &ped, FastTransform &trans, FastFit &engine, GroupFromAnnotation &group,
             SanityCheck &checkData, KinshipEmp &kin_emp);

    void
    runGenoFromPed(Pedigree &ped, FastTransform &trans, FastFit &engine, SanityCheck &checkData, KinshipEmp &kin_emp,
                   Vector &sigma2);

    void
    runGenoFromVcf(Pedigree &ped, FastTransform &trans, FastFit &engine, SanityCheck &checkData, KinshipEmp &kin_emp,
                   Vector &sigma2);

    void smallSanityCheck(double averageAF, int marker_count);


    void VCFSanityCheck();

    bool GetGenotypeVectorVCFDosage(FastTransform &trans, Pedigree &ped);

    bool GetGenotypeVectorVCF(FastTransform &trans, Pedigree &ped, SanityCheck &checkData, FILE *log);

    int GetGenotypeVectorPED(FastTransform &trans, Pedigree &ped, int markerid, FILE *log);

    //bool GetGenotypeVector(FastTransform & trans, Pedigree & ped,int markerid);
    void PrintScoreFileHeader(IFILE &output, FastFit &engine, FastTransform &trans, Pedigree &ped);

    void PrintCovFileHeader(IFILE &file);

    void PrintToCovOutput(IFILE &file, Matrix &gt, String &chr, int pos, Pedigree &ped, FastTransform &trans,
                          FastFit &engine, Vector &sigma2);

    double
    CalculateCovWithSigma(FastFit &engine, FastTransform &trans, Pedigree &ped, Matrix &genotypeAl, Vector &sigma2,
                          int m);

    double CalculatePlainCov(Matrix &gt, Vector &sigma2, int m, FastTransform &trans);

    void
    CalculateUnrelatedAssocStats(double &effSize, double &pvalue, double &chisq, double &numerator, double &denominator,
                                 Eigen::VectorXf &genotype, FastFit &engine, FastTransform &trans, Vector &sigma2);

    void CalculateAssocStats(double &effSize, double &pvalue, double &numerator, double &denominator, double &chisq,
                             Eigen::VectorXf &transGeno, FastFit &engine, FastTransform &trans, Vector &sigma2);

    void RelatedAssoc(IFILE SCOREoutput, IFILE SCOREoutput_rec, IFILE SCOREoutput_dom, Pedigree &ped, FastFit &engine,
                      FastTransform &trans, Vector &sigma2);

    void UnrelatedAssoc(Pedigree &ped, FastFit &engine, FastTransform &trans, Vector &sigma2);

    void GeneratePlots(String &filename, Vector &pvalueAll, StringArray &chr_plot, Vector &pos_plot,
                       GroupFromAnnotation &group, IFILE &SCOREoutput, Vector &chisq_before_GCcorrect, Vector &maf_GC,
                       String &model);

    void GetRealPrefix(String &file);

    void setVarList();

    void updateScoreVar(std::string &chr_str, int &current_score_var, int &current_score_index);

    void updateCovVar(std::string &chr_str, int &current_cov_var, int &current_cov_index);
};

#endif
