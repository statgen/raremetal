/*
 * Perform all calculation steps required for meta-analysis according to one or more selected methods.
 */

#include "StringArray.h"
#include "Error.h"
#include <stdio.h>
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
#include <math.h> // pow, abs

Meta::Meta() {
  summaryFiles = "";
  covFiles = "";
  popfile_name = "";
  prefix = "";
  correctGC = false;
  Burden = false;
  MB = false;
  MAB = false;
  BBeta = false;
  SKAT = false;
  VTa = false;
  outvcf = false;
  fullResult = false;
  report = false;
   report_pvalue_cutoff = 1e-06;
  founderAF = false;
  HWE = 1e-05;
  CALLRATE = 0.95;
  MAF_cutoff = 0.05;
  marker_col = 2;
  cov_col = 3;
  altMAF = false;
  tabix = false;
  dosage = false;
  RegionStatus = false;
  Region = "";
  cond = "";
  Chr = "";
  Start = -1;
  End = -1;
  useExactMetaMethod = false; // use Jingjing's exact method
  normPop = false; // correct for population stratification
  simplifyCovLoad = false; // only load cov matrix between selected marker for group test
  relateBinary = false;
  debug = false; // show debug info & debug output
  matchOnly = false; // only use matched SNP in 1000G
  matchByAbs = false; // use min abs to calculate dist
  matchDist = 0; // in --normPop, exclude variants with dist > matchDist
  minMatchMAF = 0; // MAF threshold for adjustment
  maxMatchMAF = 1;
  dosageOptionFile = "";
  sumCaseAC = false;
  bHeterogeneity = false;
  log = nullptr;
  skipOutput = false;
  Nsamples = -1;
}

Meta::~Meta() {}

void Meta::setLogFile(FILE *plog) {
  if (plog == nullptr) {
    String filename;
    if (prefix == "") {
      filename = "raremetal.log";
    }
    else if (prefix.Last() == '.' || prefix.Last() == '/') {
      filename = prefix + "raremetal.log";
    }
    else {
      filename = prefix + ".raremetal.log";
    }

    log = freopen(filename, "wt", stderr);
  }
  else {
    log = plog;
  }
}

//This function read all study names from a file
//and save the names in the StringArray files
void Meta::Prepare() {
    // Setup PDF
    if (prefix.Last() == '.' || prefix.Last() == '/') {
      pdf_filename = prefix + "meta.plots.pdf";
    }
    else {
      pdf_filename = prefix + ".meta.plots.pdf";
    }
    pdf.OpenFile(pdf_filename);

    // Parse files with paths to score statistics and covariance matrices per study
    parseScoreFiles();
    parseCovFiles();

    // Restrict to only a single region?
    if (Region != "") {
      RegionStatus = true;
      printf("Restrict analysis to region %s!\n", Region.c_str());
      StringArray tf;
      tf.AddTokens(Region, ":-");
      Chr = tf[0];
      Start = tf[1].AsInteger();
      End = tf[2].AsInteger();
    }

    SampleSize.Dimension(scorefile.Length());
    for (int s = 0; s < scorefile.Length(); s++)
    {
        SampleSize[s] = -1;
    }

    if (dosageOptionFile != "")
    { // read dosage
        IFILE file = ifopen(dosageOptionFile, "r");
        while (!ifeof(file))
        {
            String buffer;
            buffer.ReadLine(file);
            bool status;
            if (buffer == "dosage")
            {
                status = true;
            } else if (buffer == "genotype")
            {
                status = false;
            } else
            {
                error("Invalid option in dosageOptionFile! Allowed option is dosage or genotype!\n");
            }
            dosageOptions.push_back(status);
        }
        ifclose(file);
    }
    if (dosageOptionFile != "" && dosageOptions.size() != scorefile.Length())
    {
        error("DosageOptionFile should have the same line number as score files. Please check your input files!\n");
    }

    if (normPop)
    {
        FitPgamma();
    }

    //if conditional analysis says yes
    if (cond != "")
    {
        prepareConditionalAnalysis();
    }
}

/**
 * fit pgamma based on af in each study
 * need to read through all score files
 */
void Meta::FitPgamma()
{
    load1kgPopMAF();
    pgamma.Dimension(scorefile.Length(), nPop);
    for (int study = 0; study < scorefile.Length(); study++)
    {
        fitpGammaForSingleStudy(study);
    }
    if (debug)
    {
        String filename = prefix + ".debug.maf-gammas.txt";
        IFILE file = ifopen(filename, "w", InputFile::UNCOMPRESSED);
        ifprintf(file, "#STUDY Intercept Gamma\n");
        for (int study = 0; study < scorefile.Length(); study++)
        {
            ifprintf(file, "%d", study);
            for (int i = 0; i < nPop; i++)
            {
                ifprintf(file, " %g", pgamma[study][i]);
            }
            ifprintf(file, "\n");
        }
        ifclose(file);
    }
}

/**
 * Read data from the file associated with a specified study, and process for use
 *  1. Open and read the appropriate data file
 *  2. Discard lines that do not meet certain filter criteria (cutoff parameters)
 *  3. Apply dosage information (if specified)
 *  4. Calculate distances and population stratification corrections based on 1000 G (if appropriate)
 * @param study The number of the study row from the input scores file
 */
void Meta::fitpGammaForSingleStudy(int study)
{
    printf("Doing population correction for study %d...\n", study);
    Matrix X; // 1kg mafs
    Vector Y; // sample maf
//	int Nmarkers = 2000;
//	X.Dimension( Nmarkers, nPop );
    int current_index = 0;

    String match_id_file_name = prefix + ".debug.match_id_study";
    match_id_file_name += study;
    match_id_file_name += ".txt";
    IFILE match_id_file;
    if (debug)
    {
        match_id_file = ifopen(match_id_file_name, "w", InputFile::UNCOMPRESSED);
        ifprintf(match_id_file, "#ID MAF Match1KG_MAF\n");
    }

    /** load mafs from sites **/
    String filename = scorefile[study];
    IFILE file;
    file = ifopen(filename, "r");
    if (file == NULL)
    {
        error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.summary.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.summary.file.gz\n",
              scorefile[study].c_str());
    }
    bool first_line = 1;
    bool adjust;
    int dist_skipped = 0;
    int matchonly_skipped = 0;
    int match_maf_skipped = 0;

    bool use_dosage = false;
    if (dosageOptionFile != "")
    {
        use_dosage = dosageOptions[study];
    } else
    {
        use_dosage = dosage;
    }

    while (!ifeof(file))
    {
        String buffer;
        buffer.ReadLine(file);
        if (first_line)
        {
            tellRvOrRmw(buffer, adjust, marker_col, cov_col);
            first_line = 0;
            continue;
        }
        if (buffer.FindChar('#') != -1)
        {
            continue;
        }
        StringArray tokens;
        tokens.AddTokens(buffer, SEPARATORS);
        if (tokens[0].Find("chr") != -1)
        {
            tokens[0] = tokens[0].SubStr(3);
        }

        String chr_pos = tokens[0] + ":" + tokens[1];
        String markername = tokens[0] + ":" + tokens[1] + ":" + tokens[2] + ":" + tokens[3];
        bool fail = isDupMarker(tokens[0], chr_pos);
        if (fail)
        {
            continue;
        }
        int c1, c2, c3;
//printf("%s\n",markername.c_str());
        if (use_dosage)
        {
            c3 = tokens[4].AsDouble() * tokens[6].AsDouble() * tokens[6].AsDouble();
            c2 = tokens[4].AsDouble() * 2.0 * tokens[6].AsDouble() * (1.0 - tokens[6].AsDouble());
            c1 = tokens[4].AsDouble() * (1.0 - tokens[6].AsDouble()) * (1.0 - tokens[6].AsDouble());
        } else
        {
            c1 = tokens[10 - adjust].AsDouble();
            c2 = tokens[11 - adjust].AsDouble();
            c3 = tokens[12 - adjust].AsDouble();
        }
        if ((tokens[2] == "." && tokens[3] == ".") || (tokens[2] == tokens[3] && c1 + c2 != 0 && c2 + c3 != 0))
        {
            continue;
        }
        if (tokens[8 - adjust].AsDouble() < CALLRATE || tokens[9 - adjust].AsDouble() < HWE)
        {
            continue;
        }
        double current_AC;
        int current_N;
        current_N = c1 + c2 + c3;
        current_AC = 2 * c3 + c2;
//		if (current_AC==0)
//			continue;

        /** save to maf matrix **/
        if (X.rows < current_index + 1)
        {
            X.Dimension(current_index + 1, nPop);
        }

//		X[current_index][0] = 1;
        bool status = add1kgMafToX(X, markername, current_index);
        // now check distatnce
        double dist = 0;
        double maf = current_AC / current_N / 2;

        if (!status)
        { // not found in 1kg marker, use 1/4N or skip
            if (matchOnly)
            {
                if (matchonly_skipped == 0)
                {
                    printf("With --matchOnly flag, the following markers cannot be found in 1000G VCF and will be excluded in population stratification correction:\n");
                }
                printf("    %s \n", markername.c_str());
                matchonly_skipped++;
                continue;
            }
            for (int i = 0; i < nPop; i++)
            {
                X[current_index][i] = (double) 0.25 / current_N;
            }
        } else
        { // check matchDist and MAF
            // check MAF first
            if (maf < minMatchMAF || maf > maxMatchMAF)
            {
                if (match_maf_skipped == 0)
                {
                    printf("With --minMatchMAF=%g,--maxMatchMAF=%g, markers with MAF out of this range will be excluded in population correction!\n",
                           minMatchMAF, maxMatchMAF);
                }
                match_maf_skipped++;
                continue;
            }
            // then check dist
            if (matchByAbs)
            { // dist = min( |MAF_study - MAF_pop| )
                for (int i = 0; i < nPop; i++)
                {
                    double di = abs(maf - X[current_index][i]);
                    if (dist == 0)
                    {
                        dist = di;
                    } else
                    {
                        if (di < dist)
                        {
                            dist = di;
                        }
                    }
                }
            } else
            { // dist = sqrt( sum(MAF_study - MAF_pop)^2 )
                for (int i = 0; i < nPop; i++)
                {
                    dist += (maf - X[current_index][i]) * (maf - X[current_index][i]);
                }
                dist = sqrt(dist / nPop);
            }
            if (matchDist > 0 && dist > matchDist)
            {
                if (dist_skipped == 0)
                {
                    printf("With --matchDist %g flag, the following markers with dist>%g in 1000G VCF will be excluded in population stratification correction:\n",
                           matchDist, matchDist);
                }
                if (debug)
                {
                    printf("    %s, dist=%g\n", markername.c_str(), dist);
                }
                if (dist_skipped == 0)
                {
                    if (matchByAbs)
                    {
                        printf("Distance calculation method: dist = min( |MAF_study - MAF_pop| )\n\n");
                    } else
                    {
                        printf("Distance calculation method: dist = sqrt( sum(MAF_study - MAF_pop)^2 )\n\n");
                    }
                }
                dist_skipped++;
                continue;
            }
        }

        Y.Push(current_AC / current_N / 2);
        if (debug)
        {
            ifprintf(match_id_file, "%s %g", markername.c_str(), current_AC / current_N / 2);
            for (int i = 0; i < nPop; i++)
            {
                ifprintf(match_id_file, " %g", X[current_index][i]);
            }
            ifprintf(match_id_file, "\n");
        }
        current_index++;
//	printf("index=%d\n",current_index);
//		if (current_index >= Nmarkers)
//			break;
    }
    ifclose(file);
    printf("Finish loading and match to 1000G VCF!\n");
    if (matchOnly)
    {
        printf("In pop correction, total %d variants that do not have a match in 1000G were skipped.\n",
               matchonly_skipped);
    }
    if (matchDist > 0)
    {
        printf("In pop correction, total %d variants with dist>%g were skipped.\n", dist_skipped, matchDist);
    }
    if (match_maf_skipped > 0)
    {
        printf("In pop correction, total %d variants were skipped because maf is not in range of %g~%g",
               match_maf_skipped, minMatchMAF, maxMatchMAF);
    }
    if (debug)
    {
        ifclose(match_id_file);
    }
    if (matchOnly && X.rows == Y.Length() + 1)
    { // if last element is not added
        X.Dimension(Y.Length(), nPop);
    }
//	if (current_index < Nmarkers)
//		X.Dimension( current_index,nPop );
    hashToRemoveDuplicates.Clear();

    if (X.rows > 10000)
    {
        printf("Warning: regressing %d markers for population stratification. This process might be very slow...\n\n",
               X.rows);
    }
    // perform regression & save to fk
    setGammaFromRegression(study, X, Y);
}

bool Meta::add1kgMafToX(Matrix &X, String &markername, int current_index)
{
//printf("markername=%s\n",markername.c_str());
    // match allele first
    std::map<String, std::vector<double> >::iterator ptr;
    ptr = af1KG.find(markername);
    if (ptr == af1KG.end())
    { // check if flip
        StringArray tokens;
        tokens.AddTokens(markername, ":");
        String flip_markername = tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" + tokens[2];
        ptr = af1KG.find(flip_markername);
        if (ptr == af1KG.end())
        { // does not exist...
            return false;
        } else
        { // make a flip vector
            for (int i = 0; i < nPop; i++)
            {
                X[current_index][i] = 1 - ptr->second[i];
            }
        }
    } else
    {
        for (int i = 0; i < nPop; i++)
        {
            X[current_index][i] = ptr->second[i];
        }
    }
    return true;
}

// do regression, set coefficient
// fit linear model with beta >=0
void Meta::setGammaFromRegression(int study, Matrix &X, Vector &Y)
{
    if (X.rows != Y.Length())
    {
        error("[Meta::setGammaFromRegression] X.row=%d, Y length=%d. Something is wrong...\n", X.rows, Y.Length());
    }
    if (Y.Length() == 0)
    {
        error("No available non-zero maf in study #%d. Please check your input file!\n", study + 1);
    }

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
    inv.Product(transX, X);
    SVD svd;
    svd.InvertInPlace(inv);
    // X'Y
    Vector XY;
    XY.Dimension(nPop);
    for (int i = 0; i < nPop; i++)
    {
        XY[i] = Y.InnerProduct(transX[i]);
    }
//	printf("study=%d",study);
    for (int i = 0; i < nPop; i++)
    {
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

/**
 * Read through summary statistics of each study and pool the information.
 * At the end, single variant meta-analysis will be completed.
 * Information will be stored for gene-based analysis
 * @param group
 */
void Meta::PoolSummaryStat(GroupFromAnnotation &group)
{
    total_N = 0;
    flipCount = 0;
    skip_count = 0;

    if (relateBinary)
    {
        Const_binary.Dimension(scorefile.Length());
        for (int s = 0; s < scorefile.Length(); s++)
        {
            double sigma = Ysigma2g[s] + Ysigma2e[s];
//			Const_binary[s] = getCorrectC( Ymean[s],sigma );
        }
    }

    if (useExactMetaMethod)
    {
        Ydelta.Dimension(scorefile.Length());
        Ysigma2.Dimension(scorefile.Length());
        for (int s = 0; s < scorefile.Length(); s++)
        { // no need to initialize ydelta
            Ysigma2[s] = -1;
        }
        for (int s = 0; s < scorefile.Length(); s++)
        {
            bool y_status = updateYstat(s);
            if (!y_status)
            {
                error("cannot update y stat at study #%d\n", s + 1);
            }
        }
        // calculate overall mean & adjust
        int n = 0;
        for (int s = 0; s < scorefile.Length(); s++)
        {
            n += SampleSize[s];
        }
        double ymean = 0;
        for (int s = 0; s < scorefile.Length(); s++)
        {
            ymean += (double) SampleSize[s] / (double) n * Ydelta[s];
        }
        for (int s = 0; s < scorefile.Length(); s++)
        {
            Ydelta[s] -= ymean;
        }
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
    for (int study = 0; study < scorefile.Length(); study++)
    {
        int duplicateSNP = 0;
        //Set up data structure to help remove duplicates
        StringIntHash hashToRemoveDuplicates; //hash for skip duplicate markers
        varCount = 0;
        dupcheck_chr = "";

        Vector chisq_study_i;

        printf("Pooling summary statistics from study %d ...\n", study + 1);

        //read in summary statistics.
        //maf and summary stat are saved. SNPpos in the file are hashed.
        SummaryFileReader covReader;
        covReader.marker_col = this->marker_col;
        covReader.cov_col = this->cov_col;
        if (cond != "" && cond_status[study])
        {
            String covFilename = covfile[study];
            bool cov_status = covReader.ReadTabix(covFilename);
            if (!cov_status)
            {
                error("Cannot open cov file: %s\n", covFilename.c_str());
            }
        }

        String filename = scorefile[study];
        IFILE file;
        file = ifopen(filename, "r");
        if (file == NULL)
        {
            error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.summary.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.summary.file.gz\n",
                  scorefile[study].c_str());
        }

        String buffer;
        buffer.ReadLine(file);
        bool adjust;
        tellRvOrRmw(buffer, adjust, marker_col, cov_col);
//		if (!status)
//			error("File %s is neither RMW or rvtest score file!\n", filename.c_str());
        ifclose(file);
        FormatAdjust.Push(adjust);

        file = ifopen(filename, "r");
        bool pass_header = 0;
        while (!ifeof(file))
        {
            buffer.ReadLine(file);
            if (buffer.Length() < 2)
            {
                continue;
            }
            if (!pass_header)
            {
                if (SampleSize[study] == -1)
                {
                    if (buffer.Find("##AnalyzedSamples") != -1)
                    {
                        StringArray tokens;
                        tokens.AddTokens(buffer, "=");
                        SampleSize[study] = tokens[1].AsInteger();
                        continue;
                    }
                }
                if (buffer.FindChar('#') == -1 && buffer.Find("CHROM") == -1)
                {
                    pass_header = 1;
                }
            }
            if (!pass_header)
            {
                continue;
            }
            double current_chisq;
            bool status = poolSingleRecord(study, current_chisq, duplicateSNP, adjust, buffer, covReader);
            if (status)
            {
                chisq_study_i.Push(current_chisq);
            }
        }
        ifclose(file);

        total_N += SampleSize[study];
        //calculate GC
        chisq_study_i.Sort();
        if (chisq_study_i.size == 0)
        {
            printf("\nWarning: no GC calculated in chisq_study_i!\n");
            GCbyStudy.Push(0);
        } else
        {
            GCbyStudy.Push(chisq_study_i[0.5] / 0.456);
        }
        if (duplicateSNP > 0)
        {
            printf("\nWarning: %d variants are skipped from analysis due to duplicate records in study %s. Please check log for details.\n\n",
                   duplicateSNP, scorefile[study].c_str());
            fprintf(log,
                    "\nWarning: %d variants are skipped from analysis due to duplicate records in study %s. Please check log for details.\n\n",
                    duplicateSNP, scorefile[study].c_str());
        }
        printf("  done\n");
    }


    if (skip_count > 1)
    {
        printf("Warning: %d variants have been excluded from analysis due to un-matched alleles. Please check log file for details.\n",
               skip_count);
    }

    //finalize direction array
    for (int i = 0; i < directions.Length(); i++)
    {
        while (directions[i].Length() < scorefile.Length())
        {
            directions[i] += '?';
        }
    }

    //calculate pooled allele frequencies
    setPooledAF();

    if (bHeterogeneity) {
      printf("\nRe-processing studies for heterogeneity analysis ...\n");
      for (int study = 0; study < scorefile.Length(); study++) {
        SummaryFileReader covReader;
        covReader.marker_col = this->marker_col;
        covReader.cov_col = this->cov_col;

        // We only need the covariance matrix here if conditional analysis was requested.
        if (cond != "" && cond_status[study]) {
          String covFilename = covfile[study];
          bool cov_status = covReader.ReadTabix(covFilename);
          if (!cov_status) {
            error("Cannot open cov file: %s\n", covFilename.c_str());
          }
        }

        // Open file of score statistics again.
        String filename = scorefile[study];
        IFILE file;
        file = ifopen(filename, "r");

        // Get rid of first header line. We already know what type of file this is from earlier.
        String buffer;
        buffer.ReadLine(file);
        bool adjust = FormatAdjust[study];

        // Loop over lines of score file, updating heterogeneity statistic for each variant.
        bool pass_header = false;
        while (!ifeof(file)) {
          buffer.ReadLine(file);
          if (buffer.Length() < 2) {
            continue;
          }
          if (!pass_header) {
            if (buffer.FindChar('#') == -1 && buffer.Find("CHROM") == -1) {
              pass_header = true;
            }
          }
          if (!pass_header) {
            continue;
          }

          poolHeterogeneity(study, adjust, buffer, covReader);
        }
        ifclose(file);
      }
    }
}

void Meta::WriteSingleVariantResults(GroupFromAnnotation &group) {
  printf("\nPerforming Single variant meta analysis ...\n");

  //calculate final results here
  IFILE output;
  String filename;
  IFILE vcfout;
  String vcf_filename;

  if (prefix == "-") {
    filename = "-";
  }
  else if (prefix == "") {
    filename = "meta.singlevar.results";
  }
  else if (prefix.Last() == '.' || prefix.Last() == '/') {
    filename = prefix + "meta.singlevar.results";
  }
  else {
    filename = prefix + ".meta.singlevar.results";
  }

  output = ifopen(filename, "w", InputFile::UNCOMPRESSED);

  printSingleMetaHeader(filename, output);
  if (outvcf) {
    printOutVcfHeader(vcf_filename, vcfout);
  }

  //for annotation purpose
  target_chr = "";
  target_pos = 0, target = 0;
  target_pvalue = _NAN_;
  //Sort variants by chr and pos
  for (int i = 0; i < SNPmaf_maf.Length(); i++) {
    printSingleMetaVariant(group, i, output, vcfout);
  }

  plotSingleMetaGC(output, 1);
  if (cond != "") {
    plotSingleMetaGC(output, 0);
  }
  printf("\n  done.\n\n");

  ifclose(output);
  if (outvcf) {
    ifclose(vcfout);
    printf("\n  VCF file based on superset of variants from pooled studies has been saved \n    %s\n", vcf_filename.c_str());
  }
}

void Meta::Run(GroupFromAnnotation &group)
{
    if (outvcf)
    {
        return;
    }

    if (!SKAT && !Burden && !MB && !MAB && !BBeta && !VTa)
    {
        printf("\nWarning: none of the gene-level tests was requested; only single variant meta-analysis was done.\n");
        return;
    }

    maf = new Vector[group.annoGroups.Length()];
    stats = new Vector[group.annoGroups.Length()];

    if (cond != "")
    {
        cond_stats = new Vector[group.annoGroups.Length()];
        cond_cov = new Matrix[group.annoGroups.Length()];
    }
    cov = new Matrix[group.annoGroups.Length()];
    singlePvalue = new Vector[group.annoGroups.Length()];
    singleEffSize = new Vector[group.annoGroups.Length()];

    if (group.annoGroups.Length() <= 0)
    {
        error("No valid groups. Unable to do gene-based test. Check if your group file is empty!\n");
    }

    loadSingleStatsInGroup(group);

    printf("\nChecking if all groups are analyzed...\n");
    fprintf(log, "\nChecking if all groups are analyzed...\n");
    int geneCounter = 0;
    for (int g = 0; g < group.annoGroups.Length(); g++)
    {
        if (maf[g].Length() == 0)
        {
            geneCounter++;
            if (geneCounter == 1)
            {
                printf("\t\tPlease check log file for groups that are not analyzed.\n");
                fprintf(log, "Warning: the following groups has no qualifed variants to group and are skipped:\n\t");
            }
            fprintf(log, "%s, ", group.annoGroups[g].c_str());
        }
    }
    printf("\n\tdone!\n");
    fprintf(log, "\n\tdone!\n");

    //loop through cov matrices of all studies and update cov
    loadSingleCovInGroup(group);

    // update stats in group
    for (int g = 0; g < group.annoGroups.Length(); g++)
    {
        for (int i = 0; i < cov[g].rows; i++)
        {
            for (int j = i + 1; j < cov[g].cols; j++)
            {
                cov[g][j][i] = cov[g][i][j];
            }
        }
        if (cond != "")
        {
            for (int i = 0; i < cond_cov[g].rows; i++)
            {
                for (int j = i + 1; j < cond_cov[g].cols; j++)
                {
                    cond_cov[g][j][i] = cond_cov[g][i][j];
                }
            }
        }
    }
    String method = "";
    if (Burden)
    {
        method = "burden";
        BurdenAssoc(method, group, maf, stats, cond_stats, cov, cond_cov, singleEffSize, singlePvalue);
    }
    if (MB)
    {
        method = "MB";
        BurdenAssoc(method, group, maf, stats, cond_stats, cov, cond_cov, singleEffSize, singlePvalue);
    }
    if (MAB)
    {
        method = "MAB";
        BurdenAssoc(method, group, maf, stats, cond_stats, cov, cond_cov, singleEffSize, singlePvalue);
    }
    if (BBeta)
    {
        method = "BBeta";
        BurdenAssoc(method, group, maf, stats, cond_stats, cov, cond_cov, singleEffSize, singlePvalue);
    }

    if (VTa)
    {
        VTassoc(group);
    }

    if (SKAT)
    {
        SKATassoc(group);
    }

    pdf.CloseFile();

    //housekeeping
    if (maf)
    { delete[] maf; }
    if (stats)
    { delete[] stats; }
    if (cov)
    { delete[] cov; }
    if (singlePvalue)
    { delete[] singlePvalue; }
    if (singleEffSize)
    { delete[] singleEffSize; }
    printf("\nQQ plots and manhattan polts have been saved in %s.\n", pdf_filename.c_str());
}

//calculate the genotype covariance between G and X
void Meta::CalculateGenotypeCov(SummaryFileReader &covReader, String chr, int pos, int study, Vector &result)
{
    result.Dimension(commonVar_study[study].Length(), 0.0);
    for (int cv = 0; cv < commonVar_study[study].Length(); cv++)
    {
        if (chr != common_chr[commonVar_study[study][cv]])
        {
            continue;
        }
        if (abs(common_pos[commonVar_study[study][cv]] - pos) > 1000000)
        {
            continue;
        }
        if (common_pos[commonVar_study[study][cv]] > pos)
        {
            if (!covReader.ReadRecord(chr, pos))
            {
                continue;
            }
            StringArray tmp;
            tmp.AddTokens(covReader.marker_nearby, ",");
            String pos_str;
            pos_str = common_pos[commonVar_study[study][cv]];
            int loc = tmp.Find(pos_str);
            if (loc == -1)
            {
                continue;
            }
            tmp.Clear();
            tmp.AddTokens(covReader.marker_cov, ",");
            result[cv] = tmp[loc].AsDouble() * SampleSize[study];
        } else
        {
            int loc = commonVar_markers_in_window[study][cv].Find(pos);
            if (loc == -1)
            {
                continue;
            }
            result[cv] = commonVar_marker_cov[study][cv][loc] * SampleSize[study];
            /*
			if(pos==79760530)
			printf("chr %s pos %d GXcov %g\n",chr.c_str(),pos,commonVar_marker_cov[study][cv][loc]);
		*/    }
    }
}

double
Meta::GrabGenotypeCov(SummaryFileReader &covReader, int study, String chr1, String pos1, String chr2, String pos2,
                      String &SNP1, String &SNP2)
{
    double result = 0.0;
    if (chr1 != chr2)
    {
        return result;
    }

    int skip = SNPexclude.Integer(String(std::to_string(study).c_str()) + ":" + chr1 + ":" + pos1);
    if (skip != -1)
    {
        return result;
    }
    skip = SNPexclude.Integer(String(std::to_string(study).c_str()) + ":" + chr2 + ":" + pos2);
    if (skip != -1)
    {
        return result;
    }

    StringArray flip, tmp;
    int marker_idx;

    if (pos1.AsInteger() < pos2.AsInteger())
    {
        //Check in current hash table first
        marker_idx = covReader.markerPosHash.Integer(chr1 + ":" + pos1);
        //if this record has been hashed
        if (marker_idx != -1)
        {
            tmp.AddTokens(covReader.markerNearby[marker_idx], ",");
            int loc = tmp.Find(pos2);
            if (loc == -1)
            {
                return result;
            }
            tmp.Clear();
            tmp.AddTokens(covReader.markerNearbyCov[marker_idx], ",");
            result = tmp[loc].AsDouble() * SampleSize[study];
            return result;
        }
        if (!covReader.ReadRecord(chr1, pos1.AsInteger()))
        {
            return result;
        }
        tmp.Clear();
        tmp.AddTokens(covReader.marker_nearby, ",");
        int loc = tmp.Find(pos2);
        if (loc == -1)
        {
            return result;
        }
        tmp.Clear();
        tmp.AddTokens(covReader.marker_cov, ",");
        result = tmp[loc].AsDouble() * SampleSize[study];
    } else
    {
        //Check in current hash table first
        marker_idx = covReader.markerPosHash.Integer(chr2 + ":" + pos2);
        //if this record has been hashed
        if (marker_idx != -1)
        {
            StringArray tmp;
            tmp.AddTokens(covReader.markerNearby[marker_idx], ",");
            int loc = tmp.Find(pos1);
            if (loc == -1)
            {
                return result;
            }
            tmp.Clear();
            tmp.AddTokens(covReader.markerNearbyCov[marker_idx], ",");
            result = tmp[loc].AsDouble() * SampleSize[study];
            return result;
        }
        if (!covReader.ReadRecord(chr2, pos2.AsInteger()))
        {
            return result;
        }
        tmp.Clear();
        tmp.AddTokens(covReader.marker_nearby, ",");
        int loc = tmp.Find(pos1);
        if (loc == -1)
        {
            return result;
        }
        tmp.Clear();
        tmp.AddTokens(covReader.marker_cov, ",");
        result = tmp[loc].AsDouble() * SampleSize[study];
    }

    //if this marker is flipped then markers from the entire row
    //should have cov multiply by -1.0.
    double factor1 = 1.0, factor2 = 1.0;
    if (flipSNP.Integer(String(std::to_string(study).c_str()) + ":" + chr1 + ":" + pos1) != -1)
    {
        factor1 = -1.0;
    }
    if (flipSNP.Integer(String(std::to_string(study).c_str()) + ":" + chr2 + ":" + pos2) != -1)
    {
        factor2 = -1.0;
    }
    result *= (factor1 * factor2);

    return result;
}

void Meta::CalculateXXCov(int study, Matrix &result)
{
    for (int i = 0; i < result.rows; i++)
    {
        result[i][i] = commonVar_marker_cov[study][i][0] * SampleSize[study];
        for (int j = i + 1; j < result.cols; j++)
        {
            if (common_chr[commonVar_study[study][j]] != common_chr[commonVar_study[study][i]])
            {
                continue;
            }
            if (common_pos[commonVar_study[study][i]] < common_pos[commonVar_study[study][j]])
            {
                int loc = commonVar_markers_in_window[study][i].Find(common_pos[commonVar_study[study][j]]);
                if (loc == -1)
                {
                    continue;
                }
                result[i][j] = result[j][i] = commonVar_marker_cov[study][i][loc] * SampleSize[study];
            } else
            {
                int loc = commonVar_markers_in_window[study][j].Find(common_pos[commonVar_study[study][i]]);
                if (loc == -1)
                {
                    continue;
                }
                result[i][j] = result[j][i] = commonVar_marker_cov[study][j][loc] * SampleSize[study];
            }
        }
    }
    SVD svd;
    svd.InvertInPlace(result);
}

void Meta::parseScoreFiles() {
  if (scorefile.Length() > 0) {
    return; // We already added scorefiles directly, likely in unit tests
  }

  // summary file
  if (summaryFiles != "")
  {
    IFILE inFile = ifopen(summaryFiles, "r");
    if (inFile == NULL)
    {
      error("FATAL ERROR! Please check file name for --summaryFiles  %s\n", summaryFiles.c_str());
    }
    String buffer;
    while (!ifeof(inFile))
    {
      buffer.ReadLine(inFile);
      if (buffer.FindChar('#') != -1)
      {
        continue;
      }
      scorefile.AddTokens(buffer, "\n");
    }
    ifclose(inFile);
  } else
  {
    error("FATAL ERROR! --summaryFiles can not be empty! Usage: --summaryFiles your.list.of.summary.files\n");
  }
}

void Meta::parseCovFiles() {
  if (covfile.Length() > 0) {
    return; // We already added covariance files directly, likely in unit tests
  }

  // cov file
  if (covFiles != "")
  {
    IFILE inFile = ifopen(covFiles, "r");
    if (inFile == NULL)
    {
      error("Cannot open file %s\n", covFiles.c_str());
    }
    String buffer;
    while (!ifeof(inFile))
    {
      buffer.ReadLine(inFile);
      if (buffer.FindChar('#') != -1)
      {
        continue;
      }
      covfile.AddTokens(buffer, "\n");
    }
    ifclose(inFile);
    //check if summary files and cov files have the same length
    if (scorefile.Length() != covfile.Length())
    {
      error("There are %d summary files and %d covariance files. Please check to make sure the same number of files are included in the list.\n");
    }
  }
  else if (Burden || MB || VTa || SKAT)
  {
    error("Covariance files are essential to do gene-level tests. Please use --covFiles your.list.of.cov.files option.\n");
  }
}

/**
 * prepare for conditional analysis and ensure that the specified variants are actually present in this set of studies
 */
void Meta::prepareConditionalAnalysis()
{
    // read list of cond variants
    IFILE condFile = ifopen(cond, "r");
    if (condFile == NULL)
    {
        error("Can not open file %s.\n", cond.c_str());
    }

    String buffer;
    while (!ifeof(condFile))
    {
        buffer.ReadLine(condFile);
        if (buffer.FindChar('#') != -1 || buffer.FindChar(':') == -1)
        {
            continue;
        }
        StringArray tmpMarker;
        tmpMarker.AddTokens(buffer, " \t");
        for (int i = 0; i < tmpMarker.Length(); i++)
        {
            commonVar.Push(tmpMarker[i]);
            StringArray tmp;
            tmp.AddTokens(tmpMarker[i], ":");
            common_chr.Push(tmp[0]);
            common_pos.Push(tmp[1]);
            common_ref.Push(tmp[2]);
            common_alt.Push(tmp[3]);
            conditionVar.SetInteger(tmp[0] + ":" + tmp[1], commonVar.Length() - 1);
        }
    }
    ifclose(condFile);
    if (common_pos.Length() == 0)
    {
        error("No variant to be conditioned on. Condition file %s might be in the wrong format.\n", cond.c_str());
    }

    // set variable
    int numStudy = scorefile.Length();
    commonVar_study = new IntArray[numStudy];
    commonVar_effSize = new Vector[numStudy];
    commonVar_U = new Vector[numStudy];
    commonVar_V = new Vector[numStudy];
    commonVar_markers_in_window = new IntArray *[numStudy];
    commonVar_marker_cov = new Vector *[numStudy];
    XX_inv = new Matrix[numStudy];
    cond_status = new bool[numStudy];
    commonVar_betaHat = new Vector[numStudy];

    // set marker stats
    bool status = setCondMarkers();
    if (!status)
    {
        cond = "";
        printf("\nWarning: None of the variants to be conditioned on are found in any of the studies. Conditional analysis options will be ignored.\n\n");
    }
}


// read cond marker info from score & cov file
// if no marker info was load, return false
bool Meta::setCondMarkers()
{
    for (int s = 0; s < scorefile.Length(); s++)
    {
        SummaryFileReader statReader, covReader;
        statReader.marker_col = this->marker_col;
        statReader.cov_col = this->cov_col;
        covReader.marker_col = this->marker_col;
        covReader.cov_col = this->cov_col;

        String filename = scorefile[s];
        String covFilename = covfile[s];

        if (!statReader.ReadTabix(filename))
        {
            error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.summary.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.summary.file.gz\n",
                  scorefile[s].c_str());
        }

        bool adjust;
        String buffer;
        IFILE dup = ifopen(filename, "r");
        buffer.ReadLine(dup);
        tellRvOrRmw(buffer, adjust, marker_col, cov_col);
        ifclose(dup);

        if (SampleSize[s] == -1)
        {
            dup = ifopen(filename, "r");
            while (!ifeof(dup))
            {
                buffer.ReadLine(dup);
                //get sample size
                if (buffer.Find("##AnalyzedSamples") != -1)
                {
                    StringArray tokens;
                    tokens.AddTokens(buffer, "=");
                    SampleSize[s] = tokens[1].AsInteger();
                    break;
                }
            }
            ifclose(dup);
        }

        if (!covReader.ReadTabix(covFilename))
        {
            error("Cannot open file: %s!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip your.covariance.file; tabix -c \"#\" -s 1 -b 2 -e 2 your.covariance.file.gz.\n",
                  covfile[s].c_str());
        }

        for (int i = 0; i < commonVar.Length(); i++)
        {
            //if this variant is genotyped in this study
            bool cov_status = covReader.ReadRecord(common_chr[i], common_pos[i]);
            if (!cov_status)
            {
                continue;
            }
            if (!statReader.ReadRecord(common_chr[i], common_pos[i]))
            {
                continue;
            }
            StringArray record;
            record.AddTokens(statReader.buffer, "\t");

            if ((record[2] == common_ref[i] && record[3] == common_alt[i]) ||
                (record[3] == common_ref[i] && record[2] == common_alt[i]))
            {
                double v = record[14 - adjust].AsDouble();
                if (v > 0)
                {
                    commonVar_study[s].Push(i);
                    commonVar_V[s].Push(v * v);
                    commonVar_effSize[s].Push(record[15 - adjust].AsDouble());
                    commonVar_U[s].Push(record[13 - adjust].AsDouble());
                }
            }
        }

        int dim = commonVar_study[s].Length();
        if (dim != 0)
        {
            cond_status[s] = true;
        } else
        {
            cond_status[s] = false;
            continue;
        }
        commonVar_markers_in_window[s] = new IntArray[dim];
        commonVar_marker_cov[s] = new Vector[dim];

        bool cov_status = covReader.ReadTabix(covFilename);
        if (!cov_status)
        {
            error("Cannot open cov file: %s\n", covFilename.c_str());
        }

        for (int i = 0; i < dim; i++)
        {
            int idx = commonVar_study[s][i];
            if (!covReader.ReadRecord(common_chr[idx], common_pos[idx]))
            {
                continue;
            }
            StringArray tmp_markers;
            tmp_markers.AddTokens(covReader.marker_nearby, ",");
            for (int j = 0; j < tmp_markers.Length(); j++)
            {
                commonVar_markers_in_window[s][i].Push(tmp_markers[j].AsInteger());
            }
            tmp_markers.Clear();
            tmp_markers.AddTokens(covReader.marker_cov, ",");
            for (int j = 0; j < tmp_markers.Length(); j++)
            {
                commonVar_marker_cov[s][i].Push(tmp_markers[j].AsDouble());
            }
        }
        XX_inv[s].Dimension(commonVar_study[s].Length(), commonVar_study[s].Length(), 0.0);
        CalculateXXCov(s, XX_inv[s]);
        //calculate beta hat for common variants
        commonVar_betaHat[s].Dimension(dim);
        for (int i = 0; i < commonVar_betaHat[s].dim; i++)
        {
            commonVar_betaHat[s][i] = XX_inv[s][i].InnerProduct(commonVar_U[s]);
        }
    }

    // check status of cond markers
    bool status = false;
    for (int i = 0; i < scorefile.Length(); i++)
    {
        if (cond_status[i])
        {
            status = true;
            break;
        }
    }
    return status;
}

// read RMW output header that records y
// update y and y-mean
bool Meta::updateYstat(int study)
{
    bool found_y = 0;
    int status = 0;
    String fname = scorefile[study];
    IFILE file = ifopen(fname, "r");
    if (file == NULL)
    {
        error("Cannot open file: %s!\n", fname.c_str());
    }
    while (!ifeof(file))
    {
        String buffer;
        buffer.ReadLine(file);
        if (SampleSize[study] == -1)
        {
            if (buffer.Find("##AnalyzedSamples") != -1)
            {
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
        if (buffer.Find("#") == -1)
        {
            break;
        }
        if (found_y)
        { // name min 25% median 75th max mean(7th) variance(8th)
            StringArray tokens;
            tokens.AddTokens(buffer, "\t ");
            Ydelta[study] = tokens[6].AsDouble();
            Ysigma2[study] = tokens[7].AsDouble();
            status++;
            found_y = 0;
            if (status == 1)
            {
                break;
            }
        }
        if (buffer.Find("##TraitSummaries") != -1)
        {
            found_y = 1;
            continue;
        }
    }
    if (status == 1)
    {
        return 1;
    } else
    {
        return 0;
    }
}

void Meta::UpdateDirection(int &direction_idx, int study, char marker_direction, String &chr_pos, bool exclude)
{
    //if this chr_pos is not shown before, but study>0, then assume studies before this study are all monomorphic at this site.
    if (direction_idx == -1)
    {
        if (study > 0)
        {
            String tmp = "";
            for (int i = 0; i < study; i++)
            {
                tmp += '?';
            }
            tmp += marker_direction;
            directions.Push(tmp);
            directionByChrPos.SetInteger(chr_pos, directions.Length() - 1);
        } else if (study == 0)
        {
            String tmp;
            tmp = marker_direction;
            directions.Push(tmp);
            directionByChrPos.SetInteger(chr_pos, directions.Length() - 1);
        }
        if (exclude)
        {
            refalt.Push(".:.");
        }
    } else
    {
        int l = directions[direction_idx].Length();
        if (l == study)
        {
            directions[direction_idx] += marker_direction;
        }
        if (l < study)
        {
            for (int i = 0; i < study - l; i++)
            {
                directions[direction_idx] += '?';
            }
            directions[direction_idx] += marker_direction;
        }
        if (l > study)
        {
            printf("WARNING:there is duplicate marker in some of the studies at %s.\n", chr_pos.c_str());
            fprintf(log, "WARNING:there is duplicate marker in some of the studies at %s.\n", chr_pos.c_str());
            directions[direction_idx] += marker_direction;
        }
    }
}


void Meta::UpdateExcludedMarker(int &study, String &chr_pos, int filter, String markername)
{
    skip_count++;
    String skip_SNP;
    skip_SNP = study;
    skip_SNP += ":";
    skip_SNP += chr_pos;

    SNPexclude.SetInteger(skip_SNP, skip_count);
    if (filter == 0)
    {
        fprintf(log,
                "Warning: variant %s from study %d failed to pass HWE or CALLRATE filter and it is excluded from meta analysis.\n",
                markername.c_str(), study);
        return;
    }
    if (filter == 1)
    {
        fprintf(log,
                "Warning: variant %s from study %d has at least one allele missing but is polymorphic; and is excluded from meta-analysis.\n",
                markername.c_str(), study);
        return;
    }
    if (filter == 2)
    {
        int idx = directionByChrPos.Integer(chr_pos);
        fprintf(log, "Warning: variant %s from study %d is excluded due to unmatched alleles (expecting to see %s).\n",
                markername.c_str(), study, refalt[idx].c_str());
        return;
    }
}

void Meta::UpdateStrIntHash(String &chr_pos, int val, StringIntHash &sihash)
{
    int idx = sihash.Find(chr_pos);
    if (idx == -1)
    {
        sihash.SetInteger(chr_pos, val);
    } else
    {
        int old_N = sihash.GetCount(idx);
        old_N += val;
        sihash.SetInteger(idx, old_N);
    }
}

void Meta::UpdateStrDoubleHash(String &chr_pos, double val, StringDoubleHash &sdhash)
{
    int idx = sdhash.Find(chr_pos);
    if (idx == -1)
    {
        sdhash.SetDouble(chr_pos, val);
    } else
    {
        double old_N = sdhash.Double(chr_pos);
        old_N += val;
        sdhash.SetDouble(chr_pos, old_N);
    }
}

void Meta::UpdateACInfo(String &chr_pos, double AC)
{
    int idx = usefulAC.Find(chr_pos);
    if (idx == -1)
    {
        usefulAC.SetDouble(chr_pos, AC);
    } else
    {
        int old_N = usefulAC.Double(idx);
        old_N += AC;
        usefulAC.SetDouble(idx, old_N);
    }
}


// simply add up u and v
// for exact method, v passed here is vk/sigmak. V will be further adjusted in print meta record
void Meta::UpdateStats(int study, String &markerName, double stat, double vstat, bool flip)
{
    //update SNPstat
    double flip_factor = 1.0;
    if (flip)
    {
        flip_factor = -1.0;
    }
    vstat = vstat * vstat;
    stat *= flip_factor;
    if (useExactMetaMethod)
    { // differnt from u & v RMW uses
        stat *= Ysigma2[study];
        vstat *= Ysigma2[study]; // note we're recording Vexact/sigma2 here
    }
    int stat_idx = SNPstat.Find(markerName);
    if (stat_idx < 0)
    {
        SNPstat.SetDouble(markerName, stat);
    } else
    {
        double prev = SNPstat.Double(stat_idx);
        prev += stat;
        SNPstat.SetDouble(markerName, prev);
    }

    //update SNP_Vstat
    stat_idx = SNP_Vstat.Find(markerName);
    if (stat_idx == -1)
    {
        SNP_Vstat.SetDouble(markerName, vstat);
    } else
    {
        double prev = SNP_Vstat.Double(stat_idx);
        prev += vstat;
        SNP_Vstat.SetDouble(markerName, prev);
    }
}

/**
 * Update heterogeneity statistic for a particular variant given the score statistic and variance from a particular study.
 * @param study
 * @param markerName Variant in chr:pos format :(
 * @param stat Variant score statistic
 * @param sqrt_v Variance of score statistic
 * @param flip Whether the score statistic should be flipped because effect allele did not match the first analyzed study's effect allele
 */
void Meta::UpdateHetStats(int study, String &markerName, double stat, double sqrt_v, bool flip)
{
  double flip_factor = 1.0;
  if (flip) {
    flip_factor = -1.0;
  }
  stat *= flip_factor;

  // Get the meta-analysis U and V
  double u_meta = SNPstat.Double(markerName);
  double v_meta = SNP_Vstat.Double(markerName);

  // Calculate het statistic
  double v = sqrt_v * sqrt_v;
  double z = stat / v;
  double e = u_meta / v_meta;
  double het = (z - e) * (z - e) * v;

  // Add het stat onto running total
  int stat_idx = SNP_heterog_stat.Find(markerName);
  if (stat_idx < 0) {
    SNP_heterog_stat.SetDouble(markerName, het);
  } else {
    double prev = SNP_heterog_stat.Double(stat_idx);
    prev += het;
    SNP_heterog_stat.SetDouble(markerName, prev);
  }

  // Increase degrees of freedom of het stat by 1 (each study adds 1 to the total df)
  int df_idx = SNP_heterog_df.Find(markerName);
  if (df_idx < 0) {
    SNP_heterog_df.SetInteger(markerName, 1);
  } else {
    int prev = SNP_heterog_df.Integer(df_idx);
    prev += 1;
    SNP_heterog_df.SetInteger(markerName, prev);
  }
}

void Meta::UpdateHetCondStats(int study, bool flip, int adjust, String &markerName, StringArray &tokens, SummaryFileReader &covReader) {
  Vector GX;
  CalculateGenotypeCov(covReader, tokens[0], tokens[1].AsInteger(), study, GX);

  // Get the meta-analysis conditional U and V
  double u_meta = SNPstat_cond.Double(markerName);
  double v_meta = SNP_Vstat_cond.Double(markerName);

  // Calculate this study's conditional U and V
  double cond_u = tokens[13 - adjust].AsDouble() - GX.InnerProduct(commonVar_betaHat[study]);
  Vector tmp;
  for (int i = 0; i < GX.dim; i++) {
    tmp.Push(GX.InnerProduct(XX_inv[study][i]));
  }
  double v = tokens[14 - adjust].AsDouble();
  double cond_v_part = tmp.InnerProduct(GX);
  double cond_v = v * v - cond_v_part;

  // Calculate het statistic
  double z = cond_u / cond_v;
  double e = u_meta / v_meta;
  double het = (z - e) * (z - e) * cond_v;

  // Add het stat onto running total
  int stat_idx = SNP_heterog_cond_stat.Find(markerName);
  if (stat_idx < 0) {
    SNP_heterog_cond_stat.SetDouble(markerName, het);
  } else {
    double prev = SNP_heterog_cond_stat.Double(stat_idx);
    prev += het;
    SNP_heterog_cond_stat.SetDouble(markerName, prev);
  }

  // Increase degrees of freedom of het stat by 1 (each study adds 1 to the total df)
  int df_idx = SNP_heterog_cond_df.Find(markerName);
  if (df_idx < 0) {
    SNP_heterog_cond_df.SetInteger(markerName, 1);
  } else {
    int prev = SNP_heterog_cond_df.Integer(df_idx);
    prev += 1;
    SNP_heterog_cond_df.SetInteger(markerName, prev);
  }
}

/**
 * Update the heterogeneity statistic (Cochran's Q) for a variant given the score statistic in this study.
 *
 * This function adds its individual contribution to the running total for the heterogeneity statistic for this variant
 * over all studies.
 *
 * See https://git.io/JfLSM in METAL.
 *
 * @param study integer representing the study (index into the summaryFiles array of summary statistic files)
 * @param adjust bool representing whether to shift columns left by 1 depending on RAREMETAL or rvtest format
 * @param buffer string current line to be parsed
 * @param covReader extracts covariances for a variant if conditional analysis is required
 */
void Meta::poolHeterogeneity(int study, bool adjust, String &buffer, SummaryFileReader &covReader) {
  StringArray tokens;
  tokens.AddTokens(buffer, SEPARATORS);

  if (tokens[0].Find("#") != -1) {
    return;
  }

  if (tokens[0].Find("chr") != -1) {
    tokens[0] = tokens[0].SubStr(3);
  }

  String chr_pos = tokens[0] + ":" + tokens[1];
  int direction_idx = directionByChrPos.Integer(chr_pos);
  char direction = directions[direction_idx][study];

  if (direction == '?' || direction == '!') {
    // This variant was previously skipped due to QC or monomorphic
    return;
  }

  // poolSingleRecord() by this point has already figured out whether this variant's score statistic & effect size
  // need to be flipped to match the alleles used in the meta-analysis
  bool flip = flipSNP.Integer(String(std::to_string(study).c_str()) + ":" + chr_pos) != -1;

  // Get the score statistic and its variance
  double u_study = tokens[13 - adjust].AsDouble();
  double v_study = tokens[14 - adjust].AsDouble();

  // Update the heterogeneity statistic for this variant
  UpdateHetStats(study, chr_pos, u_study, v_study, flip);

  // If we're doing conditional analysis, then we also need to calculate the heterogeneity of the conditional meta-analysis score statistic
  if (cond != "" && v_study > 0.0) {
    //if this variant is not the one to be conditioned upon
    if (conditionVar.Integer(chr_pos) == -1) {
      UpdateHetCondStats(study, flip, adjust, chr_pos, tokens, covReader);
    }
  }
}

char Meta::GetDirection(String &chr_pos, double effsize, bool flip)
{
    char direction = '+';
    if (flip)
    {
        if (effsize > 0)
        {
            direction = '-';
        }
    } else
    {
        if (effsize < 0)
        {
            direction = '-';
        }
    }
    return direction;
}

int Meta::MatchTwoAlleles(String refalt_current, int &idx, String &chr_pos)
{
    //if matched an allele pair
    String refalt_hashed = refalt[idx];
    if (refalt_hashed == ".:.")
    {
        refalt[idx] = refalt_current;
        return 0;
    }
    StringArray tmp, tmp_current;
    tmp.AddTokens(refalt_hashed, ":");
    tmp_current.AddTokens(refalt_current, ":");
    if (tmp[1] == ".")
    {
        if (tmp_current[0] != tmp[0] && tmp_current[1] != tmp[0])
        {
            return 2;
        } else
        {
            refalt[idx] = refalt_current;
            return 0;
        }
    }
    if (refalt_current == refalt_hashed)
    {
        return 0;
    }
    if (tmp[1] + ":" + tmp[0] == refalt_current)
    {
        return 1;
    }
    //if previously saved alleles are from monomorphic sites
    if (tmp[1] == tmp[0])
    {
        StringArray tmp2;
        tmp2.AddTokens(refalt_current, ":");
        int idx = directionByChrPos.Integer(chr_pos);
        if (tmp2[0] == tmp[0])
        {
            refalt[idx] = refalt_current;
            return 0;
        }
        if (tmp2[1] == tmp[0])
        {
            refalt[idx] = tmp2[1] + ":" + tmp2[0];
            return 1;
        }
    }
    return 2;
}

int Meta::MatchOneAllele(String ref_current, int &idx)
{
    //if match ref allele return 0
    // if match alt allele return 1
    // if match nothing, return 2
    if (refalt[idx] == ".:.")
    {
        String name;
        name = ref_current + ":" + ref_current;
        refalt[idx] = name;
        return 0;
    }
    StringArray tmp;
    tmp.AddTokens(refalt[idx], ":");
    if (tmp[0] == "." || tmp[1] == ".")
    {
        if (ref_current == tmp[0])
        {
            return 0;
        }
        if (ref_current == tmp[1])
        {
            return 0;
        } else
        {
            return 2;
        }
    }

    if (tmp[1] == ref_current)
    {
        return 0;
    }
    if (tmp[0] == ref_current)
    {
        return 1;
    }
    return 2;
}


/**************************************************/


/*** pool summary stat related ****/

// pool single record from score file
// if u2/v2 is generated, return true
bool Meta::poolSingleRecord(int study, double &current_chisq, int &duplicateSNP, bool adjust, String &buffer,
                            SummaryFileReader &covReader)
{
    StringArray tokens;
    tokens.AddTokens(buffer, SEPARATORS);
    if (tokens[0].Find("#") != -1)
    {
        return 0;
    }

    if (tokens[0].Find("chr") != -1)
    {
        tokens[0] = tokens[0].SubStr(3);
    }

    String chr_pos = tokens[0] + ":" + tokens[1];
    int direction_idx = directionByChrPos.Integer(chr_pos);

    //CHECK duplicate markers
    bool is_dup = isDupMarker(tokens[0], chr_pos);
    if (is_dup)
    {
        duplicateSNP++;
        fprintf(log, "Warning: variant %s from study %s is skipped because of duplicate records in the same study.\n",
                chr_pos.c_str(), scorefile[study].c_str());
        return 0;
    }

    //POOLING STEP1: if fail HWE or CALLRATE then skip this record
    //if a variant has a missing allele but not monomorphic then exclude this variant without updating the total sample size
    //if(((tokens[2]=="." || tokens[3]==".") && (c1+c2!=0 && c2+c3!=0)) || (tokens[2]=="." && tokens[3]==".") || (tokens[2]==tokens[3] && c1+c2!=0 && c2+c3!=0))

    // An allele should never be zero, but this code checks for it anyway
    if (tokens[2] == "0")
    {
        tokens[2] = ".";
    }
    if (tokens[3] == "0")
    {
        tokens[3] = ".";
    }

    //check allele counts to see if the site is monomorphic
    double c1, c2, c3;
    bool use_dosage = false;
    if (dosageOptionFile != "")
    {
        use_dosage = dosageOptions[study];
    } else
    {
        use_dosage = dosage;
    }
    if (use_dosage)
    {
        c3 = tokens[4].AsDouble() * tokens[6].AsDouble() * tokens[6].AsDouble();
        c2 = tokens[4].AsDouble() * 2.0 * tokens[6].AsDouble() * (1.0 - tokens[6].AsDouble());
        c1 = tokens[4].AsDouble() * (1.0 - tokens[6].AsDouble()) * (1.0 - tokens[6].AsDouble());
    } else
    {
        c1 = tokens[10 - adjust].AsDouble();
        c2 = tokens[11 - adjust].AsDouble();
        c3 = tokens[12 - adjust].AsDouble();
    }

    double current_AC;
    int current_N;
    current_N = c1 + c2 + c3;
    current_AC = 2 * c3 + c2;
    String refalt_current = tokens[2] + ':' + tokens[3];
    String markername = chr_pos + ":" + refalt_current;
    if (useExactMetaMethod)
    {
        String snp_no_allele = tokens[0] + ":" + tokens[1];
        double raw_af;
        if (current_N == 0)
        {
            printf("Warning: effective sample size of %s:%s is zero. Are you using the correct --dosage or --dosageOptionFile?\n\n",
                   tokens[0].c_str(), tokens[1].c_str());
            raw_af = 0;
        } else
        {
            raw_af = (double) current_AC / current_N / 2;
        }
        if (normPop)
        { // adjust for population stratification by regression
            double fk = getAFtilda(markername, raw_af, study);
            addToMapStrVec(variant_fk, study, chr_pos, fk, scorefile.Length() + 1, 0);
            addToMapStrVec(variant_nk, study, chr_pos, (double) current_N, scorefile.Length() + 1, current_N);
        } else
        {
            addToMapStrVec(variant_fk, study, chr_pos, raw_af, scorefile.Length() + 1, 0);
            addToMapStrVec(variant_nk, study, chr_pos, (double) current_N, scorefile.Length() + 1, current_N);
        }
    }

    bool is_fail = 0;
    int filter_type = 1;
    if ((tokens[2] == "." && tokens[3] == ".") || (tokens[2] == tokens[3] && c1 + c2 != 0 && c2 + c3 != 0))
    {
        is_fail = 1;
    }
    if (tokens[8 - adjust].AsDouble() < CALLRATE || tokens[9 - adjust].AsDouble() < HWE)
    {
        // hwe & call rate filter
        is_fail = 1;
        filter_type = 0;
    }

    if (is_fail)
    {
        char direct = '!';
        UpdateDirection(direction_idx, study, direct, chr_pos, true);
        UpdateStrIntHash(chr_pos, SampleSize[study], usefulSize);
        UpdateExcludedMarker(study, chr_pos, filter_type,
                             tokens[0] + ":" + tokens[1] + ":" + tokens[2] + ":" + tokens[3]);
        fprintf(log,"Warning: variant %s from study %s is skipped because of call rate, HWE, or missing alleles.\n",chr_pos.c_str(),scorefile[study].c_str());
        return 0;
    }

    //STEP2: check if this position has been hashed. If yes, match alleles; if not, hash position and ref alt alleles.
    int marker_idx = directionByChrPos.Integer(chr_pos);
    bool flip = false;
    double u = tokens[13 - adjust].AsDouble();
    double v = tokens[14 - adjust].AsDouble();
    if (relateBinary)
    {
        u *= Const_binary[study];
        v *= Const_binary[study];
    }
    double raw_v = v;
    double raw_u = u;
    if (marker_idx == -1) //if this position is never hashed before
    {
        UpdateStrIntHash(chr_pos, current_N, recSize); // add size to --altMAF when site is observed
        if (c2 + c3 == 0 || c1 + c2 == 0)
        { //if monomorphic, then hash in ref=alt allele for this position and update direction with '?' usefulAC as count of ref allele
            int count = 0;
            setRefAltHashKey(refalt_current, tokens, c1, c2);
            refalt.Push(refalt_current);
            char direct = '?';
            UpdateDirection(direction_idx, study, direct, chr_pos, false);
            //should only hash the count of alternative allele
            UpdateACInfo(chr_pos, count);
            if (SampleSize[study] != current_N)
            {
                UpdateStrIntHash(chr_pos, SampleSize[study] - current_N, usefulSize);
            }
        } else
        { //if not monomorphic, then hash in ref:alt allele and update direction with alt direction, update usefulSize if count is different from sampleSize, update usefulAC with alt allele count.
            refalt.Push(refalt_current);
            char direct = GetDirection(chr_pos, tokens[15 - adjust].AsDouble(), false);
            UpdateDirection(direction_idx, study, direct, chr_pos, false);
            UpdateACInfo(chr_pos, current_AC);
            if (SampleSize[study] != current_N)
            {
                UpdateStrIntHash(chr_pos, SampleSize[study] - current_N, usefulSize);
            }
            UpdateStats(study, chr_pos, u, v, flip);
        }
    } else //if this position has been seen from previous studies
    {
        int match;
        //if the previous records are all monomorphic
        if (c2 + c3 == 0 || c1 + c2 == 0)
        {
            //if monomorphic, then check the following: ref?=alt,if yes, check if ref of this study is the same as the one of the alleles previously hashed for this position (if not the same, then exclude marker;if they are the same update counts and statistics).
            int count = 0;
            setPolyMatch(match, chr_pos, refalt_current, tokens, marker_idx);
            // set count from match
            if (tokens[2] == "." || tokens[3] == ".")
            {
                if (match == 1)
                {
                    count = 0;
                } else if (match == 0)
                {
                    count = 0;
                }
            } else
            {
                if (match == 0)
                {
                    count = 2 * c3 + c2;
                } else if (match == 1)
                {
                    count = 2 * c1 + c2;
                }
            }

            if (match == 2)
            { //if allels do not match, then exclude this variant and continue
                UpdateExcludedMarker(study, chr_pos, 2,
                                     tokens[0] + ":" + tokens[1] + ":" + tokens[2] + ":" + tokens[3]);
                UpdateStrIntHash(chr_pos, SampleSize[study], usefulSize);
                char direct = '!';
                UpdateDirection(direction_idx, study, direct, chr_pos, true);
                fprintf(log,"Warning: variant %s from study %s is skipped because alleles did not match previous study.\n",chr_pos.c_str(),scorefile[study].c_str());
                return 0;
            }
            if (match == 1)
            {
                flip = true;
            }

            UpdateStrIntHash(chr_pos, current_N, recSize); // add --altMAF
            UpdateACInfo(chr_pos, count);
            if (SampleSize[study] != current_N)
            {
                UpdateStrIntHash(chr_pos, SampleSize[study] - current_N, usefulSize);
            }
            char direct = '?';
            UpdateDirection(direction_idx, study, direct, chr_pos, false);
        } else //if site is not monomorphic,then need to match alleles and upate statistics
        {
            match = MatchTwoAlleles(refalt_current, marker_idx, chr_pos);
            if (match == 2)
            {
                UpdateExcludedMarker(study, chr_pos, 2,
                                     tokens[0] + ":" + tokens[1] + ":" + tokens[2] + ":" + tokens[3]);
                UpdateStrIntHash(chr_pos, SampleSize[study], usefulSize);
                char direct = '!';
                UpdateDirection(direction_idx, study, direct, chr_pos, true);
                fprintf(log,"Warning: variant %s from study %s is skipped because alleles did not match previous study.\n",chr_pos.c_str(),scorefile[study].c_str());
                return 0;
            }
            UpdateStrIntHash(chr_pos, current_N, recSize);
            if (match == 1)
            {
                flip = true;
                if (SampleSize[study] != current_N)
                {
                    UpdateStrIntHash(chr_pos, SampleSize[study] - current_N, usefulSize);
                }
                UpdateACInfo(chr_pos, 2 * c1 + c2);
                char direct = GetDirection(chr_pos, tokens[15 - adjust].AsDouble(), true);
                UpdateDirection(direction_idx, study, direct, chr_pos, false);
                UpdateStats(study, chr_pos, u, v, flip);
            }
            if (match == 0)
            {
                if (SampleSize[study] != current_N)
                {
                    UpdateStrIntHash(chr_pos, SampleSize[study] - current_N, usefulSize);
                }
                UpdateACInfo(chr_pos, current_AC);
                char direct = GetDirection(chr_pos, tokens[15 - adjust].AsDouble(), false);
                UpdateDirection(direction_idx, study, direct, chr_pos, false);
                UpdateStats(study, chr_pos, u, v, flip);
            }
        }
    }
    if (useExactMetaMethod)
    { // get it back for calculating gc
        v = raw_v;
        u = raw_u;
    }

    if (flip)
    {
        flipCount++;
        flipSNP.SetInteger(String(std::to_string(study).c_str()) + ":" + chr_pos, flipCount);
    }

    if (sumCaseAC && tokens.Length() == 19)
    {
        if (caseAC.Integer(chr_pos) == -1)
        {
            if (flip)
            {
                caseAC.SetInteger(chr_pos, current_N * 2 - tokens[17].AsInteger());
                controlAC.SetInteger(chr_pos, current_N * 2 - tokens[18].AsInteger());
            } else
            {
                caseAC.SetInteger(chr_pos, tokens[17].AsInteger());
                controlAC.SetInteger(chr_pos, tokens[18].AsInteger());
            }
        } else
        {
            int c1 = caseAC.Integer(chr_pos);
            int c2 = controlAC.Integer(chr_pos);
            if (flip)
            {
                caseAC.SetInteger(chr_pos, c1 + current_N * 2 - tokens[17].AsInteger());
                controlAC.SetInteger(chr_pos, c2 + current_N * 2 - tokens[18].AsInteger());
            } else
            {
                caseAC.SetInteger(chr_pos, c1 + tokens[17].AsInteger());
                controlAC.SetInteger(chr_pos, c2 + tokens[18].AsInteger());
            }
        }
    }

    //if a variant is monomorphic then update the count of sample size and generate warning
    if (c1 + c2 == 0 || c2 + c3 == 0)
    {
        fprintf(log,"Warning: variant %s from study %s is skipped because variant is monomorphic.\n",chr_pos.c_str(),scorefile[study].c_str());
        return 0;
    }

    //push the chisq statistics for GC calculation
    if (tokens[14 - adjust].AsDouble() > 0)
    {
        current_chisq = u * u / (v * v);
    }

    //update SNP_cond_stat and SNP_cond_V
    //if(cond!="" && cond_status[study] && tokens[14-adjust].AsDouble()>0.0)
    if (cond != "" && tokens[14 - adjust].AsDouble() > 0.0)
    {
        //if this variant is not the one to be conditioned upon
        if (conditionVar.Integer(tokens[0] + ":" + tokens[1]) == -1)
        {
            updateSNPcond(study, flip, adjust, chr_pos, tokens, covReader);
        }
    }
    return 1;
}

// if dup marker, return true
bool Meta::isDupMarker(String &chr_str, String &chr_pos)
{
    bool is_dup = 0;
    if (chr_str != dupcheck_chr)
    {
        hashToRemoveDuplicates.Clear();
        dupcheck_chr = chr_str;
    }
    if (hashToRemoveDuplicates.Integer(chr_pos) == -1)
    {
        varCount++;
        hashToRemoveDuplicates.SetInteger(chr_pos, varCount);
    } else
    {
        is_dup = 1;
    }
    return is_dup;
}

void Meta::setRefAltHashKey(String &refalt_current, StringArray &tokens, int c1, int c2)
{
    if (tokens[2] == "." || tokens[3] == ".")
    {
        if (tokens[2] == ".")
        {
            refalt_current = tokens[3] + ":.";
        } else
        {
            refalt_current = tokens[2] + ":.";
        }
    } else if (tokens[2] == tokens[3])
    {
        refalt_current = tokens[2] + ":.";
    } else
    {
        refalt_current = tokens[2] + ':' + tokens[3];
        if (c1 + c2 == 0)
        {
            refalt_current = tokens[3] + ':' + tokens[2];
        }
    }
}

void Meta::setPolyMatch(int &match, String &chr_pos, String &refalt_current, StringArray &tokens, int marker_idx)
{
    if (tokens[2] != tokens[3])
    {
        if (tokens[2] == ".")
        {
            match = MatchOneAllele(tokens[3], marker_idx);
        } else if (tokens[3] == ".")
        {
            match = MatchOneAllele(tokens[2], marker_idx);
        } else
        {
            match = MatchTwoAlleles(refalt_current, marker_idx, chr_pos);
        }
    } else
    {
        match = MatchOneAllele(tokens[3], marker_idx);
    }
}

//update SNP_cond_stat and SNP_cond_V
void Meta::updateSNPcond(int study, bool flip, int adjust, String &chr_pos, StringArray &tokens,
                         SummaryFileReader &covReader)
{
    Vector GX;
    CalculateGenotypeCov(covReader, tokens[0], tokens[1].AsInteger(), study, GX);
    double cond_u = tokens[13 - adjust].AsDouble() - GX.InnerProduct(commonVar_betaHat[study]);
    int idx;
    idx = SNPstat_cond.Find(chr_pos);
    double stat = 0.0;
    if (idx != -1)
    {
        stat = SNPstat_cond.Double(idx);
        if (!flip)
        {
            stat += cond_u;
        } else
        {
            stat += -1.0 * cond_u;
        }
        SNPstat_cond.SetDouble(idx, stat);
    } else
    {
        if (!flip)
        {
            stat += cond_u;
        } else
        {
            stat += -1.0 * cond_u;
        }
        SNPstat_cond.SetDouble(chr_pos, stat);
    }

    Vector tmp;
    for (int i = 0; i < GX.dim; i++)
    {
        tmp.Push(GX.InnerProduct(XX_inv[study][i]));
    }
    double v = tokens[14 - adjust].AsDouble();
    double cond_v_part = tmp.InnerProduct(GX);
    double cond_V = v * v - cond_v_part;
    idx = SNP_Vstat_cond.Find(chr_pos);
    if (idx != -1)
    {

        double cond_Vstat = SNP_Vstat_cond.Double(idx);
        cond_Vstat += cond_V;
        SNP_Vstat_cond.SetDouble(idx, cond_Vstat);
    } else
    {
        SNP_Vstat_cond.SetDouble(chr_pos, cond_V);
    }
}


// after reading & pooling stat from score file, set overall af of each marker
// results stored in SNPmaf
void Meta::setPooledAF()
{
    StringArray chr_AC, unique_chr, SNPname_AC;
    IntArray pos_AC;
    //get the unique chromosomes
    for (int i = 0; i < directionByChrPos.Capacity(); i++)
    {
        if (!directionByChrPos.SlotInUse(i))
        {
            continue;
        }
        String SNPname = directionByChrPos[i];
        int idx = directionByChrPos.Integer(SNPname);
        StringArray tmp;
        tmp.AddTokens(SNPname, ":");
        chr_AC.Push(tmp[0]);
        pos_AC.Push(tmp[1].AsInteger());
        if (directions[idx].Length() < scorefile.Length())
        {
            for (int l = directions[idx].Length(); l < scorefile.Length(); l++)
            {
                directions[idx] += '?';
            }
        }
        SNPname += ':';
        SNPname += refalt[idx];
        SNPname_AC.Push(SNPname);
        if (unique_chr.Find(tmp[0]) == -1)
        {
            unique_chr.Push(tmp[0]);
        }
    }
    QuickIndex chr_AC_idx(chr_AC);
    unique_chr.Sort();
    StringArray chr_cp, character_chr;

    // pool
    for (int i = 0; i < unique_chr.Length(); i++)
    {
        if (unique_chr[i].AsInteger() <= 22 && unique_chr[i].AsInteger() >= 1)
        {
            chr_cp.Push(unique_chr[i]);
        } else
        {
            character_chr.Push(unique_chr[i]);
        }
    }
    for (int i = 0; i < character_chr.Length(); i++)
    {
        chr_cp.Push(character_chr[i]);
    }
    unique_chr = chr_cp; //now unique_chr are sorted as 1,2,...,22,X,Y,M...
    chr_cp.Clear();
    character_chr.Clear();
    for (int i = 0; i < unique_chr.Length(); i++)
    {
        IntArray pos_i;
        StringArray SNPname_i;
        for (int j = 0; j < chr_AC.Length(); j++)
        {
            if (chr_AC[chr_AC_idx[j]] == unique_chr[i])
            {
                pos_i.Push(pos_AC[chr_AC_idx[j]]);
                SNPname_i.Push(SNPname_AC[chr_AC_idx[j]]);
            }
        }
        QuickIndex pos_i_idx(pos_i);
        for (int j = 0; j < pos_i.Length(); j++)
        {
            StringArray tmp;
            tmp.AddTokens(SNPname_i[pos_i_idx[j]], ":");
            double AC = usefulAC.Double(tmp[0] + ":" + tmp[1]);

            // usefulSize has the # of samples to be excluded
            // recSize has # of samples truly added from vcf/ped. Use this when altMAF is toggled
            int N;
            if (this->altMAF)
            {
                N = recSize.Integer(tmp[0] + ":" + tmp[1]);
                if (N == -1)
                {
                    N = 0;
                }
            } else
            { // use default
                N = usefulSize.Integer(tmp[0] + ":" + tmp[1]);
                if (N != -1)
                {
                    N = total_N - N;
                } else
                {
                    N = total_N;
                }
            }

            double maf;
            if (founderAF)
            {
                maf = AC / (2.0 * N);
            } else
            {
                maf = AC / (2.0 * N);
            }

            int idx = directionByChrPos.Integer(tmp[0] + ":" + tmp[1]);
            if (directions[idx].FindChar('+') == -1 && directions[idx].FindChar('-') == -1)
            {
                maf = 0.0;
            }

            SNPmaf_maf.Push(maf);
            SNPmaf_name.Push(SNPname_i[pos_i_idx[j]]);
            SNP_effect_N.Push(N);
            SNPmaf.SetInteger(SNPname_i[pos_i_idx[j]], SNPmaf_maf.Length() - 1);
        }
    }
}


void Meta::printSingleMetaHeader(String &filename, IFILE &output)
{
    ifprintf(output, "##Method=SinglevarScore\n");
    ifprintf(output, "##STUDY_NUM=%d\n", scorefile.Length());
    ifprintf(output, "##TotalSampleSize=%d\n", total_N);

    String header = "#CHROM\tPOS\tREF\tALT\tN\tPOOLED_ALT_AF\tDIRECTION_BY_STUDY\tEFFECT_SIZE\tEFFECT_SIZE_SD\tH2";
    header += logP ? "\tLOG_PVALUE" : "\tPVALUE";

    if (cond != "") {
      header += "\tCOND_EFFSIZE\tCOND_EFFSIZE_SD\tCOND_H2\tCOND_PVALUE";
      header += logP ? "\tCOND_LOG_PVALUE" : "\tCOND_PVALUE";
    }

    if (sumCaseAC) {
      header += "\tsumCaseAC\tsumControlAC";
    }

    if (bHeterogeneity) {
      header += "\tHET_I2\tHET_CHISQ";
      header += logP ? "\tHET_LOG_PVALUE" : "\tHET_PVALUE";
    }

    if (bHeterogeneity && cond != "") {
      header += "\tHET_COND_I2\tHET_COND_CHISQ";
      header += logP ? "\tHET_COND_LOG_PVALUE" : "\tHET_COND_PVALUE";
    }

    ifprintf(output, header.c_str());
    ifprintf(output, "\n");
}

void Meta::printOutVcfHeader(String &vcf_filename, IFILE &vcfout)
{
    if (prefix == "")
    {
        vcf_filename = "pooled.variants.vcf";
    } else if (prefix.Last() == '.' || prefix.Last() == '/')
    {
        vcf_filename = prefix + "pooled.variants.vcf";
    } else
    {
        vcf_filename = prefix + ".pooled.variants.vcf";
    }

    vcfout = ifopen(vcf_filename, "w", InputFile::UNCOMPRESSED);
    ifprintf(vcfout, "#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}


// print ith record  of single var result
void Meta::printSingleMetaVariant(GroupFromAnnotation &group, int i, IFILE &output, IFILE &vcfout)
{
    if (Nsamples == -1)
    {
        for (int i = 0; i < SNP_effect_N.Length(); i++)
        {
            if (SNP_effect_N[i] > Nsamples)
            {
                Nsamples = SNP_effect_N[i];
            }
        }
        if (Nsamples == -1)
        {
            error("No available SNPs to analysis!");
        }
    }

    String SNPname = SNPmaf_name[i];
    StringArray tmp;
    tmp.AddTokens(SNPname, ":");
    String SNPname_noallele = tmp[0] + ":" + tmp[1];
    int N = SNP_effect_N[i];
    double U = SNPstat.Double(SNPname_noallele);
    double V = SNP_Vstat.Double(SNPname_noallele);
    double maf = SNPmaf_maf[i];
    if (useExactMetaMethod)
    {
        std::map<String, std::vector<double> >::iterator pfk = variant_fk.find(SNPname_noallele);
        std::map<String, std::vector<double> >::iterator pnk = variant_nk.find(SNPname_noallele);
        if (pfk == variant_fk.end())
        {
            error("[Meta::printSingleMetaVariant] Cannot find %s!\n", SNPname_noallele.c_str());
        }
        double exact_maf;
        if (normPop)
        {
            // fill in variant_f
            double f = 0;
            for (int i = 0; i < scorefile.Length(); i++)
            {
                f += pfk->second[i] / N * pnk->second[i];
            }
            exact_maf = f;
        } else
        {
            exact_maf = maf;
        }
        // here V = sum(vk / sigmak)
        double nkdeltak = 0;
        for (int i = 0; i < scorefile.Length(); i++)
        {
            nkdeltak += pnk->second[i] * Ydelta[i];
        }
        double nkdeltakfk = 0;
        for (int i = 0; i < scorefile.Length(); i++)
        {
            nkdeltakfk += pnk->second[i] * pfk->second[i] * Ydelta[i];
        }
        U = U - 2 * exact_maf * nkdeltak + 2 * nkdeltakfk;
        double v2 = 0;
        for (int i = 0; i < scorefile.Length(); i++)
        {
            v2 += pnk->second[i] * pfk->second[i] * pfk->second[i];
        }
        double new_r = 0;
        for (int i = 0; i < scorefile.Length(); i++)
        {
            new_r += (pnk->second[i] - 1) * Ysigma2[i] + pnk->second[i] * Ydelta[i] * Ydelta[i];
        }
        new_r /= (N - 1);
        V = new_r * (V + 4 * v2 - 4 * N * exact_maf * exact_maf);
        // update U and V for gene based test
        pfk->second[scorefile.Length()] = exact_maf;
        pnk->second[scorefile.Length()] = new_r;
        SNPstat.SetDouble(SNPname_noallele, U);
        SNP_Vstat.SetDouble(SNPname_noallele, V);
    }
    int direction_idx = directionByChrPos.Integer(SNPname_noallele);

    String direction = directions[direction_idx];
    IntArray pvalue1_idx, pvalue5_idx;


    if (maf == 0.0 || maf == 1.0)
    {
        if (cond != "")
        {
            ifprintf(output, "%s\t%s\t%s\t%s\t%d\t%g\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", tmp[0].c_str(),
                     tmp[1].c_str(), tmp[2].c_str(), tmp[3].c_str(), N, maf, direction.c_str());
        } else
        {
            ifprintf(output, "%s\t%s\t%s\t%s\t%d\t%g\t%s\tNA\tNA\tNA\tNA\n", tmp[0].c_str(), tmp[1].c_str(),
                     tmp[2].c_str(), tmp[3].c_str(), N, maf, direction.c_str());
        }
    }

    if (maf <= 0.0 || maf >= 1.0 || V == 0.0)
    {
        return;
    }

    SingleVariantResult result(U, V, N);

    singleVarPvalue.SetDouble(SNPname, result.pvalue);
    singleVarEff.SetDouble(SNPname, result.effSize);

    pvalueAll.Push(result.pvalue);
    if (maf < 0.01)
    {
        pvalue1.Push(result.pvalue);
        pvalue1_idx.Push(pvalueAll.Length() - 1);
    }
    if (maf < 0.05)
    {
        pvalue5.Push(result.pvalue);
        pvalue5_idx.Push(pvalueAll.Length() - 1);
    }
    chr_plot.Push(tmp[0]);
    pos_plot.Push(tmp[1].AsInteger());

    SingleVariantResult cond_result;
    if (cond != "")
    { // print conditional analysis results
        double cond_U, cond_V;

        if (conditionVar.Integer(SNPname_noallele) == -1)
        {
            cond_U = SNPstat_cond.Double(SNPname_noallele);
            cond_V = SNP_Vstat_cond.Double(SNPname_noallele);
            cond_result.setStats(cond_U, cond_V, N);
            cond_result.calculate();
        } else
        {
            cond_result.effSize = 0.0;
            cond_result.pvalue = 1.0;
            cond_result.h2 = 0.0;
            cond_result.effSize_se = 0.0;
        }
        ifprintf(output, "%s\t%s\t%s\t%s\t%d\t%g\t%s\t%g\t%g\t%g\t%s%g\t%g\t%g\t%g\t%s%g", tmp[0].c_str(),
                 tmp[1].c_str(), tmp[2].c_str(), tmp[3].c_str(), N, maf, direction.c_str(), result.effSize, result.effSize_se, result.h2,
                 result.disect ? "<" : "", logP ? result.log_pvalue : result.pvalue, cond_result.effSize, cond_result.effSize_se, cond_result.h2, cond_result.disect ? "<" : "",
                 logP ? cond_result.log_pvalue : cond_result.pvalue);

        pvalueAll_cond.Push(cond_result.pvalue);
        if (maf < 0.01)
        {
            pvalue1_cond.Push(cond_result.pvalue);
        }
        if (maf < 0.05)
        {
            pvalue5_cond.Push(cond_result.pvalue);
        }
    } else
    {
        ifprintf(output, "%s\t%s\t%s\t%s\t%d\t%g\t%s\t%g\t%g\t%g\t%s%g", tmp[0].c_str(), tmp[1].c_str(), tmp[2].c_str(),
                 tmp[3].c_str(), N, maf, direction.c_str(), result.effSize, result.effSize_se, result.h2, result.disect ? "<" : "", logP ? result.log_pvalue : result.pvalue);
    }

    if (sumCaseAC)
    {
        int c1 = caseAC.Integer(SNPname_noallele);
        int c2 = controlAC.Integer(SNPname_noallele);
        if (c1 == -1 || c2 == -1)
        {
            printf("Warning: variant doesn't have caseAC or controlAC count!\n");
        }
        ifprintf(output, "\t%d\t%d", c1, c2);
    }

    if (bHeterogeneity) {
      /**
       * Calculate heterogeneity p-value and I2 statistic.
       * See https://git.io/JfLS9 in METAL.
       */
      double het_stat = SNP_heterog_stat.Double(SNPname_noallele);
      int het_df = SNP_heterog_df.Integer(SNPname_noallele);
      double I2 = (het_stat <= het_df - 1) || (het_df <= 1) ? 0.0 : (het_stat - het_df + 1) / het_stat * 100.0;

      double het_pval;
      if (logP) {
        het_pval = (het_stat < 1e-7) || (het_df <= 1) ? 0.0 : -pchisq(het_stat, het_df - 1, 0, 1) / LN_10;
      }
      else {
        het_pval = (het_stat < 1e-7) || (het_df <= 1) ? 1.0 : pchisq(het_stat, het_df - 1, 0, 0);
      }

      ifprintf(output, "\t%g\t%g\t%g", I2, het_stat, het_pval);
    }

    if (bHeterogeneity && cond != "") {
      double cond_het_stat = SNP_heterog_cond_stat.Double(SNPname_noallele);
      int cond_het_df = SNP_heterog_cond_df.Integer(SNPname_noallele);
      double cond_I2 = (cond_het_stat <= cond_het_df - 1) || (cond_het_df <= 1) ? 0.0 : (cond_het_stat - cond_het_df + 1) / cond_het_stat * 100.0;

      double cond_het_pval;
      if (logP) {
        cond_het_pval = (cond_het_stat < 1e-7) || (cond_het_df <= 1) ? 0.0 : -pchisq(cond_het_stat, cond_het_df - 1, 0, 1) / LN_10;
      }
      else {
        cond_het_pval = (cond_het_stat < 1e-7) || (cond_het_df <= 1) ? 1.0 : pchisq(cond_het_stat, cond_het_df - 1, 0, 0);
      }

      ifprintf(output, "\t%g\t%g\t%g", cond_I2, cond_het_stat, cond_het_pval);
    }

    ifprintf(output, "\n");

    if (outvcf)
    {
        ifprintf(vcfout, "%s\t%s\t%s\t%s\t%s\t.\t.\tALT_AF=;\n", tmp[0].c_str(), tmp[1].c_str(), tmp[1].c_str(),
                 tmp[2].c_str(), tmp[3].c_str());
    }

    // annotation if needed
    if (group.labelHits)
    {
        annotateSingleVariantToGene(group, result.pvalue, cond_result.pvalue, tmp);
    }
}

//Annotate single variants to gene
void Meta::annotateSingleVariantToGene(GroupFromAnnotation &group, double pvalue, double cond_pvalue, StringArray &tmp)
{
    // initialize
    if (geneLabel.Length() == 0)
    {
        target_chr = "";
        target_pos = 0;
        target = 0;
        target_pvalue = _NAN_;
    }

    bool skip = 0;
    if (pvalue > 0.05 / SNPmaf_maf.Length())
    {
        if (cond != "")
        {
            if (cond_pvalue > 0.05 / SNPmaf_maf.Length())
            {
                skip = 1;
            }
        } else
        {
            skip = 1;
        }
    }
    if (skip)
    {
        geneLabel.Push("");
        return;
    }

    // only annotate those significant signals
    if (target_chr == tmp[0] && target_pos > tmp[1].AsInteger() - 1000000)
    { //annotated already
        if (pvalue < target_pvalue)
        { //if this is the higher peak then update the annotation of this locus
            String current_anno = group.AnnotateSingleVar(tmp[0], tmp[1].AsInteger());
            int distance = 1000;
            while (current_anno == "" && distance <= 1000000)
            {
                current_anno = group.AnnotateSingleVar(tmp[0], tmp[1].AsInteger() + distance);
                if (current_anno == "")
                {
                    current_anno = group.AnnotateSingleVar(tmp[0], tmp[1].AsInteger() - distance);
                }
                distance += 1000;
            }
            if (geneLabel[target] != "" && current_anno != "" && geneLabel[target].Find(current_anno) == -1)
            {
                current_anno = geneLabel[target] + "/" + current_anno;
            }
            if (geneLabel[target] != "" && current_anno == "")
            {
                current_anno = geneLabel[target];
            }
            geneLabel.Push(current_anno);
            geneLabel[target] = "";
            target = geneLabel.Length() - 1;
            target_pvalue = pvalue;
            target_pos = tmp[1].AsInteger();
            target_chr = tmp[0];
        } else
        { //otherwise, leave the original annotation of this locus
            geneLabel.Push("");
        }
    } else
    {
        geneLabel.Push(group.AnnotateSingleVar(tmp[0], tmp[1].AsInteger()));
        target_chr = tmp[0];
        target_pos = tmp[1].AsInteger();
        target_pvalue = pvalue;
        target = geneLabel.Length() - 1;
    }
}

// single variant GC plot by MAF
void Meta::plotSingleMetaGC(IFILE &output, bool calc_gc)
{
    String title, demo1, demo2, demo3;
    title = "single variant analysis";
    if (pvalueAll.size <= 0)
    {
        error("No available p values. Something goes wrong?\n");
    }
    double GC1 = GetGenomicControlFromPvalue(pvalueAll);
    demo1 = "GC=";
    demo1 += GC1;
    double GC2 = 0;
    double GC3 = 0;
    if (pvalue1.size <= 0)
    {
        printf("\nWarning: no MAF<0.01 variants. This group GC = 0!\n");
    } else
    {
        GC2 = GetGenomicControlFromPvalue(pvalue1);
    }
    demo2 = "GC=";
    demo2 += GC2;
    if (pvalue5.size <= 0)
    {
        printf("\nWarning: no MAF<0.05 variants. This group GC = 0!\n");
    } else
    {
        GC3 = GetGenomicControlFromPvalue(pvalue5);
    }
    demo3 = "GC=";
    demo3 += GC3;
    writepdf.Draw(pdf, geneLabel, pvalueAll, pvalue1, pvalue5, chr_plot, pos_plot, title, demo1, demo2, demo3, false);

    //Calculate genomic control
    if (calc_gc)
    {
        ifprintf(output, "#Genomic Control for pooled sample is: %g\n", GC1);
        printf("  Genomic Control for all studies are listed below:\n");
        for (int s = 0; s < scorefile.Length(); s++)
        {
            printf("  %g\t", GCbyStudy[s]);
            ifprintf(output, "#Genomic Control for study %d is: %g\n", s, GCbyStudy[s]);
        }
    }
}

/****************************************************************************************/


/*********** for Meta::Run **************/
void Meta::loadSingleStatsInGroup(GroupFromAnnotation &group)
{
    for (int g = 0; g < group.annoGroups.Length(); g++)
    {
        IntArray del; //del has the SNPs to be deleted from SNPlist
        maf[g].Dimension(0);
        stats[g].Dimension(0);
        int count = group.SNPlist[g].Length();
        singlePvalue[g].Dimension(0);
        singleEffSize[g].Dimension(0);

        for (int m = 0; m < count; m++)
        {
            String newSNP;
            bool flipStatus = false;
            double af = 0.0;
            double singleP = singleVarPvalue.Double(group.SNPlist[g][m]);

            if (singleP != _NAN_)
            {
                af = SNPmaf_maf[SNPmaf.Integer(group.SNPlist[g][m])];
            }

            double singleEff = singleVarEff.Double(group.SNPlist[g][m]);
            if (singleP == _NAN_)
            {
                flipStatus = true;
                StringArray SNPflip;
                SNPflip.AddTokens(group.SNPlist[g][m], ":");
                newSNP = SNPflip[0] + ":" + SNPflip[1] + ":" + SNPflip[3] + ":" + SNPflip[2];
                singleP = singleVarPvalue.Double(newSNP);

                if (singleP != _NAN_)
                {
                    af = SNPmaf_maf[SNPmaf.Integer(newSNP)];
                }

                singleEff = singleVarEff.Double(newSNP);
            }

            if (singleP == _NAN_)
            {
                del.Push(m);
                fprintf(log,
                        "Warning: variant %s is excluded from group %s for the following reasons:monomorphic or failed QC.\n",
                        group.SNPlist[g][m].c_str(), group.annoGroups[g].c_str());
                continue;
            }
            if (af == 0.0 || af == 1.0)
            {
                del.Push(m);
                fprintf(log, "Warning: variant %s is excluded from group %s for the following reasons: monomorphic.\n",
                        group.SNPlist[g][m].c_str(), group.annoGroups[g].c_str());
                continue;
            }

            if (af > 0.5 && 1.0 - af > MAF_cutoff)
            {
                del.Push(m);
                fprintf(log, "Warning: variant %s is excluded from group %s for the following reason: maf>cutoff.\n",
                        group.SNPlist[g][m].c_str(), group.annoGroups[g].c_str());
                continue;
            }

            if (af <= 0.5 && af > MAF_cutoff)
            {
                del.Push(m);
                fprintf(log, "Warning: variant %s is excluded from group %s for the following reason: maf>cutoff.\n",
                        group.SNPlist[g][m].c_str(), group.annoGroups[g].c_str());
                continue;
            }

            double tmp;
            StringArray tmp_name;
            tmp_name.AddTokens(group.SNPlist[g][m], ":");
            if (flipStatus)
            {
                group.SNPlist[g][m] = newSNP;
                StringArray tmp_name;
                tmp_name.AddTokens(newSNP, ":");
                tmp = SNPstat.Double(tmp_name[0] + ":" + tmp_name[1]);
            } else
            {
                tmp = SNPstat.Double(tmp_name[0] + ":" + tmp_name[1]);
            }
            if (af > 0.5)
            {
                tmp *= -1.0;
                af = 1.0 - af;
                singleEff *= -1.0;
                StringArray tokens;
                tokens.AddTokens(group.SNPlist[g][m], ":_/");
                group.SNPlist[g][m] = tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" + tokens[2];
                //update flipSNP to get the covaraince right
                for (int s = 0; s < scorefile.Length(); s++)
                {
                    String markername = String(std::to_string(s).c_str()) + ":" + tokens[0] + ":" + tokens[1];
                    if (flipSNP.Integer(markername) != -1)
                    {
                        flipSNP.Delete(markername);
                    } else
                    {
                        flipCount++;
                        flipSNP.SetInteger(markername, flipCount);
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
        for (int i = del.Length() - 1; i >= 0; i--)
        {
            group.SNPlist[g].Delete(del[i]);
            group.SNPNoAllele[g].Delete(del[i]);
        }
        count = group.SNPlist[g].Length();
        cov[g].Dimension(count, count, 0.0);
        if (cond != "")
        {
            cond_stats[g].Dimension(count, 0.0);
            cond_cov[g].Dimension(count, count, 0.0);
        }
    }
}

//loop through cov matrices of all studies and update cov
void Meta::loadSingleCovInGroup(GroupFromAnnotation &group)
{
    if (simplifyCovLoad)
    { // only load markers in group
        for (int g = 0; g < group.annoGroups.Length(); g++)
        {
            int gvar_count = group.SNPlist[g].Length();
            for (int i = 0; i < gvar_count; i++)
            {
                groupAnchor[group.SNPNoAllele[g][i]] = 0;
            }
        }
    }

    for (int study = 0; study < covfile.Length(); study++)
    {
        SummaryFileReader covReader;
        covReader.marker_col = this->marker_col;
        covReader.cov_col = this->cov_col;

        String covFilename = covfile[study];
        setFromRvOrRmwAdjust(FormatAdjust[study], marker_col, cov_col);
        printf("Reading cov matrix from study %d ...\n", study + 1);
        String filename = covfile[study];
        IFILE covfile_;
        covfile_ = ifopen(filename, "r");
        if (covfile_ == NULL)
        {
            error("ERROR! Cannot open file: %s! Input cov file has to be bgzipped and tabix indexed using the following command:\n bgzip yourfile.singlevar.cov.txt; tabix -c \"#\" -s 1 -b 2 -e 2 yourprefix.singlevar.cov.txt.gz\n",
                  filename.c_str());
        }
        Tabix covtabix;
        String tabix_name = filename + ".tbi";
        StatGenStatus::Status libstatus = covtabix.readIndex(tabix_name.c_str());
        if (RegionStatus)
        {
            if (libstatus != StatGenStatus::SUCCESS)
            {
                error("Cannot open tabix file %s!\n", tabix_name.c_str());
            }
            bool status = SetIfilePosition(covfile_, covtabix, Chr, Start);
            if (!status)
            {
                error("Cannot find position %s:%d-%d in cov file %s!\n", Chr.c_str(), Start, End, filename.c_str());
            }
        }

        int m = 0;
        bool pass_header = 0;
        bool newFormat = 0;
        while (!ifeof(covfile_))
        {
            String buffer;
            buffer.ReadLine(covfile_);
            if (!pass_header)
            {
                if (buffer.Find("CHROM") == -1)
                {
                    continue;
                }
                // now check new or old format
                StringArray tokens;
                tokens.AddTokens(buffer, "\t ");
                if (tokens[2] == "MARKERS_IN_WINDOW" && tokens[3] == "COV_MATRICES")
                {
                    newFormat = 0;
                } else if (tokens[2] == "EXP" && tokens[3] == "COV_MATRICES")
                {
                    newFormat = 1;
                } else if (tokens[3] == "NUM_MARKER" && tokens[4] == "MARKER_POS" && tokens[5] == "COV")
                {
                    newFormat = 0;
                } else
                {
                    error("Covariance matrix is neither new or old format...are you using the right file?\n");
                }
                pass_header = 1;
                continue;
            }
            StringArray tokens;
            tokens.AddTokens(buffer, "\t ");
            if (RegionStatus)
            {
                if (tokens[1].AsInteger() > End || tokens[0] != Chr)
                { // out of this region or into another chromosome
                    break;
                }
            }
            if (FormatAdjust[study])
            { // rvtest
                readCovOldFormatLine(study, tokens, m);
            } else
            {
                if (newFormat)
                {
                    readCovNewFormatLine(study, tokens, m);
                } else
                {
                    readCovOldFormatLine(study, tokens, m);
                }
            }
        }
        ifclose(covfile_);
        printf("done\n");

        printf("Updating group stats ...\n");
        //update group statistics
        for (int g = 0; g < group.annoGroups.Length(); g++)
        {
            updateGroupStats(group, study, g, newFormat);
        }

        // clear after filling up each study
        markerPosHash.Clear();
        markersExp.Clear(); // for new format
        markersInWindow.Clear(); // for old format
        markersCov.Clear();
    }
}

void Meta::readCovOldFormatLine(int study, StringArray &tokens, int &m)
{
    removeChrFromString(tokens[0]);
    String SNP = tokens[0] + ":" + tokens[1];
    //exclude variants that have different REF and ALT
    if (simplifyCovLoad && groupAnchor.Integer(SNP) != -1)
    {
        return;
    }
    if (SNPexclude.Integer(String(std::to_string(study).c_str()) + ":" + SNP) != -1)
    {
        return;
    }
    m++;
    markerPosHash.SetInteger(SNP, m);
    markersInWindow.Push(tokens[marker_col]);
    markersCov.Push(tokens[cov_col]);
}


// do not have markersInWindow. Instead, read from markerList
void Meta::readCovNewFormatLine(int study, StringArray &tokens, int &m)
{
    removeChrFromString(tokens[0]);
    String SNP = tokens[0] + ":" + tokens[1];
    //exclude variants that have different REF and ALT
    if (simplifyCovLoad && groupAnchor.Integer(SNP) != -1)
    {
        return;
    }
    if (SNPexclude.Integer(String(std::to_string(study).c_str()) + ":" + SNP) != -1)
    {
        return;
    }
    m++;
    markerPosHash.SetInteger(SNP, m);
//	markerChrList.Push(tokens[0]);
//	markerPosList.Push(tokens[1]);
    markersExp.Push(tokens[2]);
    markersCov.Push(tokens[cov_col]);
}

// for new format
// read marker cov then add to vector
void Meta::addNewFormatCov(int mexp, String &cov_str, Vector &covs)
{
    StringArray commas;
    commas.AddTokens(cov_str, ',');
    // length of covs
    int n = commas.Length();
    if (n < 1)
    {
        error("At line: %s:...,no index of covariance matrices! Are you using the right cov file?\n", cov_str.c_str());
    }
    bool index_exist = 0;
    // now check if index exists
    StringArray first_tokens;
    first_tokens.AddTokens(commas[0], ":");
    if (first_tokens.Length() == 2)
    {
        index_exist = 1;
    } else if (first_tokens.Length() != 1)
    {
        error("At line: %s:...,abnormal index of covariance matrices! Are you using the right cov file?\n",
              cov_str.c_str());
    }

    // if no index, add directly
    if (index_exist)
    {
        int cov_len = -1;
        for (int i = n - 1; i >= 0; i--)
        {
            StringArray tokens;
            tokens.AddTokens(commas[i], ':');
            if (tokens.Length() == 1)
            {
                continue;
            }
            cov_len = tokens[1].AsInteger() + 1;
            break;
        }
        if (cov_len == -1)
        {
            error("At line: %s:...,no index of covariance matrices! Are you using the right cov file?\n",
                  commas[0].c_str());
        }
        covs.Dimension(cov_len, 0);
        // now add
        int last_index = 0;
        for (int i = 0; i < n; i++)
        {
            StringArray tokens;
            tokens.AddTokens(commas[i], ':');
            if (tokens.Length() > 2 || tokens.Length() < 1)
            {
                error("At line: ...:%s:...,abnormal separation!\n", commas[i].c_str());
            }
            if (tokens.Length() == 1)
            {
                last_index++;
            } else
            { // tokens ==2
                last_index = tokens[1].AsInteger();
            }
            //printf("covlen=%d,last_index=%d\n",cov_len,last_index);
            covs[last_index] = tokens[0].AsDouble() * pow(10, mexp);
        }
    } else
    {
        covs.Dimension(n, 0);
        for (int i = 0; i < n; i++)
        {
            covs[i] = commas[i].AsDouble() * pow(10, mexp);
        }
    }
}

void Meta::BurdenAssoc(String method, GroupFromAnnotation &group, Vector *&maf, Vector *&stats, Vector *&cond_stats,
                       Matrix *&cov, Matrix *&cond_cov, Vector *&singleEffSize, Vector *&singlePvalue)
{
    printf("Performing %s tests ...\n", method.c_str());
    //calculate final results here

    Vector pvalue_burden, pvalue_burden_cond;

    IFILE output;
    String filename;
    openMetaResultFile(prefix, filename, output, method);

    String method_out = method;
    method_out += "_";
    method_out += MAF_cutoff;

    IFILE reportOutput;
    String reportFile;
    if (report)
    {
        if (prefix == "")
        {
            reportFile = "meta.tophits." + method + ".tbl";
        } else if (prefix.Last() == '.' || prefix.Last() == '/')
        {
            reportFile = prefix + "meta.tophits." + method + ".tbl";
        } else
        {
            reportFile = prefix + ".meta.tophits." + method + ".tbl";
        }
        reportOutput = ifopen(reportFile, "w", InputFile::UNCOMPRESSED);
        if (cond != "")
        {
            ifprintf(reportOutput,
                     "GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\tCOND_EFFSIZES\tCOND_PVALUES\n");
        } else
        {
            ifprintf(reportOutput,
                     "GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\n");
        }
    }
    ifprintf(output, "##Method=Burden\n");
    ifprintf(output, "##STUDY_NUM=%d\n", scorefile.Length());
    ifprintf(output, "##TotalSampleSize=%d\n", total_N);
    if (cond != "")
    {
        if (fullResult)
        {
            ifprintf(output,
                     "#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
        } else
        {
            ifprintf(output,
                     "#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
        }
    } else if (fullResult)
    {
        ifprintf(output,
                 "#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
    } else
    {
        ifprintf(output, "#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
    }
    double numerator = _NAN_, denominator = _NAN_, chisq = 0.0, pvalue = _NAN_, cond_num = _NAN_, cond_denom = _NAN_;
    StringArray chr_plot, geneLabels;
    Vector pos_plot;
    for (int g = 0; g < group.annoGroups.Length(); g++)
    {
        if (maf[g].Length() == 0)
        {
            continue;
        }

        double average_af = maf[g].Average();
        double min_af = maf[g].Min();
        double max_af = maf[g].Max();

        String var;
        for (int i = 0; i < maf[g].Length() - 1; i++)
        {
            var += group.SNPlist[g][i] + ";";
        }
        var += group.SNPlist[g][maf[g].Length() - 1];

        Vector weight;
        weight.Dimension(maf[g].Length());
        SetWeight(method, weight, maf[g]);
        // for burden test, need to 1/w
        for (int w = 0; w < weight.Length(); w++)
        {
            weight[w] = 1 / weight[w];
        }

        numerator = weight.InnerProduct(stats[g]);
        Vector tmp;
        tmp.Dimension(group.SNPlist[g].Length());

        for (int i = 0; i < tmp.Length(); i++)
        {
            tmp[i] = weight.InnerProduct(cov[g][i]);
        }
        denominator = tmp.InnerProduct(weight);

        if (cond != "")
        {
            cond_num = weight.InnerProduct(cond_stats[g]);
            for (int i = 0; i < tmp.Length(); i++)
            {
                tmp[i] = weight.InnerProduct(cond_cov[g][i]);
            }
            cond_denom = tmp.InnerProduct(weight);
        }

        if (denominator == 0.0)
        {
            continue;
        }

        chisq = numerator * numerator / denominator;
        pvalue = pchisq(chisq, 1, 0, 0);
        double effSize = numerator / denominator;
        double cond_chisq = _NAN_, cond_effSize = _NAN_, cond_pvalue = _NAN_;

        bool disect = false;
        while (pvalue == 0.0)
        {
            disect = true;
            chisq *= 0.999;
            pvalue = pchisq(chisq, 1, 0, 0);
        }
        bool cond_disect = false;
        if (cond != "")
        {
            if (cond_denom == 0)
            {
                cond_effSize = 0.0;
                cond_pvalue = 1.0;
            } else
            {
                cond_chisq = cond_num * cond_num / cond_denom;
                cond_effSize = cond_num / cond_denom;
                cond_pvalue = pchisq(cond_chisq, 1, 0, 0);
                while (cond_pvalue == 0.0)
                {
                    cond_disect = true;
                    cond_chisq *= 0.999;
                    cond_pvalue = pchisq(cond_chisq, 1, 0, 0);
                }
            }
            pvalue_burden_cond.Push(cond_pvalue);
        }

        if (fullResult)
        {
            ifprintf(output, "%s\t%d\t%s\t", group.annoGroups[g].c_str(), group.SNPlist[g].Length(), var.c_str());

            for (int i = 0; i < maf[g].Length() - 1; i++)
            {
                ifprintf(output, "%g,", maf[g][i]);
            }
            ifprintf(output, "%g\t", maf[g][maf[g].Length() - 1]);

            for (int i = 0; i < singleEffSize[g].Length() - 1; i++)
            {
                ifprintf(output, "%g,", singleEffSize[g][i]);
            }
            ifprintf(output, "%g\t", singleEffSize[g][singleEffSize[g].Length() - 1]);

            for (int i = 0; i < singlePvalue[g].Length() - 1; i++)
            {
                ifprintf(output, "%g,", singlePvalue[g][i]);
            }
            ifprintf(output, "%g\t", singlePvalue[g][singlePvalue[g].Length() - 1]);

            if (cond != "")
            {
                ifprintf(output, "%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n", average_af, min_af, max_af, effSize,
                         disect ? "<" : "", pvalue, cond_effSize, cond_disect ? "<" : "", cond_pvalue);
            } else
            {
                ifprintf(output, "%g\t%g\t%g\t%g\t%s%g\n", average_af, min_af, max_af, effSize, disect ? "<" : "",
                         pvalue);
            }
        } else
        {
            if (cond != "")
            {
                ifprintf(output, "%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n", group.annoGroups[g].c_str(),
                         group.SNPlist[g].Length(), var.c_str(), average_af, min_af, max_af, effSize, disect ? "<" : "",
                         pvalue, cond_effSize, cond_disect ? "<" : "", cond_pvalue);
            } else
            {
                ifprintf(output, "%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\n", group.annoGroups[g].c_str(),
                         group.SNPlist[g].Length(), var.c_str(), average_af, min_af, max_af, effSize, disect ? "<" : "",
                         pvalue);
            }
        }
        if (pvalue < report_pvalue_cutoff && report)
        {
            StringArray variants;
            variants.AddTokens(var, ";");
            if (cond != "")
            {
                for (int v = 0; v < maf[g].Length(); v++)
                {
                    ifprintf(reportOutput, "%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n", group.annoGroups[g].c_str(),
                             method.c_str(), disect ? "<" : "", pvalue, cond_disect ? "<" : "", cond_pvalue, MAF_cutoff,
                             MAF_cutoff, variants[v].c_str(), maf[g][v], singleEffSize[g][v], singlePvalue[g][v]);
                }
            } else
            {
                for (int v = 0; v < maf[g].Length(); v++)
                {
                    ifprintf(reportOutput, "%s\t%s\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n", group.annoGroups[g].c_str(),
                             method.c_str(), disect ? "<" : "", pvalue, MAF_cutoff, MAF_cutoff, variants[v].c_str(),
                             maf[g][v], singleEffSize[g][v], singlePvalue[g][v]);
                }
            }
        }
        pvalue_burden.Push(pvalue);
        geneLabels.Push(group.annoGroups[g]);
        StringArray tmp_SNPname;
        tmp_SNPname.AddTokens(group.SNPlist[g][0], ":");
        chr_plot.Push(tmp_SNPname[0]);
        pos_plot.Push(tmp_SNPname[1].AsDouble());
    }
    String name = method;
    name += " (maf<";
    name += MAF_cutoff;
    name += ")";
    String extraname = "";
    String demo;
    if (pvalue_burden.Length() > 0)
    {
        //Calculate genomic control
        double GC = GetGenomicControlFromPvalue(pvalue_burden);
        demo = "GC = ";
        demo += GC;
        writepdf.Draw(pdf, geneLabels, pvalue_burden, chr_plot, pos_plot, name, extraname, demo, true);
        if (cond != "")
        {
            name += "conditional analysis";
            double GC = GetGenomicControlFromPvalue(pvalue_burden_cond);
            demo = "GC = ";
            demo += GC;
            writepdf.Draw(pdf, geneLabels, pvalue_burden_cond, chr_plot, pos_plot, name, extraname, demo, true);
        }
    }
    ifclose(output);
    if (report)
    {
        ifclose(reportOutput);
    }
    printf("  done.\n\n");
}

void Meta::VTassoc(GroupFromAnnotation &group)
{
    printf("Performing Variable Threshold tests ...\n");
    //calculate final results here
    Vector pvalue_VT, pos_plot, cond_pvalue_VT;
    StringArray chr_plot, geneLabels;

    IFILE output;
    String filename;
    String method = "VT_";
    openMetaResultFile(prefix, filename, output, method);

    method += MAF_cutoff;
    IFILE reportOutput;
    if (report)
    {
        String reportFile;
        if (prefix == "")
        {
            reportFile = "meta.tophits.VT.tbl";
        } else if (prefix.Last() == '.' || prefix.Last() == '/')
        {
            reportFile = prefix + "meta.tophits.VT.tbl";
        } else
        {
            reportFile = prefix + ".meta.tophits.VT.tbl";
        }
        reportOutput = ifopen(reportFile, "w", InputFile::UNCOMPRESSED);
        ifprintf(reportOutput, "GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\n");
    }

    ifprintf(output, "##Method=VT\n");
    ifprintf(output, "##STUDY_NUM=%d\n", scorefile.Length());
    ifprintf(output, "##TotalSampleSize=%d\n", total_N);
    if (fullResult)
    {
        ifprintf(output,
                 "#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\t");
    } else
    {
        ifprintf(output, "#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\t");
    }

    if (cond != "")
    {
        ifprintf(output, "EFFECT_SIZE\tMAF_CUTOFF\tCOND_PVALUE\n");
    } else
    {
        ifprintf(output, "\n");
    }

    for (int g = 0; g < group.annoGroups.Length(); g++)
    {
        if (g > 1 && g % 1000 == 1)
        {
            printf("Finished analyzing %d genes.\n", g - 1);
        }

        if (maf[g].Length() == 0)
        {
            continue;
        }

        if (maf[g].Length() == 1)
        {
            if (fullResult)
            {
                ifprintf(output, "%s\t1\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t", group.annoGroups[g].c_str(),
                         group.SNPlist[g][0].c_str(), maf[g][0], singleEffSize[g][0], singlePvalue[g][0], maf[g][0],
                         maf[g][0], maf[g][0], singleEffSize[g][0], maf[g][0], singlePvalue[g][0]);
            } else
            {
                ifprintf(output, "%s\t1\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t", group.annoGroups[g].c_str(),
                         group.SNPlist[g][0].c_str(), maf[g][0], maf[g][0], maf[g][0], singleEffSize[g][0], maf[g][0],
                         singlePvalue[g][0]);
            }
            pvalue_VT.Push(singlePvalue[g][0]);

            if (cond != "")
            {
                String SNPname_noallele;
                StringArray tmp;
                tmp.AddTokens(group.SNPlist[g][0], ":");
                SNPname_noallele = tmp[0] + ":" + tmp[1];
                double cond_pvalue_ = _NAN_, cond_U = _NAN_, cond_V = _NAN_, chisq = _NAN_, cond_effSize_ = _NAN_;
                bool disect = false;
                if (conditionVar.Integer(SNPname_noallele) == -1)
                {
                    cond_U = SNPstat_cond.Double(SNPname_noallele);
                    cond_V = SNP_Vstat_cond.Double(SNPname_noallele);
                    chisq = cond_U * cond_U / cond_V;
                    cond_pvalue_ = pchisq(chisq, 1, 0, 0);
                    cond_effSize_ = cond_U / cond_V;
                    while (cond_pvalue_ == 0.0)
                    {
                        disect = true;
                        chisq *= 0.999;
                        cond_pvalue_ = pchisq(chisq, 1, 0, 0);
                    }
                } else
                {
                    cond_effSize_ = 0.0;
                    cond_pvalue_ = 1.0;
                }
                cond_pvalue_VT.Push(cond_pvalue_);
                ifprintf(output, "%g\t%g\t%s%g", cond_effSize_, maf[g][0], disect ? "<" : "", cond_pvalue_);
            }
            ifprintf(output, "\n");
            StringArray tmp_SNPname;
            tmp_SNPname.AddTokens(group.SNPlist[g][0], ":");
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

        for (int i = 1; i < cp_maf.Length(); i++)
        {
            if (cp_maf[i] > maf_cutoff[maf_cutoff.Length() - 1])
            {
                maf_cutoff.Push(cp_maf[i]);
            }
        } //now unique maf cutoffs are saved in maf_cutoff.
        double pvalue = _NAN_, cond_pvalue = _NAN_;
        bool condition_status = false;
        pvalue = VTassocSingle(group, maf_cutoff, reportOutput, output, g, condition_status, method);
        if (cond != "")
        {
            condition_status = true;
            cond_pvalue = VTassocSingle(group, maf_cutoff, reportOutput, output, g, condition_status, method);
        }
        pvalue_VT.Push(pvalue);
        if (cond != "")
        {
            cond_pvalue_VT.Push(cond_pvalue);
        }

        StringArray tmp_SNPname;
        tmp_SNPname.AddTokens(group.SNPlist[g][0], ":");
        chr_plot.Push(tmp_SNPname[0]);
        pos_plot.Push(tmp_SNPname[1].AsDouble());
        geneLabels.Push(group.annoGroups[g]);
    }

    String name = "VT (maf<";
    name += MAF_cutoff;
    name += ")";
    String extraname = "";
    String demo = "";
    double GC = GetGenomicControlFromPvalue(pvalue_VT);
    demo = "GC=";
    demo += GC;
    writepdf.Draw(pdf, geneLabels, pvalue_VT, chr_plot, pos_plot, name, extraname, demo, true);
    if (cond != "")
    {
        name += " Conditional Analysis";
        double GC = GetGenomicControlFromPvalue(cond_pvalue_VT);
        demo = "GC=";
        demo += GC;
        writepdf.Draw(pdf, geneLabels, cond_pvalue_VT, chr_plot, pos_plot, name, extraname, demo, true);
    }

    ifclose(output);
    if (report)
    {
        ifclose(reportOutput);
    }
    printf("Done.\n\n");
}

double Meta::VTassocSingle(GroupFromAnnotation &group, Vector &maf_cutoff, IFILE reportOutput, IFILE output, int &g,
                           bool condition, String &method)
{
    double pvalue = _NAN_, chosen_cutoff = _NAN_, chosen_effSize = _NAN_;
    double numerator = 0.0, denominator = 0.0, t_max = _NAN_;
    Vector weight, tmp, chosen_weight, score;
    Matrix cov_weight;
    weight.Dimension(maf[g].Length());
    tmp.Dimension(group.SNPlist[g].Length());
    cov_weight.Dimension(maf_cutoff.Length(), maf[g].Length());
    score.Dimension(maf_cutoff.Length());
    for (int i = 0; i < maf_cutoff.Length(); i++)
    {
        for (int w = 0; w < weight.Length(); w++)
        {
            if (maf[g][w] <= maf_cutoff[i])
            {
                weight[w] = 1.0;
            } else
            {
                weight[w] = 0.0;
            }
            cov_weight[i][w] = weight[w];
        }
        if (condition)
        {
            numerator = weight.InnerProduct(cond_stats[g]);
        } else
        {
            numerator = weight.InnerProduct(stats[g]);
        }
        for (int d = 0; d < tmp.Length(); d++)
        {
            if (condition)
            {
                tmp[d] = weight.InnerProduct(cond_cov[g][d]);
            } else
            {
                tmp[d] = weight.InnerProduct(cov[g][d]);
            }
        }
        denominator = tmp.InnerProduct(weight);

        if (denominator != 0.0)
        {
            double t_stat = fabs(numerator / sqrt(denominator));
            score[i] = t_stat;
            if (t_max == _NAN_)
            {
                t_max = t_stat;
                chosen_cutoff = maf_cutoff[i];
                chosen_weight.Copy(weight);
                chosen_effSize = numerator / denominator;
            } else
            {
                if (t_stat > t_max)
                {
                    t_max = t_stat;
                    chosen_cutoff = maf_cutoff[i];
                    chosen_weight.Copy(weight);
                    chosen_effSize = numerator / denominator;
                }
            }
        } else
        {
            score[i] = 0.0;
        }
    }
    if (score.Max() == 0.0)
    {
        printf("Warning: group %s does not have qualified variants to group.\n", group.annoGroups[g].c_str());
        fprintf(log, "Warning: group %s does not have qualified variants to group.\n", group.annoGroups[g].c_str());
        return pvalue;
    }
    Vector tmp_maf, tmp_eff, tmp_pvalue;
    for (int i = 0; i < maf[g].Length(); i++)
    {
        if (chosen_weight[i] == 1.0)
        {
            tmp_maf.Push(maf[g][i]);
        }
    }
    for (int i = 0; i < maf[g].Length(); i++)
    {
        if (chosen_weight[i] == 1.0)
        {
            tmp_eff.Push(singleEffSize[g][i]);
            tmp_pvalue.Push(singlePvalue[g][i]);
        }
    }

    double average_af = tmp_maf.Average();
    double min_af = tmp_maf.Min();
    double max_af = tmp_maf.Max();
    String var;
    for (int i = 0; i < maf[g].Length() - 1; i++)
    {
        if (chosen_weight[i] == 1.0)
        {
            var += group.SNPlist[g][i] + ";";
        }
    }
    if (chosen_weight[maf[g].Length() - 1] == 1.0)
    {
        var += group.SNPlist[g][maf[g].Length() - 1];
    }

    //STEP3: calculate covariance matrix for (U_1 ... U_#cutoff)
    Matrix cov_U, cov_U_tmp;
    if (condition)
    {
        cov_U_tmp.Product(cov_weight, cond_cov[g]);
    } else
    {
        cov_U_tmp.Product(cov_weight, cov[g]);
    }
    Matrix cov_weight_trans;
    cov_weight_trans.Transpose(cov_weight);
    cov_U.Product(cov_U_tmp, cov_weight_trans); //now, cov(U) is saved in cov_U
    //Calculate covariance matrix for (T_1 ... T_#cutoff)
    Matrix cov_T;
    cov_T.Dimension(cov_U.rows, cov_U.cols);
    cov2cor(cov_U, cov_T);
    //STEP4: calculate VT pvalue and report.
    int cutoff = maf_cutoff.Length();
    double *lower = new double[cutoff];
    double *upper = new double[cutoff];
    double *mean = new double[cutoff];

    for (int i = 0; i < cutoff; i++)
    {
        mean[i] = 0.0;
        lower[i] = -t_max;
        upper[i] = t_max;
    }

    //Use pmvnorm to calculate the asymptotic p-value
    Vector result;
    pmvnorm(lower, upper, mean, cov_T, false, result);

    if (result[0] == -1.0)
    {
        if (!condition)
        {
            if (cond != "")
            {
                if (fullResult)
                {
                    ifprintf(output, "%s\t%d\t%s\t", group.annoGroups[g].c_str(), tmp_maf.Length(), var.c_str());

                    for (int i = 0; i < tmp_maf.Length() - 1; i++)
                    {
                        ifprintf(output, "%g,", tmp_maf[i]);
                    }
                    ifprintf(output, "%g\t", tmp_maf[tmp_maf.Length() - 1]);

                    for (int i = 0; i < tmp_eff.Length() - 1; i++)
                    {
                        ifprintf(output, "%g,", tmp_eff[i]);
                    }
                    ifprintf(output, "%g\t", tmp_eff[tmp_eff.Length() - 1]);

                    for (int i = 0; i < tmp_pvalue.Length() - 1; i++)
                    {
                        ifprintf(output, "%g,", tmp_pvalue[i]);
                    }
                    ifprintf(output, "%g\t", tmp_pvalue[tmp_pvalue.Length() - 1]);

                    ifprintf(output, "%g\t%g\t%g\t%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\t", average_af, min_af, max_af,
                             chosen_effSize, chosen_cutoff);
                } else
                {
                    ifprintf(output, "%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\t",
                             group.annoGroups[g].c_str(), tmp_maf.Length(), var.c_str(), average_af, min_af, max_af,
                             chosen_effSize, chosen_cutoff);
                }
            } else
            {
                ifprintf(output, "\n");
            }
        } else
        {
            ifprintf(output, "%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\n", chosen_effSize, chosen_cutoff);
        }
    } else
    {
        if (1.0 - result[0] == 0.0)
        {
            //           printf("gene %s has result %g\n",group.annoGroups[g].c_str(),1.0-result[0]);
            printf("Using Shuang's algorithm to calculate MVN pvalue for gene %s ... ", group.annoGroups[g].c_str());
            if (maf_cutoff.Length() > 20)
            {
                while (maf_cutoff.Length() > 20)
                {
                    maf_cutoff.Delete(0);
                }

                double numerator, denominator, t_max = _NAN_;
                Vector weight, tmp, chosen_weight, score;
                Matrix cov_weight;
                weight.Dimension(maf[g].Length());
                tmp.Dimension(group.SNPlist[g].Length());
                cov_weight.Dimension(maf_cutoff.Length(), maf[g].Length());
                for (int i = 0; i < maf_cutoff.Length(); i++)
                {
                    for (int w = 0; w < weight.Length(); w++)
                    {
                        if (maf[g][w] <= maf_cutoff[i])
                        {
                            weight[w] = 1.0;
                        } else
                        {
                            weight[w] = 0.0;
                        }
                        cov_weight[i][w] = weight[w];
                    }
                    if (condition)
                    {
                        numerator = weight.InnerProduct(cond_stats[g]);
                    } else
                    {
                        numerator = weight.InnerProduct(stats[g]);
                    }

                    for (int d = 0; d < tmp.Length(); d++)
                    {
                        if (condition)
                        {
                            tmp[d] = weight.InnerProduct(cond_cov[g][d]);
                        } else
                        {
                            tmp[d] = weight.InnerProduct(cov[g][d]);
                        }
                    }
                    denominator = tmp.InnerProduct(weight);
                    if (denominator != 0)
                    {
                        double t_stat = fabs(numerator / sqrt(denominator));
                        score.Push(t_stat);
                        if (t_max == _NAN_)
                        {
                            t_max = t_stat;
                            chosen_cutoff = maf_cutoff[i];
                            chosen_weight.Copy(weight);
                            chosen_effSize = numerator / denominator;
                        } else
                        {
                            if (t_stat > t_max)
                            {
                                t_max = t_stat;
                                chosen_cutoff = maf_cutoff[i];
                                chosen_weight.Copy(weight);
                                chosen_effSize = numerator / denominator;
                            }
                        }
                    } else
                    {
                        score.Push(0.0);
                    }
                }
                if (score.Max() == 0.0)
                {
                    printf("Warning: group %s does not have qualified variants to group.\n",
                           group.annoGroups[g].c_str());
                    fprintf(log, "Warning: group %s does not have qualified variants to group.\n",
                            group.annoGroups[g].c_str());
                    return pvalue;
                    printf("completed!\n");
                }
                Vector tmp_maf, tmp_eff, tmp_pvalue;
                for (int i = 0; i < maf[g].Length(); i++)
                {
                    if (chosen_weight[i] == 1.0)
                    {
                        tmp_maf.Push(maf[g][i]);
                    }
                }

                for (int i = 0; i < maf[g].Length(); i++)
                {
                    if (chosen_weight[i] == 1.0)
                    {
                        tmp_eff.Push(singleEffSize[g][i]);
                        tmp_pvalue.Push(singlePvalue[g][i]);
                    }
                }
                average_af = tmp_maf.Average();
                min_af = tmp_maf.Min();
                max_af = tmp_maf.Max();

                String var;
                for (int i = 0; i < maf[g].Length() - 1; i++)
                {
                    if (chosen_weight[i] == 1.0)
                    {
                        var += group.SNPlist[g][i] + ";";
                    }
                }
                if (chosen_weight[maf[g].Length() - 1] == 1.0)
                {
                    var += group.SNPlist[g][maf[g].Length() - 1];
                }
                //STEP3: calculate covariance matrix for (U_1 ... U_#cutoff)
                Matrix cov_U, cov_U_tmp;
                if (condition)
                {
                    cov_U_tmp.Product(cov_weight, cond_cov[g]);
                } else
                {
                    cov_U_tmp.Product(cov_weight, cov[g]);
                }
                Matrix cov_weight_trans;
                cov_weight_trans.Transpose(cov_weight);
                cov_U.Product(cov_U_tmp, cov_weight_trans); //now, cov(U) is saved in cov_U
                //Calculate covariance matrix for (T_1 ... T_#cutoff)
                Matrix cov_T;
                cov_T.Dimension(cov_U.rows, cov_U.cols);
                cov2cor(cov_U, cov_T);

                pvalue = CalculateMVTPvalue(score, cov_T, t_max);
                printf("completed!\n");
            } else
            {
                pvalue = CalculateMVTPvalue(score, cov_T, t_max);
                printf("completed!\n");
            }
        } else
        {
            pvalue = 1.0 - result[0];
        }

        if ((condition && cond != "") || cond == "")
        {
            if (pvalue < report_pvalue_cutoff && report)
            {
                StringArray variants;
                variants.AddTokens(var, ";");
                for (int v = 0; v < tmp_maf.Length(); v++)
                {
                    ifprintf(reportOutput, "%s\t%s\t%g\t%g\t%g\t%s\t%g\t%g\t%g\n", group.annoGroups[g].c_str(),
                             method.c_str(), pvalue, MAF_cutoff, chosen_cutoff, variants[v].c_str(), tmp_maf[v],
                             tmp_eff[v], tmp_pvalue[v]);
                }
            }
        }

        if (cond == "" || (!condition && cond != ""))
        {
            if (fullResult)
            {
                ifprintf(output, "%s\t%d\t%s\t", group.annoGroups[g].c_str(), tmp_maf.Length(), var.c_str());

                for (int i = 0; i < tmp_maf.Length() - 1; i++)
                {
                    ifprintf(output, "%g,", tmp_maf[i]);
                }
                ifprintf(output, "%g\t", tmp_maf[tmp_maf.Length() - 1]);

                for (int i = 0; i < tmp_eff.Length() - 1; i++)
                {
                    ifprintf(output, "%g,", tmp_eff[i]);
                }
                ifprintf(output, "%g\t", tmp_eff[tmp_eff.Length() - 1]);

                for (int i = 0; i < tmp_pvalue.Length() - 1; i++)
                {
                    ifprintf(output, "%g,", tmp_pvalue[i]);
                }
                ifprintf(output, "%g\t", tmp_pvalue[tmp_pvalue.Length() - 1]);

                ifprintf(output, "%g\t%g\t%g\t%g\t%g\t%g\t", average_af, min_af, max_af, chosen_effSize, chosen_cutoff,
                         pvalue);
            } else
            {
                ifprintf(output, "%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t", group.annoGroups[g].c_str(), tmp_maf.Length(),
                         var.c_str(), average_af, min_af, max_af, chosen_effSize, chosen_cutoff, pvalue);
            }
            if (cond == "")
            {
                ifprintf(output, "\n");
            }
        }

        if (cond != "" && condition)
        {
            ifprintf(output, "%g\t%g\t%g\n", chosen_effSize, chosen_cutoff, pvalue);
        }

        if (pvalue > 1.0)
        {
            pvalue = 1.0;
        }
    }
    if (lower)
    { delete[] lower; }
    if (upper)
    { delete[] upper; }
    if (mean)
    { delete[] mean; }
    return pvalue;
}

void Meta::SKATassoc(GroupFromAnnotation &group)
{
    printf("Performing SKAT ...\n");
    //calculate Q statistics here
    Vector pvalue_SKAT, pos_plot, cond_pvalue_SKAT;
    StringArray chr_plot, geneLabels;
    IFILE output;
    String filename;
    String method = "SKAT_";
    openMetaResultFile(prefix, filename, output, method);

    method += MAF_cutoff;
    IFILE reportOutput;
    if (report)
    {
        String reportFile;
        if (prefix == "")
        {
            reportFile = "meta.tophits.SKAT.tbl";
        } else if (prefix.Last() == '.' || prefix.Last() == '/')
        {
            reportFile = prefix + "meta.tophits.SKAT.tbl";
        } else
        {
            reportFile = prefix + ".meta.tophits.SKAT.tbl";
        }
        reportOutput = ifopen(reportFile, "w", InputFile::UNCOMPRESSED);
        ifprintf(reportOutput,
                 "GENE\tMETHOD\tGENE_PVALUE_DAVIES\tGENE_PVALUE_LIU\tMAF_CUTOFF\tACTUAL_CUTOFF\tVAR\tMAF\tEFFSIZE\tPVALUE\n");
    }

    ifprintf(output, "##Method=SKAT\n");
    ifprintf(output, "##STUDY_NUM=%d\n", scorefile.Length());
    ifprintf(output, "##TotalSampleSize=%d\n", total_N);
    if (fullResult)
    {
        ifprintf(output,
                 "#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tSTATISTICS\tPVALUE_DAVIES\tPVALUE_LIU\n");
    } else
    {
        ifprintf(output, "#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tSTATISTICS\tPVALUE_DAVIES\tPVALUE_LIU\n");
    }

    double Qstat, pvalue, pvalue_liu;
    for (int g = 0; g < group.annoGroups.Length(); g++)
    {
        //      printf("Now working on group %s\n",group.annoGroups[g].c_str());
        if (g % 1000 == 1 && g > 1000)
        {
            printf("Finished analyzing %d genes.\n", ((int) g / 1000) * 1000);
        }
        int n = maf[g].Length();
        Vector weight;
        weight.Dimension(n, 0);
        if (maf[g].Length() == 0)
        {
            continue;
        }

        double average_af = maf[g].Average();
        double min_af = maf[g].Min();
        double max_af = maf[g].Max();

        String var;
        for (int i = 0; i < maf[g].Length() - 1; i++)
        {
            var += group.SNPlist[g][i] + ";";
        }
        var += group.SNPlist[g][maf[g].Length() - 1];

        //get weight based on maf
//		double alpha = 1.0;
//		double beta=25.0;
        String skat_method = "BBeta";
        SetWeight(skat_method, weight, maf[g]);
        Vector tmp, cond_tmp;
        tmp.Dimension(n);
        if (cond != "")
        {
            cond_tmp.Dimension(n);
        }
        for (int i = 0; i < n; i++)
        {
            tmp[i] = weight[i] * stats[g][i];
            if (cond != "")
            {
                cond_tmp[i] = weight[i] * cond_stats[g][i];
            }
        }
        Qstat = tmp.InnerProduct(stats[g]);

        double cond_Qstat = _NAN_;
        if (cond != "")
        {
            cond_Qstat = cond_tmp.InnerProduct(cond_stats[g]);
        }
        double *lambda = new double[n];
        CalculateLambda(cov[g], weight, lambda);
        // check lambda for dead loop
        double lambda_sum = 0;
        for (int i = 0; i < n; i++)
        {
            lambda_sum += lambda[i];
        }
        if (lambda_sum < 0.0000000001)
        {
            fprintf(log, "Gene %s lambda sum is zero. Skipped!\n", group.annoGroups[g].c_str());
            continue;
        }

        double Qstat_dav = Qstat;
        double Qstat_liu = Qstat;
        double cond_pvalue = _NAN_, cond_pvalue_liu = _NAN_;

        pvalue = MixChidist(lambda, n, Qstat, "Davies");

        bool disect_davies = false, disect_liu = false, cond_disect_davies = false, cond_disect_liu = false;
        int disect_itr = 0;
        if (debug)
        {
            printf("Qstat=%g, pvalue=%g\n", Qstat, pvalue);
        }

        if (Qstat == _NAN_)
        {
            Qstat_dav = _NAN_;
            Qstat_liu = _NAN_;
            pvalue = _NAN_;
            pvalue_liu = _NAN_;
        } else
        {
            while ((pvalue <= 0.0 || pvalue == 2.0 || std::isnan(pvalue)) && disect_itr < 10000)
            {
                disect_davies = true;
                Qstat_dav *= 0.9999;
                pvalue = MixChidist(lambda, n, Qstat_dav, "Davies");
                disect_itr++;
            }
            while ((pvalue <= 0.0 || pvalue == 2.0 || std::isnan(pvalue)))
            {
                Qstat_dav *= 0.99;
                pvalue = MixChidist(lambda, n, Qstat_dav, "Davies");
            }
            pvalue_liu = MixChidist(lambda, n, Qstat_liu, "Liu");
            disect_itr = 0;
            while ((pvalue_liu <= 0.0 || pvalue_liu == 2.0 || std::isnan(pvalue_liu)) && disect_itr < 10000)
            {
                disect_liu = true;
                Qstat_liu *= 0.9999;
                pvalue_liu = MixChidist(lambda, n, Qstat_liu, "Liu");
                disect_itr++;
            }
            while ((pvalue_liu <= 0.0 || pvalue_liu == 2.0 || std::isnan(pvalue_liu)))
            {
                Qstat_liu *= 0.99;
                pvalue_liu = MixChidist(lambda, n, Qstat_liu, "Liu");
            }
        }

        if (cond != "")
        {
            double *lambda = new double[n];
            CalculateLambda(cond_cov[g], weight, lambda);
            Qstat_dav = cond_Qstat;
            Qstat_liu = cond_Qstat;

            cond_pvalue = MixChidist(lambda, n, cond_Qstat, "Davies");

            int disect_itr = 0;
            while ((cond_pvalue <= 0.0 || cond_pvalue == 2.0 || std::isnan(cond_pvalue)) && disect_itr < 10000)
            {
                cond_disect_davies = true;
                Qstat_dav *= 0.9999;
                pvalue = MixChidist(lambda, n, Qstat_dav, "Davies");
                disect_itr++;
            }
            while ((cond_pvalue <= 0.0 || cond_pvalue == 2.0 || std::isnan(cond_pvalue)))
            {
                Qstat_dav *= 0.99;
                cond_pvalue = MixChidist(lambda, n, Qstat_dav, "Davies");
            }
            cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu, "Liu");
            disect_itr = 0;
            while ((cond_pvalue_liu <= 0.0 || cond_pvalue_liu == 2.0 || std::isnan(cond_pvalue_liu)) &&
                   disect_itr < 10000)
            {
                cond_disect_liu = true;
                Qstat_liu *= 0.9999;
                cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu, "Liu");
                disect_itr++;
            }
            while ((cond_pvalue_liu <= 0.0 || cond_pvalue_liu == 2.0 || std::isnan(cond_pvalue_liu)))
            {
                Qstat_liu *= 0.99;
                cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu, "Liu");
            }
            if (lambda)
            { delete[] lambda; }
        }
        if (std::isnan(pvalue_liu) || std::isnan(pvalue))
        {
            if (fullResult)
            {
                ifprintf(output, "%s\t%d\t%s\t", group.annoGroups[g].c_str(), group.SNPlist[g].Length(), var.c_str());

                for (int i = 0; i < maf[g].Length() - 1; i++)
                {
                    ifprintf(output, "%g,", maf[g][i]);
                }
                ifprintf(output, "%g\t", maf[g][maf[g].Length() - 1]);

                for (int i = 0; i < singleEffSize[g].Length() - 1; i++)
                {
                    ifprintf(output, "%g,", singleEffSize[g][i]);
                }
                ifprintf(output, "%g\t", singleEffSize[g][singleEffSize[g].Length() - 1]);

                for (int i = 0; i < singlePvalue[g].Length() - 1; i++)
                {
                    ifprintf(output, "%g,", singlePvalue[g][i]);
                }
                ifprintf(output, "%g\t", singlePvalue[g][singlePvalue[g].Length() - 1]);
                if (cond == "")
                {
                    ifprintf(output, "%g\t%g\t%g\t%g\t-\t-\n", average_af, min_af, max_af, Qstat);
                } else
                {
                    ifprintf(output, "%g\t%g\t%g\t%g\t-\t-\t%g\t-\t-\n", average_af, min_af, max_af, Qstat, cond_Qstat);
                }
            } else
            {
                if (cond == "")
                {
                    ifprintf(output, "%s\t%d\t%s\t%g\t%g\t%g\t%g\t-\t-\n", group.annoGroups[g].c_str(),
                             group.SNPlist[g].Length(), var.c_str(), average_af, min_af, max_af, Qstat);
                } else
                {
                    ifprintf(output, "%s\t%d\t%s\t%g\t%g\t%g\t%g\t-\t-\t%g\t-\t-\n", group.annoGroups[g].c_str(),
                             group.SNPlist[g].Length(), var.c_str(), average_af, min_af, max_af, Qstat, cond_Qstat);
                }
            }
            continue;
        }
        //tabulate top results
        if (cond != "")
        {
            if ((cond_pvalue < report_pvalue_cutoff && cond_pvalue_liu < report_pvalue_cutoff) && report)
            {
                StringArray variants;
                variants.AddTokens(var, ";");
                for (int v = 0; v < maf[g].Length(); v++)
                {
                    ifprintf(reportOutput, "%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n", group.annoGroups[g].c_str(),
                             method.c_str(), cond_disect_davies ? "<" : "", cond_pvalue, cond_disect_liu ? "<" : "",
                             cond_pvalue_liu, MAF_cutoff, MAF_cutoff, variants[v].c_str(), maf[g][v],
                             singleEffSize[g][v], singlePvalue[g][v]);
                }
            } else if ((pvalue < report_pvalue_cutoff || pvalue_liu < report_pvalue_cutoff) && report)
            {
                StringArray variants;
                variants.AddTokens(var, ";");
                for (int v = 0; v < maf[g].Length(); v++)
                {
                    ifprintf(reportOutput, "%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n", group.annoGroups[g].c_str(),
                             method.c_str(), disect_davies ? "<" : "", pvalue, disect_liu ? "<" : "", pvalue_liu,
                             MAF_cutoff, MAF_cutoff, variants[v].c_str(), maf[g][v], singleEffSize[g][v],
                             singlePvalue[g][v]);
                }
            }
        }
        if (fullResult)
        {
            ifprintf(output, "%s\t%d\t%s\t", group.annoGroups[g].c_str(), group.SNPlist[g].Length(), var.c_str());

            for (int i = 0; i < maf[g].Length() - 1; i++)
            {
                ifprintf(output, "%g,", maf[g][i]);
            }
            ifprintf(output, "%g\t", maf[g][maf[g].Length() - 1]);

            for (int i = 0; i < singleEffSize[g].Length() - 1; i++)
            {
                ifprintf(output, "%g,", singleEffSize[g][i]);
            }
            ifprintf(output, "%g\t", singleEffSize[g][singleEffSize[g].Length() - 1]);

            for (int i = 0; i < singlePvalue[g].Length() - 1; i++)
            {
                ifprintf(output, "%g,", singlePvalue[g][i]);
            }
            ifprintf(output, "%g\t", singlePvalue[g][singlePvalue[g].Length() - 1]);
            if (cond == "")
            {
                ifprintf(output, "%g\t%g\t%g\t%g\t%s%g\t%s%g\n", average_af, min_af, max_af, Qstat,
                         disect_davies ? "<" : "", pvalue, disect_liu ? "<" : "", pvalue_liu);
            } else
            {
                ifprintf(output, "%g\t%g\t%g\t%g\t%s%g\t%s%g\t%g\t%s%g\t%s%g\n", average_af, min_af, max_af, Qstat,
                         disect_davies ? "<" : "", pvalue, disect_liu ? "<" : "", pvalue_liu, cond_Qstat,
                         cond_disect_davies ? "<" : "", cond_pvalue, cond_disect_liu ? "<" : "", cond_pvalue_liu);
            }
        } else
        {
            if (cond == "")
            {
                ifprintf(output, "%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%s%g\n", group.annoGroups[g].c_str(),
                         group.SNPlist[g].Length(), var.c_str(), average_af, min_af, max_af, Qstat,
                         disect_davies ? "<" : "", pvalue, disect_liu ? "<" : "", pvalue_liu);
            } else
            {
                ifprintf(output, "%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%s%g\t%g\t%s%g\t%s%g\n",
                         group.annoGroups[g].c_str(), group.SNPlist[g].Length(), var.c_str(), average_af, min_af,
                         max_af, Qstat, disect_davies ? "<" : "", pvalue, disect_liu ? "<" : "", pvalue_liu, cond_Qstat,
                         cond_disect_davies ? "<" : "", cond_pvalue, cond_disect_liu ? "<" : "", cond_pvalue_liu);
            }
        }

        if (lambda)
        { delete[] lambda; }

        pvalue_SKAT.Push(pvalue_liu);
        if (cond != "")
        {
            cond_pvalue_SKAT.Push(cond_pvalue_liu);
        }
        StringArray tmp_SNPname;
        tmp_SNPname.AddTokens(group.SNPlist[g][0], ":");
        chr_plot.Push(tmp_SNPname[0]);
        pos_plot.Push(tmp_SNPname[1].AsDouble());
        geneLabels.Push(group.annoGroups[g]);
    }

    String name = "SKAT (maf<";
    name += MAF_cutoff;
    name += ")";
    String extraname = "";
    String demo = "";
    double GC = GetGenomicControlFromPvalue(pvalue_SKAT);
    demo = "GC=";
    demo += GC;
    writepdf.Draw(pdf, geneLabels, pvalue_SKAT, chr_plot, pos_plot, name, extraname, demo, true);
    if (cond != "")
    {
        name += " conditional analysis";
        double GC = GetGenomicControlFromPvalue(cond_pvalue_SKAT);
        demo = "GC=";
        demo += GC;
        writepdf.Draw(pdf, geneLabels, cond_pvalue_SKAT, chr_plot, pos_plot, name, extraname, demo, true);
    }
    ifclose(output);
    if (report)
    {
        ifclose(reportOutput);
    }
    printf("Done.\n\n");
}

/*
void Meta::SKAToptimized( GroupFromAnnotation & group)
{
	double rho;
}
*/

void Meta::CalculateLambda(Matrix &cov, Vector &weight, double *lambda)
{
    int n = cov.rows;
    //calculat sqrt(V)
    Eigen::MatrixXd cov_eigen(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cov_eigen(i, j) = cov[i][j];
        }
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_cov(cov_eigen, Eigen::ComputeThinU);
    Eigen::MatrixXd final_eigen(n, n);
    Eigen::MatrixXd final_eigen_rhs(n, n);
    Eigen::MatrixXd final_eigen_lhs(n, n);
    Eigen::VectorXd tmp(n);
    for (int i = 0; i < n; i++)
    {
        tmp[i] = sqrt(svd_cov.singularValues()[i]);
    }
    Eigen::VectorXd weight_eigen(n);
    for (int i = 0; i < n; i++)
    {
        weight_eigen[i] = weight[i];
    }
    final_eigen_rhs = svd_cov.matrixU() * tmp.asDiagonal() * svd_cov.matrixU().transpose();
    final_eigen_lhs = final_eigen_rhs * weight_eigen.asDiagonal();
    final_eigen = final_eigen_lhs * final_eigen_rhs;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(final_eigen, Eigen::ComputeThinU);
    const Eigen::VectorXd S = svd.singularValues();
    for (int i = 0; i < n; i++)
    {
        lambda[i] = fabs(S[i]);
    }
}

// print error to both log and std err
void Meta::ErrorToLog(const char *msg)
{
    fprintf(log, "Error [Meta.cpp]: ");
    fwrite(msg, strlen(msg), 1, log);
    fprintf(log, "\n");
    error(msg);
}


/*** pop stratification ***/

// load 1KG MAF
void Meta::load1kgPopMAF()
{
    printf("\nloading 1KG population MAF...\n");
    printf("  this may take half an hour or so...\n");
    IFILE popfile = ifopen(popfile_name, "r");
    if (popfile == NULL)
    {
        error("Cannot open 1000G population MAF file: %s\n", popfile_name.c_str());
    }
    bool header_line = 1;
    while (!ifeof(popfile))
    {
        String buffer;
        buffer.ReadLine(popfile);
        StringArray tokens;
        tokens.AddTokens(buffer, '\t');
        if (header_line)
        {
            nPop = tokens.Length() - 7;
            header_line = 0;
            if (nPop <= 0)
            {
                error("No population in %s after column 8. Wrong format?\n", popfile_name.c_str());
            }
            continue;
        }
        StringArray alts;
        alts.AddTokens(tokens[4], ',');
        for (int a = 0; a < alts.Length(); a++)
        {
            String markername = tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" + alts[a];
            af1KG[markername].resize(nPop);
            for (int i = 0; i < nPop; i++)
            {
                StringArray maf_tokens;
                maf_tokens.AddTokens(tokens[7 + i], ',');
                if (maf_tokens.Length() != alts.Length())
                {
                    error("1KG maf file at line: \n %s \n alleles and mafs do not match!\n", buffer.c_str());
                }
                af1KG[markername][i] = maf_tokens[a].AsDouble();
            }
        }
    }
    ifclose(popfile);
    printf("\ndone.\n");
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
double Meta::getAFtilda(String &markername, double raw_af, int study)
{
//	double fk = raw_af - pgamma[study][0];
    double fk = raw_af;
    // match allele first & get fk
    std::map<String, std::vector<double> >::iterator ptr;
    ptr = af1KG.find(markername);
    // check MAF first
    if (raw_af < minMatchMAF || raw_af > maxMatchMAF)
    {
        fk = 0;
    } else
    {

        bool flip = false;
        bool match = true;;
        if (ptr == af1KG.end())
        { // check if flip
            StringArray tokens;
            tokens.AddTokens(markername, ":");
            String flip_markername = tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" + tokens[2];
            ptr = af1KG.find(flip_markername);
            if (ptr == af1KG.end())
            { // does not exist, use 1/4N or do not adjust
                if (matchOnly || matchDist > 0)
                { // no adjusted MAF
                    fk = 0;
                } else
                {
                    for (int i = 0; i < nPop; i++)
                    {
                        fk -= pgamma[study][i] * 0.25;
                    }
                }
                match = false;
            } else
            {
                flip = true;
            }
        }
        // check condition and calculate fk
        if (match)
        {
            double dist = 0;
            if (matchDist > 0)
            {
                bool correct = false;
                if (matchByAbs)
                {
                    for (int i = 0; i < nPop; i++)
                    {
                        double di = flip ? abs(raw_af - (1 - (ptr->second[i]))) : abs(raw_af - ptr->second[i]);
                        if (dist == 0)
                        {
                            dist = di;
                        } else
                        {
                            if (di < dist)
                            {
                                dist = di;
                            }
                        }
                    }
                } else
                {
                    for (int i = 0; i < nPop; i++)
                    {
                        dist += flip ? (raw_af - (1 - (ptr->second[i]))) * (raw_af - (1 - (ptr->second[i]))) :
                                (raw_af - ptr->second[i]) * (raw_af - ptr->second[i]);
                    }
                    dist = sqrt(dist / nPop);
                }
                // now do correction
                if (dist <= matchDist)
                {
                    for (int i = 0; i < nPop; i++)
                    {
                        fk -= flip ? pgamma[study][i] * (ptr->second[i]) : pgamma[study][i] * (1 - (ptr->second[i]));
                    }
                } else
                {
                    fk = 0;
                }
            } else
            {
                for (int i = 0; i < nPop; i++)
                {
                    fk -= flip ? pgamma[study][i] * (ptr->second[i]) : pgamma[study][i] * (1 - (ptr->second[i]));
                }
            }
        }

    }
    return fk;
}


void Meta::addToMapStrVec(std::map<String, std::vector<double> > &variant, int study, String &markername, double fk,
                          int vsize, double initial_number)
{
    std::map<String, std::vector<double> >::iterator p = variant.find(markername);
    if (p == variant.end())
    {
        variant[markername].resize(vsize, initial_number);
        variant[markername][study] = fk;
    } else
    {
        variant[markername][study] = fk;
    }
}

/**
 * Set weights for each variant. This calculates raw weights and does not specify how they are used in the calculation
 * 	(eg it will calculate w, not 1/w)
 * @param method
 * @param weight
 * @param maf
 */
void Meta::SetWeight(String &method, Vector &weight, Vector &maf)
{
    if (method == "burden")
    {
        // All weights should be equal (set =1.0)
        for (int w = 0; w < weight.Length(); w++)
        {
            weight[w] = 1.0;
        }
    } else if (method == "MB")
    {
        // Madsen-Browning: weight[w] = sqrt( maf[w] * (1. - maf[w]) );
        //      https://doi.org/10.1371/journal.pgen.1000384
        for (int w = 0; w < weight.Length(); w++)
        {
            weight[w] = sqrt(maf[w] * (1.0 - maf[w]));
        }
    } else if (method == "MAB")
    {
        // Determine weights based on MAF : w = MAF
        for (int w = 0; w < weight.Length(); w++)
        {
            weight[w] = maf[w];
        }
    } else if (method == "BBeta")
    {
        /*
         Truncated beta. Used by, eg, SKAT.
         The parameters used here are the default values recommended in the original SKAT paper, and de-emphasize
          common variants while promoting rare ones.
           "We suggest setting a1 = 1 and a2 (beta) = 25 because it increases the weight of rare variants while still
              putting decent nonzero weights for variants with MAF 1%–5%".
         See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135811/#sec2title
              and http://gero.usc.edu/CBPH/files/4_9_2013_PAA/15_Kardia_SequenceKernelAssociationTest.pdf

         If you wish to put larger emphasis on common variants, consider using MB or burden weight method instead.
        */
        double alpha = 1.0;
        double beta = 25.0;
        double f0 = 2 / Nsamples; // truncate at 4 alleles
        for (int w = 0; w < weight.Length(); w++)
        {
            double xmaf = maf[w];

            // Ensure that the allele frequency Xmaf is no smaller than 2/Nsamples or greater than 1-2/Nsamples
            if (xmaf > 0 && xmaf < f0)
            {
                xmaf = f0;
            }
            if (xmaf > 0 && xmaf > (1 - f0))
            {
                xmaf = 1 - f0;
            }
            double beta_density = GetBetaDensity(alpha, beta, xmaf);
            weight[w] = (beta_density * beta_density);
        }
    } else
    {
        error("Invalid weight %s!\n", method.c_str());
    }
}

// update v in cov matrix in exact
void Meta::updateExactCov(int study, int m, int s, StringArray &chr, StringArray &pos, Matrix &cov_i)
{
    String markername = chr[m] + ":" + pos[m];
    String markername2 = chr[s] + ":" + pos[s];
    std::map<String, std::vector<double> >::iterator pnk = variant_nk.find(markername);
    /*if (p == variant_nk.end()) {
		markername = markerName + ":" + alt_allele[s] + ":"+ref_allele[s];
		p = variant_nk.find(markername);
	}*/
    if (pnk == variant_nk.end())
    {
        error("error2137\n");
    }
    std::map<String, std::vector<double> >::iterator pfk = variant_fk.find(markername);
    std::map<String, std::vector<double> >::iterator pnk2;
    std::map<String, std::vector<double> >::iterator pfk2;
    if (m == s)
    {
        pnk2 = pnk;
        pfk2 = pfk;
    } else
    {
        pnk2 = variant_nk.find(markername2);
        pfk2 = variant_fk.find(markername2);
    }
    int nk = pnk->second[study];
    double nkfk2 = nk * pfk->second[study] * pfk2->second[study];
    double new_r = pnk->second[scorefile.Length()];
    double ff = pfk->second[scorefile.Length()] * pfk2->second[scorefile.Length()];
    // Vexact=Vrmw*sigma4
    cov_i[m][s] = new_r * (cov_i[m][s] * Ysigma2[study] - 4 * nk * ff + 4 * nkfk2);
}


/******* update covariance matrix for group test ****/
void Meta::updateGroupStats(GroupFromAnnotation &group, int study, int g, bool newFormat)
{
    //printf("doing group %d\n",g);
    int gvar_count = group.SNPlist[g].Length();
    StringArray chr, pos;
    for (int i = 0; i < gvar_count; i++)
    {
        StringArray tmp;
        tmp.AddTokens(group.SNPNoAllele[g][i], ":");
        chr.Push(tmp[0]);
        pos.Push(tmp[1]);
    }
    //now pos has all the positions of markers in group g.
    Matrix cov_i, GX;
    cov_i.Dimension(gvar_count, gvar_count, 0.0);
    if (cond != "")
    {
        GX.Dimension(gvar_count, XX_inv[study].cols, 0.0);
    }

    for (int m = 0; m < gvar_count; m++)
    {
        updateSingleVariantGroupStats(group, study, g, cov_i, GX, chr, pos, m, gvar_count, newFormat);
    }

    cov[g].Add(cov_i);
    if (cond != "")
    {
        cond_cov[g].Add(cov_i);
        Matrix GX_trans, tmp, extra_cov_i;
        GX_trans.Transpose(GX);
        tmp.Product(GX, XX_inv[study]);
        extra_cov_i.Product(tmp, GX_trans);
        extra_cov_i.Multiply(-1.0);
        cond_cov[g].Add(extra_cov_i);
    }
}

void Meta::updateSingleVariantGroupStats(GroupFromAnnotation &group, int study, int g, Matrix &cov_i, Matrix &GX,
                                         StringArray &chr, StringArray &pos, int m, int gvar_count, bool newFormat)
{
    int loc = markerPosHash.Integer(group.SNPNoAllele[g][m]);
    if (loc == -1)
    {
        return;
    }
    int skip = SNPexclude.Integer(String(std::to_string(study).c_str()) + ":" + group.SNPNoAllele[g][m]);
    if (skip != -1)
    {
        return;
    }
    int flip = flipSNP.Integer(String(std::to_string(study).c_str()) + ":" + group.SNPNoAllele[g][m]);
    double multiplyFactor = 1.0;
    if (flip != -1)
    {
        multiplyFactor = -1.0;
    }
    //read through markersInWindow and find the selected markers
    if (!newFormat)
    {
        updateSingleVariantGroupStatsOldFormat(group, study, g, cov_i, GX, chr, pos, loc, m, gvar_count,
                                               multiplyFactor);
    } else
    {
        updateSingleVariantGroupStatsNewFormat(group, study, g, cov_i, GX, chr, pos, loc, m, gvar_count,
                                               multiplyFactor);
    }
}

void
Meta::updateSingleVariantGroupStatsOldFormat(GroupFromAnnotation &group, int study, int g, Matrix &cov_i, Matrix &GX,
                                             StringArray &chr, StringArray &pos, int loc, int m, int gvar_count,
                                             double multiplyFactor)
{
    StringArray markers, markerscov;
    markers.AddTokens(markersInWindow[loc - 1], ",");
    markerscov.AddTokens(markersCov[loc - 1], ",");
    for (int s = m; s < gvar_count; s++)
    {
        int p = markers.Find(pos[s]);
        if (p == -1)
        {
            return;
        }
        String markerName = String(std::to_string(study).c_str()) + ":" + chr[s] + ":" + pos[s];
        //String markername = markerName + ":"+ref_allele[s] + ":"+alt_allele[s];
        //if the marker in window is supposed to be excluded due to non-consistent ref/alt allele
        //then skip
        int skip = SNPexclude.Integer(markerName);
        if (skip != -1)
        {
            return;
        }
        int flip = flipSNP.Integer(markerName);
        double factor = 1.0;
        if (flip != -1)
        {
            factor = -1.0;
        }
        cov_i[m][s] = multiplyFactor * factor * markerscov[p].AsDouble() * SampleSize[study];
        if (useExactMetaMethod)
        {
            updateExactCov(study, m, s, chr, pos, cov_i);
        }
    }
    // fill in GX
    if (cond != "")
    {
        StringArray name;
        name.AddTokens(group.SNPlist[g][m], ":");
        cond_stats[g][m] = SNPstat_cond.Double(name[0] + ":" + name[1]);
        String pos_str;
        for (int s = 0; s < commonVar_study[study].Length(); s++)
        {
            if (pos[m].AsInteger() < common_pos[commonVar_study[study][s]])
            { // direct fill in
                pos_str = common_pos[commonVar_study[study][s]];
                int p = markers.Find(pos_str);
                if (p == -1)
                {
                    continue;
                }
                GX[m][s] = markerscov[p].AsDouble() * SampleSize[study];
            } else
            { // need to put the smaller one first
                int loc = markerPosHash.Integer(
                        common_chr[commonVar_study[study][s]] + ":" + common_pos[commonVar_study[study][s]]);
                //If this SNP is not genotpyed in this study, then skip
                if (loc == -1)
                {
                    continue;
                }
                //read through markersInWindow and find the selected markers
                pos_str = pos[m];
                StringArray new_markers, new_markerscov;
                new_markers.AddTokens(markersInWindow[loc - 1], ",");
                new_markerscov.AddTokens(markersCov[loc - 1], ",");
                int p = new_markers.Find(pos_str);
                if (p == -1)
                {
                    continue;
                }
                GX[m][s] = new_markerscov[p].AsDouble() * SampleSize[study];
            }
        }
    }
}

void
Meta::updateSingleVariantGroupStatsNewFormat(GroupFromAnnotation &group, int study, int g, Matrix &cov_i, Matrix &GX,
                                             StringArray &chr, StringArray &pos, int loc, int m, int gvar_count,
                                             double multiplyFactor)
{
    Vector markerscov;
    addNewFormatCov(markersExp[loc - 1], markersCov[loc - 1], markerscov);
    for (int s = m; s < gvar_count; s++)
    {
        String mkname = chr[s] + ":" + pos[s];
        int p = markerPosHash.Integer(mkname);
        if (p == -1)
        {
            return;
        }
        p--;
        p -= m;
        String markerName = String(std::to_string(study).c_str()) + ":" + chr[s] + ":" + pos[s];
        int skip = SNPexclude.Integer(markerName);
        if (skip != -1)
        {
            return;
        }
        int flip = flipSNP.Integer(markerName);
        double factor = 1.0;
        if (flip != -1)
        {
            factor = -1.0;
        }
        if (p >= markerscov.Length())
        { // zeros in the tail
            cov_i[m][s] = 0;
        } else
        {
            cov_i[m][s] = multiplyFactor * factor * markerscov[p] * SampleSize[study];
        }
        if (useExactMetaMethod)
        {
            updateExactCov(study, m, s, chr, pos, cov_i);
        }
    }
    // fill in GX
    if (cond != "")
    {
        String pos_str;
        StringArray name;
        name.AddTokens(group.SNPlist[g][m], ":");
        cond_stats[g][m] = SNPstat_cond.Double(name[0] + ":" + name[1]);
        for (int s = 0; s < commonVar_study[study].Length(); s++)
        {
            if (pos[m].AsInteger() < common_pos[commonVar_study[study][s]])
            { // direct fill in
                pos_str = common_pos[commonVar_study[study][s]];
                String mkname = chr[s] + ":" + pos[s];
                int p = markerPosHash.Integer(mkname);
                if (p == -1)
                {
                    continue;
                }
                p -= m;
                GX[m][s] = markerscov[p] * SampleSize[study];
            } else
            { // need to put the smaller one first
                int loc = markerPosHash.Integer(
                        common_chr[commonVar_study[study][s]] + ":" + common_pos[commonVar_study[study][s]]);
                //If this SNP is not genotpyed in this study, then skip
                if (loc == -1)
                {
                    continue;
                }
                //read through markersInWindow and find the selected markers
                pos_str = pos[m];
                StringArray new_markers;
                Vector new_markerscov;
                String mkname = chr[s] + ":" + pos[s];
                int p = markerPosHash.Integer(mkname);
                p -= m;
                addNewFormatCov(markersExp[loc - 1], markersCov[loc - 1], new_markerscov);
                if (p == -1)
                {
                    continue;
                }
                GX[m][s] = new_markerscov[p] * SampleSize[study];
            }
        }
    }
}


