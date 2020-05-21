#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <stdio.h>
#include "StringBasics.h"
#include "../src/GroupFromAnnotation.h"
#include "../src/Meta.h"
#include "RMSingleVariantReader.h"
#include "RMGroupTestReader.h"
using namespace std;

TEST_CASE("Program Arguments") {
  SECTION("--range") {
    Meta meta;
    meta.prefix = "test.range";
    meta.setLogFile();
    meta.Region = "1:1-87";
    meta.SKAT = true;
    meta.scorefile.Add("tests/datasets/simulated/region/test.smallchunk.MetaScore.assoc.gz");
    meta.covfile.Add("tests/datasets/simulated/region/test.smallchunk.MetaCov.assoc.gz");

    GroupFromAnnotation group;
    group.groupFile = "tests/datasets/simulated/region/test.smallchunk.mask.tab";

    meta.Prepare();
    group.Run("", meta.log);
    meta.PoolSummaryStat(group);

    // This is pretty wild. If you don't run the printing routine, the single variant p-values aren't stored
    // internally, so group-based tests can't look them up. TODO: refactor this...
    meta.WriteSingleVariantResults(group);
    meta.Run(group);

    // Given the range above, only ZSYH2 should be tested and not the other gene.
    auto group_reader = RMGroupTestReader("test.range.meta.SKAT_.results");
    auto group_rec = group_reader.get_record("ZSYH2");
    auto num_group_rec = group_reader.get_num_records();

    REQUIRE(num_group_rec == 1);
    REQUIRE(group_rec->pvalue_liu == Approx(1.28628e-09));

    // Given the range above, the single variant results should only contain records from position 1:2 to 1:87.
    auto sv_reader = RMSingleVariantReader("test.range.meta.singlevar.results");
    auto num_sv_rec = sv_reader.get_num_records();
    auto sv_rec_first = *sv_reader.begin();
    auto sv_rec_last = *(--sv_reader.end());

    REQUIRE(num_sv_rec == 86);

    REQUIRE(sv_rec_first->chrom == "1");
    REQUIRE(sv_rec_first->pos == 2);
    REQUIRE(sv_rec_first->pvalue == Approx(0.1487579));

    REQUIRE(sv_rec_last->chrom == "1");
    REQUIRE(sv_rec_last->pos == 87);
    REQUIRE(sv_rec_last->pvalue == Approx(0.7183580));

    remove("test.range.raremetal.log");
    remove("test.range.meta.plots.pdf");
    remove("test.range.meta.singlevar.results");
    remove("test.range.meta.SKAT_.results");
    remove("raremetal.log");
  }
}

TEST_CASE("P-value precision") {
  SECTION("Score statistic resulting in very small p-value") {
    Meta meta;

    meta.setLogFile();

    meta.scorefile.Add("tests/raremetal/test_tut_rm/inputs/STUDY1.QT1.singlevar.score.txt.gz");
    meta.scorefile.Add("tests/raremetal/test_tut_rm/inputs/STUDY2.QT1.singlevar.score.txt.gz");

    meta.covfile.Add("tests/raremetal/test_tut_rm/inputs/STUDY1.QT1.singlevar.cov.txt.gz");
    meta.covfile.Add("tests/raremetal/test_tut_rm/inputs/STUDY2.QT1.singlevar.cov.txt.gz");

    GroupFromAnnotation group;
    meta.PoolSummaryStat(group);

    // This is actually an error in the tutorial dataset, where the alt allele is missing and the row is read incorrectly.
    // It ends up being convenient, though, because it creates a record with a very extremely small p-value.
    // The :432 is because raremetal ends up shifting the entire row left by 1 because of the missing alt entry.
    int idx = meta.SNPmaf_name.Find("9:494428375:G:432");

    String snp_name = meta.SNPmaf_name[idx];
    StringArray tmp;
    tmp.AddTokens(snp_name, ":");
    String SNPname_noallele = tmp[0] + ":" + tmp[1];
    int N = meta.SNP_effect_N[idx];
    double U = meta.SNPstat.Double(SNPname_noallele);
    double V = meta.SNP_Vstat.Double(SNPname_noallele);

    SingleVariantResult result(U, V, N);
    REQUIRE(result.effSize == Approx(1074.71));
    REQUIRE(result.log_pvalue == Approx(1727.694));

    // Catch2 can't test for approximate long doubles, have to log10 it
    double p_from_log = -static_cast<double>(log10(result.pvalue));
    REQUIRE(p_from_log == Approx(1727.694));
  }
}

TEST_CASE("Allele frequencies") {
  SECTION("Average and min/max") {
    Meta meta;
    meta.prefix = "test.allelefreq";
    meta.setLogFile();
    meta.averageFreq = true;
    meta.minMaxFreq = true;
    meta.scorefile.Add("tests/datasets/simulated/heterog/study0_raremetal.txt.gz");
    meta.scorefile.Add("tests/datasets/simulated/heterog/study1_raremetal.txt.gz");

    GroupFromAnnotation group;
    meta.Prepare();
    meta.PoolSummaryStat(group);
    meta.WriteSingleVariantResults(group);

    auto score_reader = RMSingleVariantReader("test.allelefreq.meta.singlevar.results");
    auto rec1 = score_reader.get_record("8:875238_G/C");

    REQUIRE(rec1->alt_af_mean == Approx(0.486337));
    REQUIRE(rec1->alt_af_se == Approx(0.0519997));
    REQUIRE(rec1->alt_af_min == Approx(0.4345));
    REQUIRE(rec1->alt_af_max == Approx(0.5385));

    // Note: when given the test files in order of study0, then study1, metal will select "A" as the
    // effect allele. However, raremetal will always report towards the alt allele, which is "G".
    // Remember this when interpreting allele frequencies/effects. A mean of 0.515486 is approx 1 - 0.48.
    auto rec2 = score_reader.get_record("3:1291852_A/G");
    REQUIRE(rec2->alt_af_mean == Approx(0.515486));
    REQUIRE(rec2->alt_af_se == Approx(0.0149234));
    REQUIRE(rec2->alt_af_min == Approx(0.502));
    REQUIRE(rec2->alt_af_max == Approx(0.532));

    remove("test.allelefreq.meta.singlevar.results");
    remove("test.allelefreq.meta.plots.pdf");
    remove("test.allelefreq.raremetal.log");
    remove("raremetal.log");
  }
}

TEST_CASE("File I/O") {
  SECTION("Simple meta-analysis") {
    Meta meta;
    meta.prefix = "test.fileio.simple";
    meta.setLogFile();

    meta.scorefile.Add("tests/raremetal/test_tut_rm/inputs/STUDY1.QT1.singlevar.score.txt.gz");
    meta.scorefile.Add("tests/raremetal/test_tut_rm/inputs/STUDY2.QT1.singlevar.score.txt.gz");

    meta.covfile.Add("tests/raremetal/test_tut_rm/inputs/STUDY1.QT1.singlevar.cov.txt.gz");
    meta.covfile.Add("tests/raremetal/test_tut_rm/inputs/STUDY2.QT1.singlevar.cov.txt.gz");

    GroupFromAnnotation group;
    meta.Prepare();
    meta.PoolSummaryStat(group);
    meta.WriteSingleVariantResults(group);

    auto score_reader = RMSingleVariantReader("test.fileio.simple.meta.singlevar.results");
    auto rec1 = score_reader.get_record("9:44001280_G/A");
    REQUIRE(rec1->pvalue == Approx(0.348151));
    REQUIRE(rec1->effect_size == Approx(0.423643));
    REQUIRE(rec1->effect_stderr == Approx(0.451557));
    REQUIRE(rec1->h2 == Approx(0.000860396));
    REQUIRE(rec1->pooled_alt_af == Approx(0.00244379));
    REQUIRE(score_reader.get_nstudies() == 2);

    auto rec2 = score_reader.get_record("9:494428375_G/432");
    double p_from_log = -static_cast<double>(log10(rec2->pvalue));
    REQUIRE(p_from_log == Approx(1727.694));

    remove("test.fileio.simple.meta.singlevar.results");
    remove("test.fileio.simple.meta.plots.pdf");
    remove("test.fileio.simple.raremetal.log");
    remove("raremetal.log");
  }
}

TEST_CASE("Heterogeneity statistics") {
  SECTION("Should be correct on simple example") {
    Meta meta;

    meta.setLogFile();
    meta.bHeterogeneity = true;
    meta.scorefile.Add("tests/datasets/simulated/heterog/study0_raremetal.txt.gz");
    meta.scorefile.Add("tests/datasets/simulated/heterog/study1_raremetal.txt.gz");

    GroupFromAnnotation group;
    meta.PoolSummaryStat(group);

    // High heterogeneity
    REQUIRE(meta.SNP_heterog_stat.Double("3:1291852") == Approx(109.990000));

    // Low heterogeneity
    REQUIRE(meta.SNP_heterog_stat.Double("8:875238") == Approx(0.357466));

    for (int i = 0; i < meta.SNPmaf_name.Length(); i++) {
      String snp_name = meta.SNPmaf_name[i];
      StringArray tmp;
      tmp.AddTokens(snp_name, ":");
      String snp_chrpos = tmp[0] + ":" + tmp[1];

      double het_chisq = meta.SNP_heterog_stat.Double(snp_chrpos);
      REQUIRE(het_chisq != _NAN_);
      if (het_chisq != _NAN_) { REQUIRE(het_chisq >= 0.0); }
    }

    // Check size match
    REQUIRE(meta.SNPstat.Entries() == meta.SNP_heterog_stat.Entries());

    remove("raremetal.log");
  }

  SECTION("Check heterogeneity limited to --range") {
    Meta meta;
    meta.prefix = "test.range.heterog";
    meta.setLogFile();
    meta.bHeterogeneity = true;
    meta.Region = "3:1291852-1955783";
    meta.scorefile.Add("tests/datasets/simulated/heterog/study0_raremetal.txt.gz");
    meta.scorefile.Add("tests/datasets/simulated/heterog/study1_raremetal.txt.gz");

    GroupFromAnnotation group;

    meta.Prepare();
    meta.PoolSummaryStat(group);
    meta.WriteSingleVariantResults(group);

    // Given the range above, the single variant results should only contain records from position 1:2 to 1:87.
    auto sv_reader = RMSingleVariantReader("test.range.heterog.meta.singlevar.results");
    auto num_sv_rec = sv_reader.get_num_records();
    auto sv_rec_first = *sv_reader.begin();
    auto sv_rec_last = *(--sv_reader.end());

    REQUIRE(num_sv_rec == 3);
    REQUIRE(sv_rec_first->chrom == "3");
    REQUIRE(sv_rec_first->pos == 1291852);
    REQUIRE(sv_rec_first->het_pvalue == Approx(9.84832e-26));
    REQUIRE(sv_rec_last->chrom == "3");
    REQUIRE(sv_rec_last->pos == 1955783);
    REQUIRE(sv_rec_last->het_pvalue == Approx(3.55586e-18));

    // All variants should have a heterogeneity p-value available
    for (int i = 0; i < meta.SNPmaf_name.Length(); i++) {
      String snp_name = meta.SNPmaf_name[i];
      StringArray tmp;
      tmp.AddTokens(snp_name, ":");
      String snp_chrpos = tmp[0] + ":" + tmp[1];

      double het_chisq = meta.SNP_heterog_stat.Double(snp_chrpos);
      REQUIRE(het_chisq != _NAN_);
      if (het_chisq != _NAN_) { REQUIRE(het_chisq >= 0.0); }
    }

    // Check size match
    REQUIRE(meta.SNPstat.Entries() == meta.SNP_heterog_stat.Entries());

    remove("test.range.heterog.raremetal.log");
    remove("test.range.heterog.meta.plots.pdf");
    remove("test.range.heterog.meta.singlevar.results");
    remove("test.range.heterog.meta.SKAT_.results");
    remove("raremetal.log");
  }
}

TEST_CASE("Tutorial datasets") {
  SECTION("tut_rm") {
    Meta meta;
    meta.prefix = "test.tut_rm";
    meta.SKAT = true;
    meta.Burden = true;
    meta.VTa = true;
    meta.setLogFile();

    meta.scorefile.Add("tests/raremetal/test_tut_rm/inputs/STUDY1.QT1.singlevar.score.txt.gz");
    meta.scorefile.Add("tests/raremetal/test_tut_rm/inputs/STUDY2.QT1.singlevar.score.txt.gz");

    meta.covfile.Add("tests/raremetal/test_tut_rm/inputs/STUDY1.QT1.singlevar.cov.txt.gz");
    meta.covfile.Add("tests/raremetal/test_tut_rm/inputs/STUDY2.QT1.singlevar.cov.txt.gz");

    GroupFromAnnotation group;
    group.groupFile = "tests/raremetal/test_tut_rm/inputs/group.file";

    meta.Prepare();
    group.Run("", meta.log);
    meta.PoolSummaryStat(group);
    meta.WriteSingleVariantResults(group);
    meta.Run(group);

    auto reader_sv_tested = RMSingleVariantReader("test.tut_rm.meta.singlevar.results");
    auto reader_sv_expect = RMSingleVariantReader("tests/raremetal/test_tut_rm/expected/COMBINED.QT1.meta.singlevar.results");
    REQUIRE(reader_sv_tested == reader_sv_expect);

    auto reader_skat_tested = RMGroupTestReader("test.tut_rm.meta.SKAT_.results");
    auto reader_skat_expect = RMGroupTestReader("tests/raremetal/test_tut_rm/expected/COMBINED.QT1.meta.SKAT_.results");
    REQUIRE(reader_skat_tested == reader_skat_expect);

    auto reader_burden_tested = RMGroupTestReader("test.tut_rm.meta.burden.results");
    auto reader_burden_expect = RMGroupTestReader("tests/raremetal/test_tut_rm/expected/COMBINED.QT1.meta.burden.results");
    REQUIRE(reader_burden_tested == reader_burden_expect);

    auto reader_vt_tested = RMGroupTestReader("test.tut_rm.meta.VT_.results");
    auto reader_vt_expect = RMGroupTestReader("tests/raremetal/test_tut_rm/expected/COMBINED.QT1.meta.VT_.results");
    REQUIRE(reader_vt_tested == reader_vt_expect);

    remove("test.tut_rm.meta.singlevar.results");
    remove("test.tut_rm.meta.SKAT_.results");
    remove("test.tut_rm.meta.burden.results");
    remove("test.tut_rm.meta.VT_.results");
    remove("test.tut_rm.meta.plots.pdf");
    remove("test.tut_rm.raremetal.log");
    remove("raremetal.log");
  }

  SECTION("tut_rm_rvt_cov") {
    // Run meta-analysis with scores/cov from RAREMETALWORKER.
    Meta meta_rm;
    meta_rm.prefix = "test.tut_rm_rvt_cov.raremetal";
    meta_rm.Burden = true;
    meta_rm.setLogFile();
    meta_rm.scorefile.Add("tests/raremetal/test_tut_rm_rvt_cov/inputs/STUDY1.QT1.singlevar.score.txt.gz");
    meta_rm.scorefile.Add("tests/raremetal/test_tut_rm_rvt_cov/inputs/STUDY2.QT1.singlevar.score.txt.gz");
    meta_rm.covfile.Add("tests/raremetal/test_tut_rm_rvt_cov/inputs/STUDY1.QT1.singlevar.cov.txt.gz");
    meta_rm.covfile.Add("tests/raremetal/test_tut_rm_rvt_cov/inputs/STUDY2.QT1.singlevar.cov.txt.gz");
    GroupFromAnnotation group_rm;
    group_rm.groupFile = "tests/raremetal/test_tut_rm_rvt_cov/inputs/group.file";
    meta_rm.Prepare();
    group_rm.Run("", meta_rm.log);
    meta_rm.PoolSummaryStat(group_rm);
    meta_rm.WriteSingleVariantResults(group_rm);
    meta_rm.Run(group_rm);

    // Run meta-analysis with scores/cov from RVTEST.
    Meta meta_rvtest;
    meta_rvtest.prefix = "test.tut_rm_rvt_cov.rvtest";
    meta_rvtest.Burden = true;
    meta_rvtest.setLogFile();
    meta_rvtest.scorefile.Add("tests/raremetal/test_tut_rm_rvt_cov/inputs/STUDY1.rvtests.MetaScore.assoc.gz");
    meta_rvtest.scorefile.Add("tests/raremetal/test_tut_rm_rvt_cov/inputs/STUDY2.rvtests.MetaScore.assoc.gz");
    meta_rvtest.covfile.Add("tests/raremetal/test_tut_rm_rvt_cov/inputs/STUDY1.rvtests.MetaCov.assoc.gz");
    meta_rvtest.covfile.Add("tests/raremetal/test_tut_rm_rvt_cov/inputs/STUDY2.rvtests.MetaCov.assoc.gz");
    GroupFromAnnotation group_rvtest;
    group_rvtest.groupFile = "tests/raremetal/test_tut_rm_rvt_cov/inputs/group.file";
    meta_rvtest.Prepare();
    group_rvtest.Run("", meta_rvtest.log);
    meta_rvtest.PoolSummaryStat(group_rvtest);
    meta_rvtest.WriteSingleVariantResults(group_rvtest);
    meta_rvtest.Run(group_rvtest);

    auto reader_raremetal = RMGroupTestReader("test.tut_rm_rvt_cov.raremetal.meta.burden.results");
    auto reader_rvtest = RMGroupTestReader("test.tut_rm_rvt_cov.rvtest.meta.burden.results");
    REQUIRE(reader_raremetal == reader_rvtest);

    remove("test.tut_rm_rvt_cov.raremetal.meta.singlevar.results");
    remove("test.tut_rm_rvt_cov.raremetal.meta.burden.results");
    remove("test.tut_rm_rvt_cov.raremetal.meta.plots.pdf");
    remove("test.tut_rm_rvt_cov.raremetal.raremetal.log");
    remove("test.tut_rm_rvt_cov.rvtest.meta.singlevar.results");
    remove("test.tut_rm_rvt_cov.rvtest.meta.burden.results");
    remove("test.tut_rm_rvt_cov.rvtest.meta.plots.pdf");
    remove("test.tut_rm_rvt_cov.rvtest.raremetal.log");
    remove("raremetal.log");
  }
}