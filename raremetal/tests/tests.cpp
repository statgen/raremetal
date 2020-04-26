#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <stdio.h>
#include "StringBasics.h"
#include "../src/GroupFromAnnotation.h"
#include "../src/Meta.h"

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
  }
}