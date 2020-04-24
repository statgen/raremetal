#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <stdio.h>
#include "StringBasics.h"
#include "../src/GroupFromAnnotation.h"
#include "../src/Meta.h"

TEST_CASE("Heterogeneity statistics") {
  SECTION("Should be correct on simple example") {
    Meta meta;
    String filename;

    if (meta.prefix == "") {
      filename = "raremetal.log";
    }
    else if (meta.prefix.Last() == '.' || meta.prefix.Last() == '/') {
      filename = meta.prefix + "raremetal.log";
    }
    else {
      filename = meta.prefix + ".raremetal.log";
    }

    FILE *logFile = freopen(filename, "wt", stderr);
    meta.setLogFile(logFile);
    meta.skipOutput = true;
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