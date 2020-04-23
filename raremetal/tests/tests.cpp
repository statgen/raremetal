#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <stdio.h>
#include "StringBasics.h"
#include "../src/GroupFromAnnotation.h"
#include "../src/Meta.h"

TEST_CASE("Heterogeneity statistics") {
  SECTION("Should be correct on simple example") {
    FILE *logFile;
    String filename;
    if (Meta::prefix == "") {
      filename = "raremetal.log";
    }
    else if (Meta::prefix.Last() == '.' || Meta::prefix.Last() == '/') {
      filename = Meta::prefix + "raremetal.log";
    }
    else {
      filename = Meta::prefix + ".raremetal.log";
    }

    logFile = freopen(filename, "wt", stderr);

    String path;
    Meta meta(logFile);

    meta.scorefile.Add("tests/datasets/simulated/heterog/study0_raremetal.txt.gz");
    meta.scorefile.Add("tests/datasets/simulated/heterog/study1_raremetal.txt.gz");

    GroupFromAnnotation group;
    meta.Prepare();
    group.Run(path, logFile);
    meta.PoolSummaryStat(group);
    meta.Run(group);
  }
}