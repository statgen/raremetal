#ifndef RAREMETAL_RMGROUPTESTREADER_H
#define RAREMETAL_RMGROUPTESTREADER_H

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <exception>
#include <limits>
#include "reader_util.h"

struct RMGroupTestRecord {
  std::string group;
  uint64_t num_var;
  std::vector<std::string> variants;
  double avg_af = -1;
  double min_af = -1;
  double max_af = -1;
  double stat = -1;
  double pvalue_davies = -1;
  double pvalue_liu = -1;
  double pvalue = -1;
  double effect_size = -1;

  bool operator==(const RMGroupTestRecord &other) const;
  bool operator!=(const RMGroupTestRecord &other) const;
};

class RMGroupTestReader {
protected:
  uint64_t nsamples;
  uint64_t nstudies;
  std::string group_test;
  std::string trait_name;
  std::vector<std::shared_ptr<RMGroupTestRecord>> records;
  std::map<std::string, std::shared_ptr<RMGroupTestRecord>> index;
public:
  RMGroupTestReader(const std::string &file);
  void load(const std::string &file);
  double get_nsamples();
  uint64_t get_nstudies();
  std::string get_group_test();
  uint64_t get_num_records();
  std::shared_ptr<RMGroupTestRecord> get_record(const std::string &i);

  bool operator==(const RMGroupTestReader &other) const;
};

#endif //RAREMETAL_RMGROUPTESTREADER_H
