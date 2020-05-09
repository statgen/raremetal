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
  double avg_af;
  double min_af;
  double max_af;
  double stat;
  double pvalue_davies;
  double pvalue_liu;
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
};

#endif //RAREMETAL_RMGROUPTESTREADER_H
