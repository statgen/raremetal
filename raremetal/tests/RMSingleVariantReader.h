#ifndef RAREMETAL_RMSINGLEVARIANTREADER_H
#define RAREMETAL_RMSINGLEVARIANTREADER_H

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
#include "catch.hpp"

struct RMSingleVariantRecord {
  std::string chrom;
  uint64_t pos;
  std::string ref;
  std::string alt;
  uint64_t n;
  double pooled_alt_af;
  std::string direction_by_study;
  double effect_size;
  double effect_stderr;
  double h2;
  double alt_af_mean = -1;
  double alt_af_se = -1;
  double alt_af_min = -1;
  double alt_af_max = -1;
  long double pvalue;
  long double het_pvalue = -1;

  bool operator==(const RMSingleVariantRecord &other) const;
  bool operator!=(const RMSingleVariantRecord &other) const;
};

class RMSingleVariantReader {
protected:
  uint64_t nsamples;
  uint64_t nstudies;
  std::string trait_name;
  std::vector<std::shared_ptr<RMSingleVariantRecord>> records;
  std::map<std::string, std::shared_ptr<RMSingleVariantRecord>> index;
  std::vector<std::string> header;
public:
  using record_iterator = std::vector<std::shared_ptr<RMSingleVariantRecord>>::const_iterator;
  using header_iterator = std::vector<std::string>::const_iterator;

  RMSingleVariantReader(const std::string &file);
  void load(const std::string &file);
  double get_nsamples();
  uint64_t get_nstudies();
  uint64_t get_num_records();
  std::shared_ptr<RMSingleVariantRecord> get_record(const std::string &i);

  record_iterator records_begin() const;
  record_iterator records_end() const;
  header_iterator header_begin() const;
  header_iterator header_end() const;

  bool operator==(const RMSingleVariantReader &other) const;
};

#endif //RAREMETAL_RMSINGLEVARIANTREADER_H
