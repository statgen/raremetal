#ifndef RAREMETAL_RMSINGLEVARIANTREADER_H
#define RAREMETAL_RMSINGLEVARIANTREADER_H

#include <string>
#include <map>
#include <vector>
#include <memory>

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
  double pvalue;
};

class RMSingleVariantReader {
protected:
  uint64_t nsamples;
  uint64_t nstudies;
  std::string trait_name;
  std::vector<std::shared_ptr<RMSingleVariantRecord>> records;
  std::map<std::string, std::shared_ptr<RMSingleVariantRecord>> index;
public:
  RMSingleVariantReader(const std::string &file);
  void load(const std::string &file);
  double get_nsamples();
  uint64_t get_nstudies();
  std::shared_ptr<RMSingleVariantRecord> get_record(const std::string &i);
};

#endif //RAREMETAL_RMSINGLEVARIANTREADER_H
