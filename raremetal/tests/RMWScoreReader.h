#ifndef RAREMETAL_RMWSCOREREADER_H
#define RAREMETAL_RMWSCOREREADER_H

#include <string>
#include <map>
#include <vector>
#include <memory>

struct RMWScoreRecord {
  std::string chrom;
  uint64_t pos;
  std::string ref;
  std::string alt;
  uint64_t n_informative;
  double founder_af;
  double all_af;
  uint64_t informative_alt_ac;
  double call_rate;
  double hwe_pvalue;
  uint64_t n_ref;
  uint64_t n_het;
  uint64_t n_alt;
  double u_stat;
  double sqrt_vstat;
  double alt_effsize;
  double pvalue;
};

class RMWScoreReader {
protected:
  uint64_t nsamples;
  double sigma2; // sigma_e2_hat from RAREMETAL
  std::string trait_name;
  std::vector<std::shared_ptr<RMWScoreRecord>> records;
  std::map<std::string, std::shared_ptr<RMWScoreRecord>> index;
public:
  RMWScoreReader(const std::string &file);
  void load(const std::string &file);
  double get_sigma2();
  double get_nsamples();
  std::shared_ptr<RMWScoreRecord> get_record(const std::string &i);
};

#endif //RAREMETAL_RMWSCOREREADER_H
