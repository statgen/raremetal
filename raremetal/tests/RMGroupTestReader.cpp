#include "RMGroupTestReader.h"

using namespace std;

uint64_t RMGroupTestReader::get_num_records() {
  return records.size();
}

void RMGroupTestReader::load(const string &file) {
  if (!filepath_exists(file)) {
    throw std::runtime_error("Could not find file: " + file);
  }

  ifstream input_file(file);
  string line;
  auto line_separator = regex("[ \t]");

  // Line regexes
  auto regex_samples = regex("##TotalSampleSize=(\\d+)");
  auto regex_studies = regex("##STUDY_NUM=(\\d+)");
  auto regex_method = regex("##Method=(\\w+)");
  auto regex_header = regex("#GROUPNAME\tNUM_VAR.*");

  bool header_done = false;
  vector<string> header;
  while (getline(input_file, line)) {
    smatch match;
    if (!header_done) {
      if (regex_search(line, match, regex_samples) && match.size() > 1) {
        this->nsamples = stoul(match.str(1));
      }
      else if (regex_search(line, match, regex_studies) && match.size() > 1) {
        this->nstudies = stoul(match.str(1));
      }
      else if (regex_search(line, match, regex_method) && match.size() > 1) {
        this->group_test = match.str(1);
      }
      else if (regex_search(line, match, regex_header)) {
        header_done = true;
        copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(header));
      }
    }
    else {
      // Begin parsing record row
      vector<string> tokens;
      copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(tokens));

      // For some reason RAREMETAL puts a genomic control line all the way at the end of the file...
      if (line.substr(0,1) == "#") { continue; }

      // Create record
      auto rec = make_shared<RMGroupTestRecord>();
      for (int i = 0; i < tokens.size(); i++) {
        string col(header[i]);
        if (col == "#GROUPNAME") {
          rec->group = tokens.at(i);
        }
        else if (col == "NUM_VAR") {
          rec->num_var = stoul(tokens.at(i));
        }
        else if (col == "VARs") {
          auto semi = regex(";");
          copy(sregex_token_iterator(tokens.at(i).begin(), tokens.at(i).end(), semi, -1), sregex_token_iterator(), back_inserter(rec->variants));
        }
        else if (col == "PVALUE_LIU") {
          rec->pvalue_liu = extract_fp<double>(tokens.at(i));
        }
        else if (col == "PVALUE_DAVIES") {
          rec->pvalue_davies = extract_fp<double>(tokens.at(i));
        }
        else if (col == "PVALUE") {
          rec->pvalue = extract_fp<double>(tokens.at(i));
        }
        else if (col == "EFFECT_SIZE") {
          rec->effect_size = extract_fp<double>(tokens.at(i));
        }
      }

      // Insert
      records.push_back(rec);
      this->index.emplace(rec->group, rec);
    }
  }

  if (records.empty()) {
    throw std::runtime_error("No group test results were read from file: " + file);
  }
}

RMGroupTestReader::RMGroupTestReader(const string &file) {
  load(file);
}

double RMGroupTestReader::get_nsamples() {
  return nsamples;
}

uint64_t RMGroupTestReader::get_nstudies() {
  return nstudies;
}

shared_ptr<RMGroupTestRecord> RMGroupTestReader::get_record(const string &i) {
  auto it = this->index.find(i);
  if (it != this->index.end()) {
    return it->second;
  }
  else {
    auto null = shared_ptr<RMGroupTestRecord>();
    return null;
  }
}

bool RMGroupTestReader::operator==(const RMGroupTestReader &other) const {
  // Check basic stats.
  if (nsamples != other.nsamples) {
    return false;
  }

  if (nstudies != other.nstudies) {
    return false;
  }

  if (trait_name != other.trait_name) {
    return false;
  }

  if (records.size() != other.records.size()) {
    return false;
  }

  // Iterate over each record in this reader and compare to the other reader.
  for (uint64_t i = 0; i < records.size(); i++) {
    auto &rec1 = *(this->records[i]);
    auto &rec2 = *(other.records[i]);
    if (rec1 != rec2) {
      return false;
    }
  }

  return true;
}

bool RMGroupTestRecord::operator==(const RMGroupTestRecord &other) const {
  bool b_group = group == other.group;
  bool b_num_var = num_var == other.num_var;
  bool b_avg_af = approx_nan(avg_af, other.avg_af);
  bool b_min_af = approx_nan(min_af, other.min_af);
  bool b_max_af = approx_nan(max_af, other.max_af);
  bool b_stat = approx_nan(stat, other.stat);
  bool b_pval_davies = approx_nan(pvalue_davies, other.pvalue_davies);
  bool b_pval_liu = approx_nan(pvalue_liu, other.pvalue_liu);
  bool b_pval = approx_nan(pvalue, other.pvalue);
  bool b_effect = approx_nan(effect_size, other.effect_size);

  bool b_final = b_group && b_num_var && b_avg_af && b_min_af && b_max_af &&
                 b_stat && b_pval_davies && b_pval_liu && b_pval && b_effect;

  return b_final;
}

bool RMGroupTestRecord::operator!=(const RMGroupTestRecord &other) const {
  return !(*this == other);
}