#include "RMGroupTestReader.h"

using namespace std;

uint64_t RMGroupTestReader::get_num_records() {
  return records.size();
}

void RMGroupTestReader::load(const string &file) {
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
      }

      // Insert
      records.push_back(rec);
      this->index.emplace(rec->group, rec);
    }
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