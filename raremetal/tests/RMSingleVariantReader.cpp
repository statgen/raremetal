#include "RMSingleVariantReader.h"
#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <exception>
#include <limits>

using namespace std;

template <typename T> T extract_fp(const std::string &s) {
  T d;
  try {
    if (std::is_same<T, long double>::value) {
      d = stold(s);
    }
    else if (std::is_same<T, double>::value) {
      d = stod(s);
    }
    else if (std::is_same<T, float>::value) {
      d = stof(s);
    }
    else {
      throw std::invalid_argument("Invalid return type when extracting floating point type number from string");
    }
  }
  catch (const std::exception &e) {
    d = numeric_limits<T>::quiet_NaN();
  }
  return d;
}

void RMSingleVariantReader::load(const string &file) {
  ifstream input_file(file);
  string line;
  auto line_separator = regex("[ \t]");

  // Line regexes
  auto regex_samples = regex("##TotalSampleSize=(\\d+)");
  auto regex_studies = regex("##STUDY_NUM=(\\d+)");
  auto regex_header = regex("#CHROM\tPOS.*");

  bool header_done = false;
  while (getline(input_file, line)) {
    smatch match;
    if (!header_done) {
      if (regex_search(line, match, regex_samples) && match.size() > 1) {
        this->nsamples = stoul(match.str(1));
      }
      else if (regex_search(line, match, regex_studies) && match.size() > 1) {
        this->nstudies = stoul(match.str(1));
      }
      else if (regex_search(line, match, regex_header)) {
        header_done = true;
      }
    }
    else {
      // Begin parsing record row
      vector<string> tokens;
      copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(tokens));

      // For some reason RAREMETAL puts a genomic control line all the way at the end of the file...
      if (line.substr(0,1) == "#") { continue; }

      // Create record
      auto rec = make_shared<RMSingleVariantRecord>();
      rec->chrom = tokens.at(0);
      rec->pos = stoul(tokens.at(1));
      rec->ref = tokens.at(2);
      rec->alt = tokens.at(3);
      rec->n = stoul(tokens.at(4));
      rec->pooled_alt_af = stod(tokens.at(5));
      rec->direction_by_study = tokens.at(6);
      rec->effect_size = extract_fp<double>(tokens.at(7));
      rec->effect_stderr = extract_fp<double>(tokens.at(8));
      rec->h2 = extract_fp<double>(tokens.at(9));
      rec->pvalue = extract_fp<long double>(tokens.at(10));

      // Keys
      string chrpos = rec->chrom + ":" + to_string(rec->pos);
      string variant = rec->chrom + ":" + to_string(rec->pos) + "_" + rec->ref + "/" + rec->alt;

      // Insert
      records.push_back(rec);
      this->index.emplace(chrpos, rec);
      this->index.emplace(variant, rec);
    }
  }
}

RMSingleVariantReader::RMSingleVariantReader(const string &file) {
  load(file);
}

double RMSingleVariantReader::get_nsamples() {
  return nsamples;
}

uint64_t RMSingleVariantReader::get_nstudies() {
  return nstudies;
}

shared_ptr<RMSingleVariantRecord> RMSingleVariantReader::get_record(const string &i) {
  auto it = this->index.find(i);
  if (it != this->index.end()) {
    return it->second;
  }
  else {
    auto null = shared_ptr<RMSingleVariantRecord>();
    return null;
  }
}