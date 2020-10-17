#include "RMSingleVariantReader.h"

using namespace std;

uint64_t RMSingleVariantReader::get_num_records() {
  return records.size();
}

RMSingleVariantReader::record_iterator RMSingleVariantReader::records_begin() const {
  return records.begin();
}

RMSingleVariantReader::record_iterator RMSingleVariantReader::records_end() const {
  return records.end();
}

RMSingleVariantReader::header_iterator RMSingleVariantReader::header_begin() const {
  return header.begin();
}

RMSingleVariantReader::header_iterator RMSingleVariantReader::header_end() const {
  return header.end();
}

void RMSingleVariantReader::load(const string &file) {
  if (!filepath_exists(file)) {
    throw std::runtime_error("Could not find file: " + file);
  }

  ifstream input_file(file);
  string line;
  auto line_separator = regex("[ \t]");

  // Line regexes
  auto regex_samples = regex("##TotalSampleSize=(\\d+)");
  auto regex_studies = regex("##STUDY_NUM=(\\d+)");
  auto regex_header = regex("#CHROM\tPOS.*");

  bool header_done = false;
  header.clear();
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
      auto rec = make_shared<RMSingleVariantRecord>();
      for (int i = 0; i < tokens.size(); i++) {
        string col(header[i]);
        if (col == "#CHROM") {
          rec->chrom = tokens.at(i);
        }
        else if (col == "POS") {
          rec->pos = stoul(tokens.at(i));
        }
        else if (col == "REF") {
          rec->ref = tokens.at(i);
        }
        else if (col == "ALT") {
          rec->alt = tokens.at(i);
        }
        else if (col == "N_INFORMATIVE") {
          rec->n = stoul(tokens.at(i));
        }
        else if (col == "POOLED_ALT_AF") {
          rec->pooled_alt_af = extract_fp<double>(tokens.at(i));
        }
        else if (col == "DIRECTION_BY_STUDY") {
          rec->direction_by_study = tokens.at(i);
        }
        else if (col == "EFFECT_SIZE") {
          rec->effect_size = extract_fp<double>(tokens.at(i));
        }
        else if (col == "EFFECT_SIZE_SD") {
          rec->effect_stderr = extract_fp<double>(tokens.at(i));
        }
        else if (col == "H2") {
          rec->h2 = extract_fp<double>(tokens.at(i));
        }
        else if (col == "PVALUE") {
          rec->pvalue = extract_fp<long double>(tokens.at(i));
        }
        else if (col == "ALT_AF_MEAN") {
          rec->alt_af_mean = extract_fp<double>(tokens.at(i));
        }
        else if (col == "ALT_AF_SE") {
          rec->alt_af_se = extract_fp<double>(tokens.at(i));
        }
        else if (col == "ALT_AF_MIN") {
          rec->alt_af_min = extract_fp<double>(tokens.at(i));
        }
        else if (col == "ALT_AF_MAX") {
          rec->alt_af_max = extract_fp<double>(tokens.at(i));
        }
        else if (col == "HET_PVALUE") {
          rec->het_pvalue = extract_fp<long double>(tokens.at(i));
        }
      }

      // Keys
      string chrpos = rec->chrom + ":" + to_string(rec->pos);
      string variant = rec->chrom + ":" + to_string(rec->pos) + "_" + rec->ref + "/" + rec->alt;

      // Insert
      records.push_back(rec);
      this->index.emplace(chrpos, rec);
      this->index.emplace(variant, rec);
    }
  }

  if (records.empty()) {
    throw std::runtime_error("No single variant test results were read from file: " + file);
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

bool RMSingleVariantReader::operator==(const RMSingleVariantReader &other) const {
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

bool RMSingleVariantRecord::operator==(const RMSingleVariantRecord &other) const {
  bool b_chrom = chrom == other.chrom;
  bool b_pos = pos == other.pos;
  bool b_ref = ref == other.ref;
  bool b_alt = alt == other.alt;
  bool b_n = n == other.n;
  bool b_pooled_alt_af = approx_nan(pooled_alt_af, other.pooled_alt_af);
  bool b_direction = direction_by_study == other.direction_by_study;
  bool b_effect = approx_nan(effect_size, other.effect_size);
  bool b_se = approx_nan(effect_stderr, other.effect_stderr);
  bool b_h2 = approx_nan(h2, other.h2);
  bool b_pval = approx_nan(pvalue, other.pvalue);
  bool b_alt_af_mean = approx_nan(alt_af_mean, other.alt_af_mean);
  bool b_alt_af_se = approx_nan(alt_af_se, other.alt_af_se);
  bool b_alt_af_min = approx_nan(alt_af_min, other.alt_af_min);
  bool b_alt_af_max = approx_nan(alt_af_max, other.alt_af_max);
  bool b_het_pval = approx_nan(het_pvalue, other.het_pvalue);

  bool b_final = b_chrom && b_pos && b_ref && b_alt && b_n && b_pooled_alt_af && b_direction &&
    b_effect && b_se && b_h2 && b_pval && b_alt_af_mean && b_alt_af_se && b_alt_af_min &&
    b_alt_af_max && b_het_pval;

  return b_final;
}

bool RMSingleVariantRecord::operator!=(const RMSingleVariantRecord &other) const {
  return !(*this == other);
}