#ifndef RAREMETAL_READER_UTIL_H
#define RAREMETAL_READER_UTIL_H

#include <string>
#include <limits>
#include <exception>

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
    d = std::numeric_limits<T>::quiet_NaN();
  }
  return d;
}
#endif //RAREMETAL_READER_UTIL_H
