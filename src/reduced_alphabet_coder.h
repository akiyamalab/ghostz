/*
 * reduced_alphabet_coder.h
 *
 *  Created on: 2012/11/29
 *      Author: shu
 */

#ifndef REDUCED_ALPHABET_CODER_H_
#define REDUCED_ALPHABET_CODER_H_

#include <string>
#include <vector>
#include "alphabet_type.h"
#include "alphabet_coder.h"

class ReducedAlphabetCoder : public AlphabetCoder {
public:
  ReducedAlphabetCoder();
  ReducedAlphabetCoder(const AlphabetType &type);
  ReducedAlphabetCoder(const AlphabetType &type, const std::vector<std::string> &alphabet_sets);
  virtual ~ReducedAlphabetCoder();

  bool Set(const AlphabetType &type, const std::vector<std::string> &alphabet_sets);
private:
  char GetRepresentativeAlphabet(const std::vector<std::string> &alphabet_sets, const char c);
};

#endif /* REDUCED_ALPHABET_CODER_H_ */
