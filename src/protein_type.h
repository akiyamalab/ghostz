/*
 * protein.h
 *
 *  Created on: 2010/09/14
 *      Author: shu
 */

#ifndef PROTEIN_TYPE_H_
#define PROTEIN_TYPE_H_

#include "alphabet_type.h"
#include <string>

using namespace std;

class ProteinType : public AlphabetType {
public:
  string GetRegularLetters() const {
    return kRegularLetters;
  }

  string GetAmbiguousLetters() const {
    return kAmbiguousLetters;
  }

  char GetUnknownLetter() const {
    return kUnknownLetter;
  }

private:
  static const string kRegularLetters;
  static const string kAmbiguousLetters;
  static const char kUnknownLetter;
};

#endif /* PROTEIN_TYPE_H_ */
