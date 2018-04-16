/*
 * dna.h
 *
 *  Created on: 2010/09/28
 *      Author: shu
 */

#ifndef DNA_TYPE_H_
#define DNA_TYPE_H_

#include "alphabet_type.h"
#include <string>

using namespace std;

class DnaType : public AlphabetType {
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

#endif /* DNA_TYPE_H_ */
