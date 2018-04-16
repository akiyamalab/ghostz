/*
 * dna_sequence.h
 *
 *  Created on: 2011/04/11
 *      Author: shu
 */

#ifndef DNA_SEQUENCE_H_
#define DNA_SEQUENCE_H_

#include "sequence.h"
#include <string.h>

using namespace std;

class DnaSequence : public Sequence {
public:
  static string GetComplementaryStrand(const string &dna);

  DnaSequence(const std::string &name, const std::string &sequence_data)
  : Sequence(name, sequence_data)
  {}

  virtual ~DnaSequence()
  {}

  string GetComplementaryStrand() const;

private:
  struct ToComplementaryLetter{
      int operator()(int c) {
        if (c == 'A' || c == 'a') {
          return 'T';
        } else if (c == 'C' || c == 'c') {
          return 'G';
        } else if (c == 'G' || c == 'g') {
          return 'C';
        } else if (c == 'T' || c == 't') {
          return 'A';
        } else {
          return 'N';
        }
      }
  };
};

#endif /* DNA_SEQUENCE_H_ */
