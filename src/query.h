/*
 * query.h
 *
 *  Created on: 2011/02/21
 *      Author: shu
 */

#ifndef QUERY_H_
#define QUERY_H_

#include <string>
#include "alphabet_coder.h"
#include "statistics.h"
#include "sequence.h"

class Query {
public:

  Query();
  virtual ~Query();

  virtual std::string GetName() {
    return name_;
  }

  virtual AlphabetCoder::Code *GetSequence() {
    return &sequence_[0];
  }

  virtual uint32_t GetSequenceLength() {
    return sequence_.size();
  }

  virtual uint32_t GetRealSequenceLength() {
    return sequence_.size();
  }

  virtual uint32_t GetRealStart(uint32_t pos) {
    return pos;
  }

  virtual uint32_t GetRealEnd(uint32_t pos) {
    return pos;
  }

  virtual AlphabetCoder::Code GetSequenceDelimiter() {
    return sequence_delimiter_;
  }

  virtual Statistics::KarlinParameters GetUngappedKarlinParameters() {
    return ungapped_karlin_parameters_;
  }


protected:
  std::string name_;
  std::vector<AlphabetCoder::Code> sequence_;
  AlphabetCoder::Code sequence_delimiter_;
  Statistics::KarlinParameters ungapped_karlin_parameters_;
};

#endif /* QUERY_H_ */
