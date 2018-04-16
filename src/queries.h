/*
 * queries.h
 *
 *  Created on: 2009/10/09
 *      Author: shu
 */

#ifndef QUERIES_H_
#define QUERIES_H_

#include "sequence.h"
#include "fasta_sequence_reader.h"
#include "sequence_filter_interface.h"
#include "query.h"
#include "statistics.h"
#include "translator.h"
#include <vector>
#include <string>
#include <stdint.h>
#include <tr1/memory>

class Queries {
public:
  typedef struct {
	  bool filter;
    std::tr1::shared_ptr<AlphabetType> file_sequence_type_ptr;
    std::tr1::shared_ptr<AlphabetType> aligning_sequence_type_ptr;
    Statistics::KarlinParameters ungapped_karlin_parameters;
    unsigned int chunk_size;
    ScoreMatrix score_matrix;
    AlphabetCoder::Code sequence_delimiter;
  } Parameters;

  Queries(std::istream &in, const Parameters &parameters);

  virtual ~Queries()
  {}

  void Next();

  Query *GetQuery(uint32_t id) {
    return queries_[id].get();
  }
  uint32_t GetNumberSequences() {
    return queries_.size();
  }

  uint32_t GetMaxQuerySequenceLength() {
    return max_query_sequence_length_;
  }

private:
  typedef std::tr1::shared_ptr<Query> QueryPtr;
  bool SetQueries();
  QueryPtr BuildQuery(Sequence &sequence);

  Parameters parameters_;
  std::vector<QueryPtr> queries_;
  uint32_t max_query_sequence_length_;
  Translator translator_;
  std::tr1::shared_ptr<SequenceFilterInterface> sequence_filter_;
  FastaSequenceReader reader_;
  QueryPtr next_query_;
};

#endif /* QUERIES_H_ */
