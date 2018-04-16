/*
 * query.cpp
 *
 *  Created on: 2009/10/09
 *      Author: shu
 */

#include "fasta_sequence_reader.h"
#include "sequence.h"
#include "alphabet_coder.h"
#include "dna_sequence.h"
#include "dna_type.h"
//#include "dna_query.h"
#include "sequence_no_filter.h"
#include "sequence_seg_filter.h"
#include "protein_type.h"
#include "statistics.h"
#include "translator.h"
#include "translated_dna_query.h"
#include "protein_query.h"
#include "protein_type.h"
#include "queries.h"
#include "logger.h"
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <numeric>
#include <typeinfo>
#include <string>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <functional>

using namespace std;

Queries::Queries(istream &in, const Parameters &parameters) :
		parameters_(parameters), max_query_sequence_length_(0), reader_(in), next_query_(
				QueryPtr()) {
	if (parameters_.filter) {
		if (typeid(*(parameters_.aligning_sequence_type_ptr))
				== typeid(ProteinType)) {
			sequence_filter_ = std::tr1::shared_ptr<SequenceFilterInterface>(
					new SequenceSegFilter());
		} else {
			Logger *logger = Logger::GetInstance();
			logger->ErrorLog("can't use a filter for dna sequences");
			exit(1);
		}
	} else {
		sequence_filter_ = std::tr1::shared_ptr<SequenceFilterInterface>(
				new SequenceNoFilter());
	}
	Next();
}

void Queries::Next() {
	max_query_sequence_length_ = 0;
	SetQueries();
}

bool Queries::SetQueries() {
	queries_.clear();
	max_query_sequence_length_ = 0;
	unsigned int chunk_size = 0;
	string name;
	string sequence;
	bool reader_ret;
	if (!next_query_) {
		reader_ret = reader_.Read(name, sequence);
		if (reader_ret) {
			Sequence s(name, sequence);
			next_query_ = BuildQuery(s);
		}
	}

	while (1) {
		if (!next_query_) {
			break;
		}
		max_query_sequence_length_ = max(max_query_sequence_length_,
				next_query_->GetSequenceLength());
		unsigned int new_chunk_size = chunk_size
				+ next_query_->GetSequenceLength();
		if (new_chunk_size >= parameters_.chunk_size) {
			break;
		}
		queries_.push_back(next_query_);
		chunk_size = new_chunk_size;
		reader_ret = reader_.Read(name, sequence);
		if (reader_ret) {
			Sequence s(name, sequence);
			next_query_ = BuildQuery(s);
		} else {
			next_query_ = QueryPtr();
			break;
		}
	}
	return chunk_size;
}

Queries::QueryPtr Queries::BuildQuery(Sequence &sequence) {
	if (typeid(*(parameters_.file_sequence_type_ptr)) == typeid(DnaType)) {
		if (typeid(*(parameters_.aligning_sequence_type_ptr))
				!= typeid(ProteinType)) {
			Logger *logger = Logger::GetInstance();
			logger->ErrorLog("can't use aligner for query dna");
			exit(1);
		} else {
			return Queries::QueryPtr(
					new TranslatedDnaQuery(sequence,
							parameters_.sequence_delimiter, translator_,
							sequence_filter_.get(), parameters_.score_matrix,
							parameters_.ungapped_karlin_parameters));
		}
	} else {
		return Queries::QueryPtr(
				new ProteinQuery(sequence, parameters_.sequence_delimiter,
						sequence_filter_.get(), parameters_.score_matrix,
						parameters_.ungapped_karlin_parameters));
	}
	return QueryPtr();
}

