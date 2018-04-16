/*
 * chain_filter.h
 *
 *  Created on: 2011/01/27
 *      Author: shu
 */

#ifndef CHAIN_FILTER_H_
#define CHAIN_FILTER_H_

#include <vector>
#include "seed_searcher.h"
#include "alphabet_coder.h"
#include "score_matrix.h"
using namespace std;

class ChainFilter {
public:
	typedef SeedSearcher::Hit Hit;
	ChainFilter(AlphabetCoder::Code sequence_delimiter,
			ScoreMatrix &score_matrix);
	void Filter(const AlphabetCoder::Code *query_sequence,
			uint32_t query_sequence_length, int cutoff,
			const AlphabetCoder::Code *database_sequence,
			uint32_t database_sequence_length, vector<Hit> &hits);

private:
	typedef struct {
		uint32_t query_end;
		uint32_t diagonal;
	} Chain;

	class DbStartComp {
	public:
		bool operator()(const Hit &a1, const Hit &a2) const {
			return a1.database_sequence_position < a2.database_sequence_position;
		}
	};

	bool Connected(const AlphabetCoder::Code *query_sequence,
			const AlphabetCoder::Code *database_sequence, uint32_t chain_query_end,
			Hit &hit, int cutoff);

	AlphabetCoder::Code sequence_delimiter_;
	ScoreMatrix score_matrix_;
	vector<Chain> chains_;
};

#endif /* CHAIN_FILTER_H_ */
