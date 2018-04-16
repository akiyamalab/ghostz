/*
 * seed_searcher_query_parameters.h
 *
 *  Created on: 2013/04/23
 *      Author: shu
 */

#ifndef SEED_SEARCHER_QUERY_PARAMETERS_H_
#define SEED_SEARCHER_QUERY_PARAMETERS_H_

#include <stdint.h>
#include "seed_searcher_database_parameters.h"
#include "alphabet_coder.h"
#include "queries.h"

class SeedSearcherQueryParameters {
public:
	struct HashPositionData {
		bool similarity_check;
		uint32_t query_id;
		uint32_t position;
	};
	struct HashPositionDataList {
		SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash;
		std::vector<HashPositionData> position_data;
	};

	struct BuildParameters {
		AlphabetCoder::Code sequence_delimiter;
		uint32_t subsequence_length;
		std::vector<int> *ungapped_extension_cutoffs;
		std::vector<int> *gapped_extension_triggers;
		ScoreMatrix *score_matrix;
		SeedSearcherDatabaseParameters::SequencesIndex *sequences_index;
	};
	int Build(Queries &queries, BuildParameters &parameters);

	AlphabetCoder::Code GetSequenceDelimiter();
	std::vector<int> &GetUngappedExtensionCutoffs();
	std::vector<int> &GetGappedExtensionTriggers();
	ScoreMatrix &GetScoreMatrix();
	uint32_t GetHashPositionDataListSize();
	std::vector<AlphabetCoder::Code *> &GetQuerySequences();
	std::vector<HashPositionDataList> &GetHashPositionDataLists();

private:
	int BuildSubsequences(AlphabetCoder::Code *sequence, const uint32_t length,
			const AlphabetCoder::Code &sequence_delimiter,
			const uint32_t min_length);

	AlphabetCoder::Code sequence_delimiter_;
	std::vector<int> *ungapped_extension_cutoffs_ptr_;
	std::vector<int> *gapped_extension_triggers_ptr_;
	uint32_t number_using_hash_position_data_list_;
	ScoreMatrix score_matrix_;
	std::vector<AlphabetCoder::Code *> query_sequences_;
	std::vector<HashPositionDataList> hash_position_data_lists_;
};

#endif /* SEED_SEARCHER_QUERY_PARAMETERS_H_ */
