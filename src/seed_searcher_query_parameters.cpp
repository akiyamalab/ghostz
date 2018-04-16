/*
 * seed_searcher_query_parameters.cpp
 *
 *  Created on: 2013/04/23
 *      Author: shu
 */

#include "seed_searcher_common.h"
#include "seed_searcher_query_parameters.h"

int SeedSearcherQueryParameters::Build(Queries &queries,
		BuildParameters &parameters) {
	sequence_delimiter_ = parameters.sequence_delimiter;
	number_using_hash_position_data_list_ = 0;
	score_matrix_ = *parameters.score_matrix;
	ungapped_extension_cutoffs_ptr_ = parameters.ungapped_extension_cutoffs;
	gapped_extension_triggers_ptr_ = parameters.gapped_extension_triggers;
	query_sequences_.clear();
	uint32_t subsequence_length = parameters.subsequence_length;
	typedef std::tr1::unordered_map<
			SeedSearcherDatabaseParameters::SequencesIndex::HashKey, uint32_t> HashMap;
	HashMap hash_map;
	for (uint32_t query_id = 0; query_id < queries.GetNumberSequences();
			++query_id) {
		Query *query = queries.GetQuery(query_id);
		AlphabetCoder::Code *query_sequence = query->GetSequence();
		uint32_t query_sequence_length = query->GetSequenceLength();
		query_sequences_.push_back(query_sequence);
		uint32_t subsequence_offset = 0;
		while (subsequence_offset < query_sequence_length) {
			uint32_t subsequence_end = subsequence_offset;
			for (;
					(subsequence_end < query_sequence_length)
							&& (query_sequence[subsequence_end]
									!= sequence_delimiter_);
					++subsequence_end) {
			}
			if ((subsequence_end - subsequence_offset) >= subsequence_length) {
				uint32_t similarity_check_start = subsequence_offset
						+ subsequence_length / 2;
				uint32_t similarity_check_end = subsequence_end
						- subsequence_length / 2;
				for (uint32_t query_sequence_i = subsequence_offset;
						query_sequence_i < subsequence_end;
						++query_sequence_i) {
					SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash =
							0;
					uint32_t hashed_sequence_length = 0;
					if (!parameters.sequences_index->GetHashKey(
							query_sequence + query_sequence_i, &hash,
							&hashed_sequence_length)) {
						uint32_t seed_position =
								SeedSearcherCommon::CalculateSeedPosition(
										query_sequence_i,
										hashed_sequence_length);
						HashMap::iterator find_it = hash_map.find(hash);
						HashPositionData hash_position_data;
						if (similarity_check_start <= seed_position
								&& query_sequence_i < similarity_check_end) {
							hash_position_data.similarity_check = true;
						} else {
							hash_position_data.similarity_check = false;
						}
						hash_position_data.query_id = query_id;
						hash_position_data.position = seed_position;
						if (find_it == hash_map.end()) {
							hash_map.insert(
									std::make_pair(hash,
											number_using_hash_position_data_list_));
							++number_using_hash_position_data_list_;
							HashPositionDataList hash_position_data_list;
							hash_position_data_list.hash = 0;
							if (hash_position_data_lists_.size()
									<= number_using_hash_position_data_list_) {
								hash_position_data_lists_.push_back(
										hash_position_data_list);
							}
							hash_position_data_lists_[number_using_hash_position_data_list_
									- 1].hash = hash;
							hash_position_data_lists_[number_using_hash_position_data_list_
									- 1].position_data.clear();
							hash_position_data_lists_[number_using_hash_position_data_list_
									- 1].position_data.push_back(
									hash_position_data);

						} else {
							hash_position_data_lists_[find_it->second].position_data.push_back(
									hash_position_data);
						}
					}
				}
			}
			subsequence_offset = subsequence_end + 1;
		}
	}
#if 0
	for (uint32_t i = 0; i < number_using_hash_position_data_list_; ++i) {
		HashPositionDataList &hash_position_data_list =
		hash_position_data_lists_[i];
		for (uint32_t j = 0; j < hash_position_data_list.position_data.size();
				++j) {
			Query *query = queries.GetQuery(
					hash_position_data_list.position_data[j].query_id);
			AlphabetCoder::Code *query_sequence = query->GetSequence();
			uint32_t query_sequence_length = query->GetSequenceLength();
			assert(
					hash_position_data_list.position_data[j].position
					< query_sequence_length);
			SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash = 0;
			if (!parameters.sequences_index->GetHashKey(
							query_sequence
							+ hash_position_data_list.position_data[j].position,
							&hash)) {
				assert(hash == hash_position_data_list.hash);
			}
		}
	}
#endif

	return 0;
}

AlphabetCoder::Code SeedSearcherQueryParameters::GetSequenceDelimiter() {
	return sequence_delimiter_;
}
std::vector<int> &SeedSearcherQueryParameters::GetUngappedExtensionCutoffs() {
	return *ungapped_extension_cutoffs_ptr_;
}
std::vector<int> &SeedSearcherQueryParameters::GetGappedExtensionTriggers() {
	return *gapped_extension_triggers_ptr_;
}
ScoreMatrix &SeedSearcherQueryParameters::GetScoreMatrix() {
	return score_matrix_;
}

uint32_t SeedSearcherQueryParameters::GetHashPositionDataListSize() {
	return number_using_hash_position_data_list_;
}

std::vector<AlphabetCoder::Code *> &SeedSearcherQueryParameters::GetQuerySequences() {
	return query_sequences_;
}

std::vector<SeedSearcherQueryParameters::HashPositionDataList> &SeedSearcherQueryParameters::GetHashPositionDataLists() {
	return hash_position_data_lists_;
}

