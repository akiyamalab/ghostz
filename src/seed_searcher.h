/*
 * seed_searcher.h
 *
 *  Created on: 2013/02/07
 *      Author: shu
 */

#ifndef SEED_SEARCHER_H_
#define SEED_SEARCHER_H_

#include "alphabet_coder.h"
#include "reduced_alphabet_k_mer_hash_function.h"
#include "seed_searcher_database_parameters.h"
#include "seed_searcher_query_parameters.h"
#include "ungapped_extender.h"
#include <stdint.h>
#include <list>
#include <iostream>
#include <iomanip>
#include <map>
#include <tr1/memory>

class SeedSearcher {
public:
	typedef SeedSearcherCommon::Distance Similarity;
	static SeedSearcher::Similarity CalculateSimilarity(
			const AlphabetCoder::Code *query_sequence,
			const AlphabetCoder::Code *database_sequence,
			const uint32_t max_one_direction_length,
			const std::vector<AlphabetCoder::Code> &redueced_code_map);

	static SeedSearcher::Similarity GetSeedPosition(
			const AlphabetCoder::Code *subsequence,

			const AlphabetCoder::Code *database_sequence,
			const uint32_t max_one_direction_length,
			const std::vector<AlphabetCoder::Code> &redueced_code_map);

	typedef struct {
		uint32_t query_sequence_position;
		uint32_t database_sequence_position;
		int k_mer_count;
	} Hit;

	SeedSearcher();
	virtual ~SeedSearcher();

	void SetNumberThreads(uint32_t number_threads) {
		number_threads_ = number_threads;
	}

	uint32_t Reset(SeedSearcherQueryParameters &query_parameters,
			SeedSearcherDatabaseParameters &database_parameters,
			uint32_t k_mer_hit_count_threshold);

	bool Search(std::vector<std::vector<Hit> > &hits_list);
	int DumpSearchLog();

private:
	int SearchFromClusteringSubsequences(const int thread_id,
			const std::vector<AlphabetCoder::Code *> &query_sequences,
			const AlphabetCoder::Code *database_sequence,
			const AlphabetCoder::Code sequence_delimiter,
			const SeedSearcherQueryParameters::HashPositionDataList &hash_position_data_list,
			const SeedSearcherDatabaseParameters::ClusteringSequencesIndex &clustering_sequences_index,
			const uint32_t subsequence_length,
			const std::vector<AlphabetCoder::Code> &redueced_code_map,
			const Similarity representation_member_similarity,
			const Similarity hit_similarity_threshold,
			const ScoreMatrix &score_matrix,
			const std::vector<int> &ungapped_extension_cutoffs,
			const std::vector<int> &gapped_extension_triggers,
			std::vector<std::vector<Hit> > &hits_list,
			std::vector<std::pair<uint32_t, Hit> > &tmp_hits,
			std::vector<uint32_t> &hitting_queries);
	int SearchFromNonClusteringSubsequences(const int thread_id,
			const std::vector<AlphabetCoder::Code *> &query_sequences,
			const AlphabetCoder::Code *database_sequence,
			const AlphabetCoder::Code sequence_delimiter,
			const SeedSearcherQueryParameters::HashPositionDataList &hash_position_data_list,
			const SeedSearcherDatabaseParameters::SequencesIndex &sequences_index,
			const ScoreMatrix &score_matrix,
			const std::vector<int> &ungapped_extension_cutoffs,
			const std::vector<int> &gapped_extension_triggers,
			std::vector<std::vector<Hit> > &hits_list,
			std::vector<std::pair<uint32_t, Hit> > &tmp_hits);
	int UngappedExtend(const AlphabetCoder::Code *query_sequence,
			uint32_t query_position,
			const AlphabetCoder::Code *database_sequence,
			uint32_t database_position, AlphabetCoder::Code sequence_delimiter,
			const ScoreMatrix &score_matrix, int cutoff,
                     int gapped_extension_trigger); // __attribute__((noinline));
	bool CheckSimilarity(AlphabetCoder::Code *query_sequence,
			uint32_t query_position, AlphabetCoder::Code *database_sequence,
			uint32_t database_position, uint32_t subsequence_length,
			std::vector<AlphabetCoder::Code> &redueced_code_map,
			Similarity representation_member_similarity,
			Similarity hit_similarity_threshold) const;
	Hit BuildHit(uint32_t query_position, uint32_t database_position) const;
	bool IsSufficientSimilarityCluster(const int query_representation_distance,
			const int representation_member_distance,
			const int distance_threshold) const;

	uint32_t number_threads_;
	uint32_t query_sequence_i_;
	size_t query_subsequences_id_;
	size_t query_subsequences_ids_start_;
	size_t query_subsequences_ids_end_;
	SeedSearcherDatabaseParameters *database_parameters_ptr_;
	SeedSearcherQueryParameters *query_parameters_ptr_;
	std::vector<Hit> new_hits_;

#if 1
	uint64_t hash_hits_count;
	uint64_t calculate_distance_count;
	uint64_t filter_out_from_distance_count;
	uint64_t representation_hit_count;
	uint64_t representation_no_similarity_check_hit_count;
	uint64_t member_hit_count;
	uint64_t member_no_similarity_check_hit_count;
	uint64_t ungapped_extension_count;
#endif

};

inline SeedSearcher::Similarity SeedSearcher::CalculateSimilarity(
		const AlphabetCoder::Code *query_sequence,
		const AlphabetCoder::Code *database_sequence,
		const uint32_t max_one_direction_length,
		const std::vector<AlphabetCoder::Code> &redueced_code_map) {
	int mismatch_count = 0;
	int extension_offset = -1;
	int extension_end = -max_one_direction_length + extension_offset;
	for (int i = extension_offset; i > extension_end; --i) {
#if 0
		std::cout << int(redueced_code_map[query_sequence[i]]) << ","
		<< int(redueced_code_map[database_sequence[i]]) << std::endl;
#endif
		mismatch_count +=
				(redueced_code_map[query_sequence[i]]
						!= redueced_code_map[database_sequence[i]]) ? 1 : 0;
	}

	extension_offset = 0;
	extension_end = max_one_direction_length;
	for (int i = extension_offset; i < extension_end; ++i) {
#if 0
		std::cout << int(redueced_code_map[query_sequence[i]]) << ","
		<< int(redueced_code_map[database_sequence[i]]) << std::endl;
#endif
		mismatch_count +=
				(redueced_code_map[query_sequence[i]]
						!= redueced_code_map[database_sequence[i]]) ? 1 : 0;
	}

#if 0
	std::cout << "length : " << 2*max_one_direction_length << ", mismatch : " << mismatch_count
	<< std::endl;
#endif

	return max_one_direction_length * 2 - mismatch_count;

}

inline int SeedSearcher::UngappedExtend(
		const AlphabetCoder::Code *query_sequence, uint32_t query_position,
		const AlphabetCoder::Code *database_sequence,
		uint32_t database_position, AlphabetCoder::Code sequence_delimiter,
		const ScoreMatrix &score_matrix, int cutoff,
		int gapped_extension_trigger) {
	//++ungapped_extension_count;

	int sum_score = 0;
	int score = 0;
	UngappedExtender::ExtendOneSideScoreOnly(
			query_sequence + query_position - 1,
			database_sequence + database_position - 1, sequence_delimiter, true,
			score_matrix, cutoff, &score);
	sum_score += score;

	if (sum_score <= gapped_extension_trigger) {
		UngappedExtender::ExtendOneSideScoreOnly(
				query_sequence + query_position,
				database_sequence + database_position, sequence_delimiter,
				false, score_matrix, cutoff, &score);
		sum_score += score;
	}

	return sum_score;
}
#endif /* SEED_SEARCHER_H_ */
