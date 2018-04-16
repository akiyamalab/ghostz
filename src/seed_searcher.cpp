/*
 * seed_search.cpp
 *
 *  Created on: 2013/02/07
 *      Author: shu
 */

#include <vector>
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <algorithm>
#include <iterator>
#include <assert.h>
#include <float.h>
#include <cmath>
#include "protein_type.h"
#include "alphabet_coder.h"
#include "seed_searcher.h"

#include <fstream>
#include <sstream>
#include "reduced_alphabet_coder.h"
#include "reduced_alphabet_file_reader.h"

#include <omp.h>

SeedSearcher::SeedSearcher() :
		number_threads_(1), query_subsequences_id_(0), query_subsequences_ids_start_(
				0), query_subsequences_ids_end_(0), database_parameters_ptr_(
				NULL), query_parameters_ptr_(NULL), new_hits_(0)
#if 1
				, hash_hits_count(0), calculate_distance_count(0), filter_out_from_distance_count(
				0), representation_hit_count(0), representation_no_similarity_check_hit_count(
				0), member_hit_count(0), member_no_similarity_check_hit_count(
				0), ungapped_extension_count(0)
#endif
{
}
SeedSearcher::~SeedSearcher() {
}

uint32_t SeedSearcher::Reset(SeedSearcherQueryParameters &query_parameters,
		SeedSearcherDatabaseParameters &database_parameters,
		uint32_t k_mer_hit_count_threshold) {

	database_parameters_ptr_ = &database_parameters;
	query_parameters_ptr_ = &query_parameters;
	query_subsequences_ids_start_ = query_subsequences_ids_end_;
	query_sequence_i_ = 0;
	return 0;
}

bool SeedSearcher::Search(std::vector<std::vector<Hit> > &hits_list) {

	const float hit_similarity_threshold_rate = 0.8f;
	const AlphabetCoder::Code sequence_delimiter =
			query_parameters_ptr_->GetSequenceDelimiter();
	const std::vector<int> &ungapped_extension_cutoffs =
			query_parameters_ptr_->GetUngappedExtensionCutoffs();
	const std::vector<int> &gapped_extension_triggers =
			query_parameters_ptr_->GetGappedExtensionTriggers();
	const ScoreMatrix &score_matrix = query_parameters_ptr_->GetScoreMatrix();
	const AlphabetCoder::Code *database_sequence =
			database_parameters_ptr_->GetSequence();
	const SeedSearcherDatabaseParameters::SequencesIndex &sequences_index =
			database_parameters_ptr_->GetSequencesIndex();
	const vector<SeedSearcherQueryParameters::HashPositionDataList> &hash_position_data_lists =
			query_parameters_ptr_->GetHashPositionDataLists();
	const uint32_t hash_position_data_list_size =
			query_parameters_ptr_->GetHashPositionDataListSize();
	const vector<AlphabetCoder::Code *> &query_sequences =
			query_parameters_ptr_->GetQuerySequences();
	const std::vector<AlphabetCoder::Code> &redueced_code_map =
			database_parameters_ptr_->GetReduecedCodeMap();
	const uint32_t subsequence_length =
			database_parameters_ptr_->GetSubsequenceLength();
	const SeedSearcherDatabaseParameters::ClusteringSequencesIndex &clustering_sequences_index =
			database_parameters_ptr_->GetClusteringSequencesIndex();

	int max_representation_member_distance =
			database_parameters_ptr_->GetNumberMaxMismatchs();
	Similarity distance_threshold = subsequence_length
			- hit_similarity_threshold_rate * subsequence_length;
#ifdef _OPENMP
	omp_set_num_threads(number_threads_);
#endif // _OPENMP
	vector<vector<pair<uint32_t, Hit> > > tmp_hits_list(number_threads_);

#pragma omp parallel
	{
		int thread_id = 0;
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#endif // _OPENMP
		vector<uint32_t> hitting_query_position_data_i_list(0);
#pragma omp for schedule(dynamic, 10)
		for (uint32_t hash_position_data_list_i = 0;
				hash_position_data_list_i < hash_position_data_list_size;
				++hash_position_data_list_i) {

			const SeedSearcherQueryParameters::HashPositionDataList &hash_position_data_list =
					hash_position_data_lists[hash_position_data_list_i];

			SearchFromClusteringSubsequences(thread_id, query_sequences,
					database_sequence, sequence_delimiter,
					hash_position_data_list, clustering_sequences_index,
					subsequence_length, redueced_code_map,
					max_representation_member_distance, distance_threshold,
					score_matrix, ungapped_extension_cutoffs,
					gapped_extension_triggers, hits_list,
					tmp_hits_list[thread_id],
					hitting_query_position_data_i_list);

			SearchFromNonClusteringSubsequences(thread_id, query_sequences,
					database_sequence, sequence_delimiter,
					hash_position_data_list, sequences_index, score_matrix,
					ungapped_extension_cutoffs, gapped_extension_triggers,
					hits_list, tmp_hits_list[thread_id]);

		}

	}

	if (number_threads_ > 0) {
		for (size_t i = 1; i < number_threads_; ++i) {
			vector<pair<uint32_t, Hit> > &hits = tmp_hits_list[i];
			for (vector<pair<uint32_t, Hit> >::const_iterator it = hits.begin();
					it != hits.end(); ++it) {
				hits_list[it->first].push_back(it->second);
			}
		}
	}

	return 0;
}

int SeedSearcher::SearchFromClusteringSubsequences(const int thread_id,
		const vector<AlphabetCoder::Code *> &query_sequences,
		const AlphabetCoder::Code *database_sequence,
		const AlphabetCoder::Code sequence_delimiter,
		const SeedSearcherQueryParameters::HashPositionDataList &hash_position_data_list,
		const SeedSearcherDatabaseParameters::ClusteringSequencesIndex &clustering_sequences_index,
		const uint32_t subsequence_length,
		const std::vector<AlphabetCoder::Code> &redueced_code_map,
		const Similarity max_representation_member_distance,
		const Similarity distance_threshold, const ScoreMatrix &score_matrix,
		const std::vector<int> &ungapped_extension_cutoffs,
		const std::vector<int> &gapped_extension_triggers,
		std::vector<std::vector<Hit> > &hits_list,
		vector<pair<uint32_t, Hit> > &tmp_hits,
		std::vector<uint32_t> &hitting_query_position_data_i_list) {


#if 0
// debug ///////////////////
	uint32_t x_query_position = 18;
	uint32_t x_db_position = 888598284;
	//cout << "hash : " << hash_position_data_list.hash << endl;
////////////////////////////
#endif

	size_t representatives_length = 0;
	const SeedSearcherDatabaseParameters::ClusteringSequencesIndex::RepresentativeData *representatives;
	if (database_parameters_ptr_->IsClustering()
			&& !clustering_sequences_index.GetRepresentatives(
					hash_position_data_list.hash, &representatives,
					&representatives_length)) {
#if DEBUG 
		hash_hits_count += representatives_length
		* hash_position_data_list.position_data.size();
#endif

		for (size_t representatives_i = 0;
				representatives_i < representatives_length;
				++representatives_i) {
			hitting_query_position_data_i_list.clear();
			uint32_t representative_position =
					representatives[representatives_i].second;
			size_t number_query_position_data =
					hash_position_data_list.position_data.size();
			const AlphabetCoder::Code *representative_subsequence_center =
					database_sequence + representative_position;
			for (size_t query_position_data_i = 0;
					query_position_data_i < number_query_position_data;
					++query_position_data_i) {
				const SeedSearcherQueryParameters::HashPositionData &query_position_data =
						hash_position_data_list.position_data[query_position_data_i];
				int query_representation_distance = 0;
				if (query_position_data.similarity_check) {
#if DEBUG 
					++calculate_distance_count;
#endif

					query_representation_distance =
							SeedSearcherCommon::CalculateDistance(
									query_sequences[query_position_data.query_id]
											+ query_position_data.position,
									representative_subsequence_center,
									subsequence_length, redueced_code_map);

					if (IsSufficientSimilarityCluster(
							query_representation_distance,
							max_representation_member_distance,
							distance_threshold)) {
						hitting_query_position_data_i_list.push_back(
								query_position_data_i);
#if 0
					if (x_db_position == representative_position) {
						cout << "sim pass " << ": q "
								<< query_position_data.position
								<< ", d " << representative_position << " dis "
								<< query_representation_distance << endl;
					}
#endif
					}
				} else {
					hitting_query_position_data_i_list.push_back(
							query_position_data_i);
				}
				if (distance_threshold >= query_representation_distance) {
#if DEBUG 
					if (query_position_data.similarity_check) {
						++representation_hit_count;
					} else {
						++representation_no_similarity_check_hit_count;
					}
#endif

#if 0
					if (x_db_position == representative_position) {
						cout << "sim r ok" << ": q "
								<< query_position_data.position
								<< ", d " << representative_position << " dis "
								<< query_representation_distance << endl;
					}
#endif

					int score =
							UngappedExtend(
									query_sequences[query_position_data.query_id],
									query_position_data.position,
									database_sequence, representative_position,
									sequence_delimiter, score_matrix,
									ungapped_extension_cutoffs[query_position_data.query_id],
									gapped_extension_triggers[query_position_data.query_id]);
					if (score
							> gapped_extension_triggers[query_position_data.query_id]) {
						Hit new_hit = BuildHit(query_position_data.position,
								representative_position);
						if (thread_id == 0) {
							hits_list[query_position_data.query_id].push_back(
									new_hit);
						} else {
							tmp_hits.push_back(
									make_pair(query_position_data.query_id,
											new_hit));
						}
					}
				}
			}

			const SeedSearcherDatabaseParameters::ClusteringSequencesIndex::Position *members;
			size_t members_length = 0;
			if (!hitting_query_position_data_i_list.empty()
					&& !clustering_sequences_index.GetMembers(
							representatives[representatives_i].first, &members,
							&members_length)) {
				for (size_t members_i = 0; members_i < members_length;
						++members_i) {
					for (size_t i = 0;
							i < hitting_query_position_data_i_list.size();
							++i) {
						const SeedSearcherQueryParameters::HashPositionData &query_position_data =
								hash_position_data_list.position_data[hitting_query_position_data_i_list[i]];
#if DEBUG 
						if (query_position_data.similarity_check) {
							++member_hit_count;
							++hash_hits_count;
						} else {
							++member_no_similarity_check_hit_count;
							++hash_hits_count;
						}
#endif
#if 0
					if (x_db_position == members[members_i]) {
						cout << "sim m ok " << ": q "
								<< query_position_data.position
								<< ", d " << members[members_i] << endl;
					}
#endif
						int score =
								UngappedExtend(
										query_sequences[query_position_data.query_id],
										query_position_data.position,
										database_sequence, members[members_i],
										sequence_delimiter, score_matrix,
										ungapped_extension_cutoffs[query_position_data.query_id],
										gapped_extension_triggers[query_position_data.query_id]);
						if (score
								> gapped_extension_triggers[query_position_data.query_id]) {
							Hit new_hit = BuildHit(query_position_data.position,
									members[members_i]);
							if (thread_id == 0) {
								hits_list[query_position_data.query_id].push_back(
										new_hit);
							} else {
								tmp_hits.push_back(
										make_pair(query_position_data.query_id,
												new_hit));
							}
						}
					}
				}
			}

		}

	}
	return 0;
}

int SeedSearcher::SearchFromNonClusteringSubsequences(const int thread_id,
		const vector<AlphabetCoder::Code *> &query_sequences,
		const AlphabetCoder::Code *database_sequence,
		const AlphabetCoder::Code sequence_delimiter,
		const SeedSearcherQueryParameters::HashPositionDataList &hash_position_data_list,
		const SeedSearcherDatabaseParameters::SequencesIndex &sequences_index,
		const ScoreMatrix &score_matrix,
		const std::vector<int> &ungapped_extension_cutoffs,
		const std::vector<int> &gapped_extension_triggers,
		std::vector<std::vector<Hit> > &hits_list,
		vector<pair<uint32_t, Hit> > &tmp_hits) {


#if 0
// debug ///////////////////
	uint32_t x_query_position = 18;
	uint32_t x_db_position = 888598284;
	//cout << "hash : " << hash_position_data_list.hash << endl;
////////////////////////////
#endif

	const SeedSearcherDatabaseParameters::SequencesIndex::Value *values = NULL;
	size_t values_length = 0;
	if (!sequences_index.GetValues(hash_position_data_list.hash, &values,
			&values_length)) {
#if DEBUG
		hash_hits_count += values_length
		* hash_position_data_list.position_data.size();
#endif
		for (size_t values_i = 0; values_i < values_length; ++values_i) {
			for (std::vector<SeedSearcherQueryParameters::HashPositionData>::const_iterator queries_position_data_it =
					hash_position_data_list.position_data.begin();
					queries_position_data_it
							!= hash_position_data_list.position_data.end();
					++queries_position_data_it) {
#if 0
				Similarity similarity = CalculateSimilarity(
						query_sequences[queries_position_data_it->query_id]
						+ queries_position_data_it->query_position,
						database_sequence + values[values_i],
						subsequence_length / 2, redueced_code_map);
				if (similarity <= hit_similarity_threshold) {
					continue;
				}
#endif
#if 0
					if (x_db_position ==  values[values_i]) {
						cout << "v ok " << ": q "
								<< queries_position_data_it->position
								<< ", d " << values[values_i] << endl;
					}
#endif
				int score =
						UngappedExtend(
								query_sequences[queries_position_data_it->query_id],
								queries_position_data_it->position,
								database_sequence, values[values_i],
								sequence_delimiter, score_matrix,
								ungapped_extension_cutoffs[queries_position_data_it->query_id],
								gapped_extension_triggers[queries_position_data_it->query_id]);
				if (score
						> gapped_extension_triggers[queries_position_data_it->query_id]) {
					Hit new_hit = BuildHit(queries_position_data_it->position,
							values[values_i]);
					if (thread_id == 0) {
						hits_list[queries_position_data_it->query_id].push_back(
								new_hit);
					} else {
						tmp_hits.push_back(
								make_pair(queries_position_data_it->query_id,
										new_hit));
					}
				}
			}
		}
	}

	return 0;
}

SeedSearcher::Hit SeedSearcher::BuildHit(uint32_t query_position,
		uint32_t database_position) const {
	Hit new_hit;
	new_hit.query_sequence_position = query_position;
	new_hit.database_sequence_position = database_position;
	new_hit.k_mer_count = 1;
//hits_list[query_id].push_back(new_hit);
	return new_hit;
}

bool SeedSearcher::IsSufficientSimilarityCluster(
		const int query_representation_distance,
		const int max_representation_member_distance,
		const int distance_threshold) const {

	int lower_bound = max(
			query_representation_distance - max_representation_member_distance,
			0);
	return distance_threshold >= lower_bound;
}

int SeedSearcher::DumpSearchLog() {
	cout << "hash hits count " << hash_hits_count << endl;
	cout << "calculate_distance_count " << calculate_distance_count << endl;
	cout << "representation_hit_count " << representation_hit_count << endl;
	cout << "member_hit_count " << member_hit_count << endl;
	cout << "representation_no_similarity_check_hit_count "
			<< representation_no_similarity_check_hit_count << endl;
	cout << "member_no_similarity_check_hit_count "
			<< member_no_similarity_check_hit_count << endl;
	cout << "ungapped_extension_count " << ungapped_extension_count << endl;
	return 0;
}
