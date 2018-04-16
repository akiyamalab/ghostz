/*
 * aligner.cpp
 *
 *  Created on: 2012/10/05
 *      Author: shu
 */

#include "aligner.h"
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include "limits.h"
#include "logger.h"
#include "alphabet_coder.h"
#include "chain_filter.h"
#include "edit_blocks.h"
#include "score_matrix_reader.h"
#include "score_matrix.h"
#include "alphabet_type.h"
#include "protein_type.h"
#include "dna_type.h"
#include "ungapped_extender.h"
#include "gapped_extender.h"
#include "seed_searcher.h"
#include "seed_searcher_query_parameters.h"
#include "reduced_alphabet_coder.h"

#include <omp.h>

using namespace std;

void Aligner::BuildDatabase(string &input_filename, string &database_filename,
		DatabaseParameters paramters) {
	Logger *logger = Logger::GetInstance();
	Database::Parameters database_parameters;
	BuildDatabaseParameters(paramters, database_parameters);
	ifstream is(input_filename.c_str());
	if (!is) {
		logger->ErrorLog("invalid input filename");
	}
	Database database(is, database_parameters);
	if (database.GetChunkId() < 0) {
		logger->ErrorLog("invalid input file or parameters");
	} else {
		database.Save(database_filename);
		cout << "number chunks is " << database.GetNumberChunks() << endl;
		cout << "number of sequences is " << database.GetNumberTotalSequences()
				<< endl;
		cout << "total size is " << database.GetDatabaseTotalLenght() << endl;
	}
	is.close();
}

void Aligner::Align(string &queries_filename, string &database_filename,
		string &output_filename, AligningParameters &parameters) {
	Queries::Parameters queries_parameters;
	ifstream queries_is(queries_filename.c_str());
	BuildQueriesParameters(parameters, queries_parameters);
	Database database(database_filename);
	ofstream os(output_filename.c_str());
	for (Queries queries(queries_is, queries_parameters);
			queries.GetNumberSequences() != 0; queries.Next()) {
		cout << "number queries is " << queries.GetNumberSequences() << endl;
		vector<vector<Aligner::PresearchedResult> > presearch_results_list(
				queries.GetNumberSequences());
		vector<vector<Aligner::Result> > results_list(
				queries.GetNumberSequences());
		cout << "starts presearch " << endl;
		Presearch(queries, database, parameters, presearch_results_list);
		cout << "starts build results" << endl;
		BuildResults(queries, database, parameters, presearch_results_list,
				results_list);
		cout << "writes results" << endl;
		WriteOutput(os, queries, database, parameters, results_list);
	}
	queries_is.close();
	os.close();
}

int Aligner::AlignmentComp(const int a_score, const Coordinate &a_start,
		const Coordinate &a_end, const int b_score, const Coordinate &b_start,
		const Coordinate &b_end) {
	if (a_score > b_score) {
		return 1;
	} else if (a_score < b_score) {
		return -1;
	}

	uint32_t a_q_length = a_end.query_position - a_start.query_position;
	uint32_t b_q_length = b_end.query_position - b_start.query_position;

	if (a_q_length > b_q_length) {
		return 1;
	} else if (a_q_length < b_q_length) {
		return -1;
	}

	uint32_t a_s_length = a_end.database_position - a_start.database_position;
	uint32_t b_s_length = b_end.database_position - b_start.database_position;

	if (a_s_length > b_s_length) {
		return 1;
	} else if (a_s_length < b_s_length) {
		return -1;
	}

#if 1
	if (a_start.database_position < b_start.database_position) {
		return 1;
	} else {
		return -1;
	}
	if (a_start.query_position < b_start.query_position) {
		return 1;
	} else {
		return -1;
	}
#endif

	return 0;
}

bool Aligner::IsContained(Coordinate &a_start, Coordinate &a_end,
		Coordinate &b_start, Coordinate &b_end) {
	if ((a_start.query_position == b_start.query_position
			&& a_start.database_position == b_start.database_position) // overlapped start
			|| (a_end.query_position == b_end.query_position
					&& a_end.database_position == b_end.database_position) // overlapped end
			|| (a_start.query_position <= b_start.query_position // a contain b
			&& b_end.query_position <= a_end.query_position
					&& a_start.database_position <= b_start.database_position
					&& b_end.database_position <= a_end.database_position)) {
		return true;
	} else {
		return false;
	}
}

AlphabetCoder::Code Aligner::GetSequenceDelimiter(AlphabetType &type) {
	AlphabetCoder coder(type);
	return coder.GetMaxCode() + 1;
}

void Aligner::BuildQueriesParameters(AligningParameters &parameters,
		Queries::Parameters &queries_parameters) {
	Statistics statistics(*parameters.aligning_sequence_type_ptr);
	queries_parameters.filter = parameters.filter;
	queries_parameters.file_sequence_type_ptr =
			parameters.queries_file_sequence_type_ptr;
	queries_parameters.aligning_sequence_type_ptr =
			parameters.aligning_sequence_type_ptr;
	statistics.CalculateUngappedIdealKarlinParameters(parameters.score_matrix,
			&queries_parameters.ungapped_karlin_parameters);
	queries_parameters.chunk_size = parameters.queries_chunk_size;
	queries_parameters.score_matrix = parameters.score_matrix;
	queries_parameters.sequence_delimiter = GetSequenceDelimiter(
			*parameters.aligning_sequence_type_ptr);
}
void Aligner::BuildDatabaseParameters(DatabaseParameters &parameters,
		Database::Parameters &database_parameters) {
	database_parameters.chunk_size = parameters.chunk_size;
	database_parameters.chunk_build_option = parameters.chunk_build_option;
	database_parameters.sequence_type_ptr = parameters.sequence_type_ptr;
	database_parameters.seed_search_parameters_build_parameters.number_threads =
			parameters.number_threads;
	database_parameters.seed_search_parameters_build_parameters.seed_threshold =
			parameters.seed_threshold;
	database_parameters.seed_search_parameters_build_parameters.score_matrix =
			parameters.score_matrix;
	database_parameters.sequence_delimiter = GetSequenceDelimiter(
			*parameters.sequence_type_ptr);
	database_parameters.seed_search_parameters_build_parameters.sequence_delimiter =
			database_parameters.sequence_delimiter;
	database_parameters.seed_search_parameters_build_parameters.clustering =
			parameters.clustering;
	database_parameters.seed_search_parameters_build_parameters.subsequence_length =
			parameters.clustering_subsequence_length;
	AlphabetCoder coder(*parameters.sequence_type_ptr);
	if (parameters.hash_alphabet_sets.empty()) {
		database_parameters.seed_search_parameters_build_parameters.max_indexing_code =
				coder.GetMaxRegularLetterCode();
	} else {
		ReducedAlphabetCoder hash_reduced_alphabet_coder(
				*(parameters.sequence_type_ptr), parameters.hash_alphabet_sets);
		AlphabetCoder::Code max_code = coder.GetMaxCode();
		vector<AlphabetCoder::Code> &hash_reduced_code_map =
				database_parameters.seed_search_parameters_build_parameters.hash_code_map;
		hash_reduced_code_map.resize(coder.GetMaxCode() + 1);
		for (AlphabetCoder::Code code = coder.GetMinCode(); code <= max_code;
				++code) {
			char c = coder.Decode(code);
			AlphabetCoder::Code reduced_code =
					hash_reduced_alphabet_coder.Encode(c);
			hash_reduced_code_map[code] = reduced_code;
		}
		database_parameters.seed_search_parameters_build_parameters.max_indexing_code =
				hash_reduced_alphabet_coder.GetMaxRegularLetterCode();

		ReducedAlphabetCoder similarity_reduced_alphabet_coder(
				*(parameters.sequence_type_ptr),
				parameters.similarity_alphabet_sets);
		vector<AlphabetCoder::Code> &similairty_reduced_code_map =
				database_parameters.seed_search_parameters_build_parameters.similairty_code_map;
		similairty_reduced_code_map.resize(coder.GetMaxCode() + 1);
		for (AlphabetCoder::Code code = coder.GetMinCode(); code <= max_code;
				++code) {
			char c = coder.Decode(code);
			AlphabetCoder::Code reduced_code =
					similarity_reduced_alphabet_coder.Encode(c);
			similairty_reduced_code_map[code] = reduced_code;
		}
	}
}

void Aligner::Presearch(Queries &queries, Database &database,
		AligningParameters &parameters,
		std::vector<std::vector<PresearchedResult> > &results_list) {

	Statistics statistics(*(parameters.aligning_sequence_type_ptr));
	Statistics::KarlinParameters gapped_karlin_parameters;
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension,
			&gapped_karlin_parameters);
	int gapped_extension_cutoff = Statistics::Normalized2Nominal(
			parameters.normalized_presearched_gapped_extension_cutoff,
			gapped_karlin_parameters);
	SeedSearcher searcher;
	searcher.SetNumberThreads(parameters.number_threads);

	vector<vector<SeedSearcher::Hit> > hits_list(queries.GetNumberSequences());
	vector<vector<PresearchedResult> > tmp_results_list(
			parameters.number_threads);
	vector<vector<std::vector<AlignmentPosition> > > alignment_start_positions_in_database_list(
			parameters.number_threads);
	std::vector<int> ungapped_extension_cutoffs(queries.GetNumberSequences());
	std::vector<int> gapped_extension_triggers(queries.GetNumberSequences());
	for (uint32_t query_id = 0; query_id < queries.GetNumberSequences();
			++query_id) {
		Query *query = queries.GetQuery(query_id);
		ungapped_extension_cutoffs[query_id] =
				Statistics::NormalizedCutoff2NominalCutoff(
						parameters.normalized_presearched_ungapped_extension_cutoff,
						query->GetUngappedKarlinParameters());
		gapped_extension_triggers[query_id] = Statistics::Normalized2Nominal(
				parameters.normalized_presearched_gapped_extension_trigger,
				query->GetUngappedKarlinParameters());

	}
	SeedSearcherQueryParameters::BuildParameters seed_search_query_parameters_build_parameter;
	seed_search_query_parameters_build_parameter.subsequence_length =
			database.GetSeedSearcherParameters().GetSubsequenceLength();
	seed_search_query_parameters_build_parameter.sequence_delimiter =
			queries.GetQuery(0)->GetSequenceDelimiter();
	seed_search_query_parameters_build_parameter.ungapped_extension_cutoffs =
			&ungapped_extension_cutoffs;
	seed_search_query_parameters_build_parameter.gapped_extension_triggers =
			&gapped_extension_triggers;
	seed_search_query_parameters_build_parameter.score_matrix =
			&parameters.score_matrix;
	seed_search_query_parameters_build_parameter.sequences_index =
			&database.GetSeedSearcherParameters().GetSequencesIndex();
	SeedSearcherQueryParameters seed_search_query_parameter;
	seed_search_query_parameter.Build(queries,
			seed_search_query_parameters_build_parameter);

#if DEBUG 
	uint32_t total_gapped_extention_count = 0;
#endif
	for (database.ResetChunk();
			database.GetChunkId() < (int) database.GetNumberChunks();
			database.NextChunk()) {
		const AlphabetCoder::Code *database_concatenated_sequence =
				database.GetConcatenatedSequence();
		const uint32_t database_concatenated_sequence_length =
				database.GetConcatenatedSequenceLength();
		searcher.Reset(seed_search_query_parameter,
				database.GetSeedSearcherParameters(), 1);
		searcher.Search(hits_list);
#ifdef _OPENMP
		omp_set_num_threads(parameters.number_threads);
#endif // _OPENMP
#pragma omp parallel
		{

			int thread_id = 0;
#ifdef _OPENMP
			thread_id = omp_get_thread_num();
#endif // _OPENMP
			ChainFilter chain_filter(
					queries.GetQuery(0)->GetSequenceDelimiter(),
					parameters.score_matrix);
			UngappedExtender ungapped_extender;
			GappedExtender gapped_extender;
			vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database =
					alignment_start_positions_in_database_list[thread_id];
			if (alignment_start_positions_in_database.size()
					< database.GetNumberSequencesInChunk()) {
				alignment_start_positions_in_database.resize(
						database.GetNumberSequencesInChunk());
			}
			vector<PresearchedResult> &tmp_results = tmp_results_list[thread_id];
			PresearchedResult tmp_result;
			tmp_result.database_chunk_id = UINT_MAX;
			tmp_result.subject_id_in_chunk = UINT_MAX;
			uint32_t number_queries = queries.GetNumberSequences();
#pragma omp for schedule(dynamic, 10)
			for (uint32_t query_id = 0; query_id < number_queries; ++query_id) {
				Query *query = queries.GetQuery(query_id);
				vector<SeedSearcher::Hit> &hits = hits_list[query_id];
				chain_filter.Filter(query->GetSequence(),
						query->GetSequenceLength(),
						ungapped_extension_cutoffs[query_id],
						database_concatenated_sequence,
						database_concatenated_sequence_length, hits);
				tmp_results.clear();
				for (size_t hit_id = 0; hit_id < hits.size(); ++hit_id) {
#if DEBUG 
					++total_gapped_extention_count;
#endif

					int score = 0;
					int query_position;
					int database_position;
					tmp_result.score = 0;
					tmp_result.hit.query_position =
							hits[hit_id].query_sequence_position;
					tmp_result.hit.database_position =
							hits[hit_id].database_sequence_position;

					gapped_extender.ExtendOneSide(
							query->GetSequence() + tmp_result.hit.query_position
									- 1, tmp_result.hit.query_position - 1,

							database.GetConcatenatedSequence()
									+ tmp_result.hit.database_position - 1,
							query->GetSequenceDelimiter(), true,
							parameters.score_matrix, parameters.gap_open,
							parameters.gap_extension, gapped_extension_cutoff,
							&score, &query_position, &database_position, NULL);

					tmp_result.score += score;
					tmp_result.start.query_position =
							tmp_result.hit.query_position - 1 + query_position;
					tmp_result.start.database_position =
							tmp_result.hit.database_position - 1
									+ database_position;

					gapped_extender.ExtendOneSide(
							query->GetSequence()
									+ tmp_result.hit.query_position,
							query->GetSequenceLength()
									- tmp_result.hit.query_position,
							database.GetConcatenatedSequence()
									+ tmp_result.hit.database_position,
							query->GetSequenceDelimiter(), false,
							parameters.score_matrix, parameters.gap_open,
							parameters.gap_extension, gapped_extension_cutoff,
							&score, &query_position, &database_position, NULL);
					tmp_result.score += score;
					tmp_result.end.query_position += query_position + 1;
					tmp_result.end.database_position += database_position + 1;

					tmp_result.hit_count = hits[hit_id].k_mer_count;
					tmp_results.push_back(tmp_result);

				}

				AddResults(database, parameters,
						alignment_start_positions_in_database, tmp_results,
						results_list[query_id]);
				hits.clear();

#if DEBUG
				cout << query->GetName() << "\t" << k_mers_hits_count << "\t" << hits_count << "\t" <<gapped_extentio
				n_count << endl;
#endif
			}

		}

	}

#if DEBUG 
	cout << "Searcher Dumplog : " << endl;
	searcher.DumpSearchLog();
	cout << "total_gapped_extention_count\t" << total_gapped_extention_count
	<< endl;
#endif
}

void Aligner::BuildResults(Queries &queries, Database &database,
		AligningParameters &parameters,
		std::vector<std::vector<PresearchedResult> > &presearch_results_list,
		std::vector<std::vector<Result> > &results_list) {

	Statistics statistics(*(parameters.aligning_sequence_type_ptr));
	Statistics::KarlinParameters gapped_karlin_parameters;
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension,
			&gapped_karlin_parameters);
	const int gapped_extension_cutoff = Statistics::Normalized2Nominal(
			parameters.normalized_result_gapped_extension_cutoff,
			gapped_karlin_parameters);

	vector<vector<pair<uint32_t, uint32_t> > > result_ids_list(
			database.GetNumberChunks());
	if (presearch_results_list.size() < queries.GetNumberSequences()) {
		presearch_results_list.resize(queries.GetNumberSequences());
	}
	for (uint32_t query_id = 0; query_id < queries.GetNumberSequences();
			++query_id) {
		results_list[query_id].resize(presearch_results_list[query_id].size());
		for (uint32_t result_id = 0;
				result_id < presearch_results_list[query_id].size();
				++result_id) {
			result_ids_list[presearch_results_list[query_id][result_id].database_chunk_id].push_back(
					make_pair(query_id, result_id));
		}
	}

#ifdef _OPENMP
	omp_set_num_threads(parameters.number_threads);
#endif // _OPENMP
	for (database.ResetChunk();
			database.GetChunkId() < database.GetNumberChunks();
			database.NextChunk()) {
		const uint32_t database_chunk_id = database.GetChunkId();
		const AlphabetCoder::Code *database_concatenated_sequence =
				database.GetConcatenatedSequence();
		vector<pair<uint32_t, uint32_t> > &result_ids =
				result_ids_list[database_chunk_id];
		const size_t result_ids_list_size = result_ids.size();
#pragma omp parallel
		{
			EditBlocks edit_blocks;
			EditBlocks tmp_edit_blocks;
			GappedExtender gapped_extender;
#pragma omp for schedule(dynamic, 10)
			for (size_t i = 0; i < result_ids_list_size; ++i) {
				edit_blocks.Clear();
				uint32_t query_id = result_ids[i].first;
				uint32_t result_id = result_ids[i].second;
				Query *query = queries.GetQuery(query_id);
				int sum_score = 0;
				Coordinate hit;
				Coordinate start;
				Coordinate end;
				PresearchedResult *presearched_result_ptr =
						&presearch_results_list[query_id][result_id];
				hit.query_position = presearched_result_ptr->hit.query_position;
				hit.database_position =
						presearched_result_ptr->hit.database_position;
				int score;
				int query_position;
				int database_position;
				gapped_extender.ExtendOneSide(
						query->GetSequence() + hit.query_position - 1,
						hit.query_position - 1,
						database_concatenated_sequence + hit.database_position
								- 1, query->GetSequenceDelimiter(), true,
						parameters.score_matrix, parameters.gap_open,
						parameters.gap_extension, gapped_extension_cutoff,
						&score, &query_position, &database_position,
						&tmp_edit_blocks);
				sum_score += score;
				start.query_position = hit.query_position - 1 + query_position;
				start.database_position = hit.database_position - 1
						+ database_position;
				tmp_edit_blocks.Reverse();
				edit_blocks.Add(tmp_edit_blocks);
				gapped_extender.ExtendOneSide(
						query->GetSequence() + hit.query_position,
						query->GetSequenceLength() - hit.query_position,
						database_concatenated_sequence + hit.database_position,
						query->GetSequenceDelimiter(), false,
						parameters.score_matrix, parameters.gap_open,
						parameters.gap_extension, gapped_extension_cutoff,
						&score, &query_position, &database_position,
						&tmp_edit_blocks);
				sum_score += score;
				end.query_position = hit.query_position + query_position;
				end.database_position = hit.database_position
						+ database_position;
				edit_blocks.Add(tmp_edit_blocks);
				vector<EditBlocks::EditOpType> edits = edit_blocks.ToVector();
				BuildResult(*query, database,
						presearched_result_ptr->database_chunk_id,
						presearched_result_ptr->subject_id_in_chunk, sum_score,
						hit, start, end, edits,
						results_list[query_id][result_id]);
				// debug ///////////
				results_list[query_id][result_id].hit_count =
						presearched_result_ptr->hit_count;
				///////////////////////
			}

		}

	}
	for (uint32_t query_id = 0; query_id < queries.GetNumberSequences();
			++query_id) {
		sort(results_list[query_id].begin(), results_list[query_id].end(),
				ResultGreaterScore());
	}
}

void Aligner::WriteOutput(ostream &os, Queries &queries, Database &database,
		AligningParameters &parameters,
		std::vector<std::vector<Result> > &results) {
	Statistics::KarlinParameters gapped_karlin_parameters;
	double alpha;
	double beta;

	Statistics statistics(*(parameters.aligning_sequence_type_ptr));
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension,
			&gapped_karlin_parameters);
	statistics.CalculateAlphaBeta(parameters.score_matrix, parameters.gap_open,
			parameters.gap_extension, &alpha, &beta);
	for (uint32_t i = 0; i < queries.GetNumberSequences(); ++i) {
		Query *query = queries.GetQuery(i);
		string query_name = query->GetName();
		uint64_t search_space = statistics.CalculateSearchSpace(
				query->GetRealSequenceLength(),
				database.GetDatabaseTotalLenght(),
				database.GetNumberTotalSequences(), gapped_karlin_parameters,
				alpha, beta);
		for (vector<Result>::iterator it = results[i].begin();
				it != results[i].end(); ++it) {
			os << query_name << "\t";
			os << it->subject_name << "\t";
			os
					<< (1.0
							- (static_cast<float>(it->mismatches
									+ it->gap_openings + it->gap_extensions)
									/ static_cast<float>(it->alignment_length)))
							* 100 << "\t";
			os << it->alignment_length << "\t";
			os << it->mismatches << "\t";
			os << it->gap_openings << "\t";
			os << it->start.query_position << "\t";
			os << it->end.query_position << "\t";
			os << it->start.database_position << "\t";
			os << it->end.database_position << "\t";
			os
					<< Statistics::Nominal2EValue(it->score, search_space,
							gapped_karlin_parameters) << "\t";
			os
					<< Statistics::Nominal2Normalized(it->score,
							gapped_karlin_parameters);
			//os << "\t" << it->score << "\t" << it->hit_count << "\t";
			//os << "\t(" << it->hit.query_position << ", "
			//		<< it->hit.database_position << ")";
			os << '\n';
		}
	}
}

void Aligner::AddResults(Database &database, AligningParameters &parameters,
		std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database,
		std::vector<PresearchedResult> &added_results,
		std::vector<PresearchedResult> &results) {
	vector<PresearchedResult> new_results(0);
	if (alignment_start_positions_in_database.size()
			< database.GetNumberSequencesInChunk()) {
		alignment_start_positions_in_database.resize(
				database.GetNumberSequencesInChunk());
	}
	stable_sort(added_results.begin(), added_results.end(),
			PresearchedResultGreaterScore());

#if 0
// debug //////////////////////////////////////////////////
	cout << "sorted hits" << endl;
	for (uint32_t x = 0; x < hits.size(); ++x) {
		cout << hits[x].score << " "<<hits[x].start.query_position << " " << hits[x].start.database_position << " " << hits[x].end.query_position << " " << hits[x].end.database_position<< " " << hits[x].hit.query_position << " " << hits[x].hit.database_position << endl;
	}
	cout << endl;
///////////////////////////////////////////////////////////
#endif

	vector<PresearchedResult>::iterator results_iter = results.begin();
	vector<PresearchedResult>::iterator results_end = results.end();
	vector<PresearchedResult>::iterator added_results_iter =
			added_results.begin();
	vector<PresearchedResult>::iterator added_results_end = added_results.end();
	while (new_results.size() < parameters.max_number_results
			&& (added_results_iter != added_results_end
					|| results_iter != results_end)) {
		PresearchedResult added_alignment;
		if (added_results_iter == added_results_end) {
			added_alignment = *results_iter;
			assert(
					added_alignment.database_chunk_id
							< database.GetNumberChunks());
			++results_iter;
		} else if (results_iter == results_end
				|| AlignmentComp(results_iter->score, results_iter->start,
						results_iter->end, added_results_iter->score,
						added_results_iter->start, added_results_iter->end)
						< 0) {
			added_alignment.database_chunk_id = database.GetChunkId();
			added_alignment.subject_id_in_chunk = database.GetId(
					added_results_iter->start.database_position);
			added_alignment.score = added_results_iter->score;
			added_alignment.start = added_results_iter->start;
			added_alignment.end = added_results_iter->end;
			added_alignment.hit = added_results_iter->hit;
			added_alignment.hit_count = added_results_iter->hit_count;
			assert(
					added_alignment.database_chunk_id
							< database.GetNumberChunks());
			++added_results_iter;
		} else {
			added_alignment = *results_iter;
			assert(
					added_alignment.database_chunk_id
							< database.GetNumberChunks());
			++results_iter;
		}
		if (added_alignment.database_chunk_id == database.GetChunkId()) {
			if (alignment_start_positions_in_database[added_alignment.subject_id_in_chunk].size()
					< parameters.max_number_one_subject_results) {
				vector<AlignmentPosition> *start_positions =
						&alignment_start_positions_in_database[added_alignment.subject_id_in_chunk];
				AlignmentPosition searched_position;
				searched_position.result_id = 0;
				searched_position.position =
						added_alignment.start.database_position;
				vector<AlignmentPosition>::iterator start_positions_iter =
						lower_bound(start_positions->begin(),
								start_positions->end(), searched_position,
								AlignmentPositionLessPosition());
				vector<AlignmentPosition>::iterator start_positions_insert_iter =
						start_positions_iter;
				searched_position.position =
						added_alignment.end.database_position;
				vector<AlignmentPosition>::iterator start_positions_iter_end =
						lower_bound(start_positions_iter,
								start_positions->end(), searched_position,
								AlignmentPositionLessPosition());

				bool added = true;
				for (; start_positions_iter != start_positions_iter_end;
						++start_positions_iter) {
					PresearchedResult *r_ptr =
							&new_results[start_positions_iter->result_id];
					if (IsContained(r_ptr->start, r_ptr->end,
							added_alignment.start, added_alignment.end)) {
						added = false;
						break;
					}
				}

				if (added == true) {
					AlignmentPosition p;
					p.position = added_alignment.start.query_position;
					p.result_id = new_results.size();
					alignment_start_positions_in_database[added_alignment.subject_id_in_chunk].insert(
							start_positions_insert_iter, p);
					new_results.push_back(added_alignment);
				}
			}
		} else {
			new_results.push_back(added_alignment);
		}
	}
	results.clear();
	results.insert(results.begin(), new_results.begin(), new_results.end());
	for (vector<PresearchedResult>::iterator iter = new_results.begin();
			iter != new_results.end(); ++iter) {
		//if (iter->database_chunk_id == database.GetChunkId()) {
		alignment_start_positions_in_database[iter->subject_id_in_chunk].clear();
		//}
	}
#if 0
// debug //////////////////////////////////
	cout << "results" << endl;
	for (uint32_t x = 0; x < results.size(); ++x) {
		cout << results[x].score << " "<<results[x].start.query_position << " "<< results[x].start.database_position << " "<< results[x].end.query_position << " " << results[x].end.database_position<<endl;
	}
//////////////////////////////////////////
#endif
}

void Aligner::BuildResult(Query &query, Database &database,
		uint32_t database_chunk_id, uint32_t subject_id_in_chunk, int score,
		Coordinate &hit, Coordinate &start, Coordinate &end,
		vector<EditBlocks::EditOpType> &edits, Result &result) {
	AlphabetCoder::Code *query_sequence = query.GetSequence();
	AlphabetCoder::Code *database_sequence = database.GetConcatenatedSequence();
	uint32_t query_position = start.query_position;
	uint32_t database_position = start.database_position;
	uint32_t m = 0;
	uint32_t o = 0;
	uint32_t e = 0;
	EditBlocks::EditOpType prev_op = EditBlocks::kSubstitution;
	for (vector<EditBlocks::EditOpType>::iterator it = edits.begin();
			it != edits.end(); ++it) {
		switch (*it) {
		case EditBlocks::kSubstitution:
			if (query_sequence[query_position]
					!= database_sequence[database_position]) {
				++m;
			}
			++query_position;
			++database_position;
			break;
		case EditBlocks::kGapInSeq0:
			if (prev_op != EditBlocks::kGapInSeq0) {
				++o;
			} else {
				++e;
			}
			++database_position;
			break;
		case EditBlocks::kGapInSeq1:
			if (prev_op != EditBlocks::kGapInSeq1) {
				++o;
			} else {
				++e;
			}
			++query_position;
			break;
		default:
			abort();
			break;
		}
		prev_op = *it;
	}
#if 0
/// debug //////////////////////////////////////////////////
	cout << result.score << " "<<result.start.query_position << " "<< result.start.database_position << " "<< result.end.query_position << " " << result.end.database_position<<endl;
/////////////////////////////////////////////////////////////
#endif
	result.database_chunk_id = database_chunk_id;
	result.subject_id_in_chunk = subject_id_in_chunk;
	result.subject_name = database.GetName(subject_id_in_chunk);
	result.score = score;
	result.alignment_length = edits.size();
	result.mismatches = m;
	result.gap_openings = o;
	result.gap_extensions = e;
	result.hit = hit;
	result.start.query_position = query.GetRealStart(start.query_position);
	result.end.query_position = query.GetRealEnd(end.query_position);
	result.start.database_position = start.database_position
			- (database.GetOffset(subject_id_in_chunk) - 1);
	result.end.database_position = end.database_position
			- (database.GetOffset(subject_id_in_chunk) - 1);
}
