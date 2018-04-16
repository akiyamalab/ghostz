/*
 * aligner.h
 *
 *  Created on: 2012/10/05
 *      Author: shu
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <string>
#include <vector>
#include <tr1/memory>
#include "alphabet_coder.h"
#include "edit_blocks.h"
#include "alphabet_type.h"
#include "score_matrix.h"
#include "queries.h"
#include "database.h"

class Aligner {
public:
	typedef struct {
		bool clustering;
		uint32_t clustering_subsequence_length;
		uint32_t number_threads;
		std::string filename_prefix;
		unsigned int chunk_size;
		Database::ChunkBuildOption chunk_build_option;
		std::tr1::shared_ptr<AlphabetType> sequence_type_ptr;
		unsigned int seed_threshold;
		ScoreMatrix score_matrix;
		std::vector<std::string> hash_alphabet_sets;
		std::vector<std::string> similarity_alphabet_sets;
	} DatabaseParameters;

	typedef struct {
		bool filter;
		int gap_open;
		int gap_extension;
		float normalized_presearched_ungapped_extension_cutoff;
		float normalized_presearched_gapped_extension_trigger;
		float normalized_presearched_gapped_extension_cutoff;
		float normalized_result_gapped_extension_cutoff;
		uint32_t number_threads;
		unsigned int queries_chunk_size;
		unsigned int max_number_results;
		unsigned int max_number_one_subject_results;
		ScoreMatrix score_matrix;
		std::tr1::shared_ptr<AlphabetType> queries_file_sequence_type_ptr;
		std::tr1::shared_ptr<AlphabetType> aligning_sequence_type_ptr;
	} AligningParameters;

	void BuildDatabase(std::string &input_filename,
			std::string &database_filename, DatabaseParameters paramters);
	void Align(std::string &queries_filename, std::string &database_filename,
			std::string &output_filename, AligningParameters &parameters);

private:
	typedef struct {
		uint32_t query_position;
		uint32_t database_position;
	} Coordinate;

	typedef struct {
		uint32_t database_chunk_id;
		uint32_t subject_id_in_chunk;
		int score;
		Coordinate hit;
		Coordinate start;
		Coordinate end;
		uint32_t hit_count;
	} PresearchedResult;

	typedef struct {
		uint32_t result_id;
		uint32_t position;
	} AlignmentPosition;

	typedef struct {
		uint32_t database_chunk_id;
		uint32_t subject_id_in_chunk;
		std::string subject_name;
		int score;
		uint32_t alignment_length;
		uint32_t mismatches;
		uint32_t gap_openings;
		uint32_t gap_extensions;
		Coordinate hit;
		Coordinate start;
		Coordinate end;
		uint32_t hit_count;
	} Result;

	static int AlignmentComp(const int a_score, const Coordinate &a_start,
			const Coordinate &a_end, const int b_score,
			const Coordinate &b_start, const Coordinate &b_end);

	static bool IsContained(Coordinate &a_start, Coordinate &a_end,
			Coordinate &b_start, Coordinate &b_end);

	struct AlignmentPositionLessPosition {
		bool operator()(const AlignmentPosition &a1,
				const AlignmentPosition &a2) const {
			return a1.position < a2.position;
		}
	};

	struct PresearchedResultGreaterScore {
		bool operator ()(const PresearchedResult &r1,
				const PresearchedResult &r2) {
			return AlignmentComp(r1.score, r1.start, r1.end, r2.score, r2.start,
					r2.end) > 0;
		}
	};

	struct ResultGreaterScore {
		bool operator ()(const Result &r1, const Result &r2) {
			return AlignmentComp(r1.score, r1.start, r1.end, r2.score, r2.start,
					r2.end) > 0;
		}
	};

	AlphabetCoder::Code GetSequenceDelimiter(AlphabetType &type);
	void BuildQueriesParameters(AligningParameters &parameters,
			Queries::Parameters &queries_parameters);
	void BuildDatabaseParameters(DatabaseParameters &parameters,
			Database::Parameters &database_parameters);
	void Presearch(Queries &queries, Database &database,
			AligningParameters &parameters,
			std::vector<std::vector<PresearchedResult> > &results_list);
	void BuildResults(Queries &queries, Database &database,
			AligningParameters &parameters,
			std::vector<std::vector<PresearchedResult> > &presearch_results_list,
			std::vector<std::vector<Result> > &results_list);
	void WriteOutput(std::ostream &os, Queries &queries, Database &database,
			AligningParameters &parameters,
			std::vector<std::vector<Result> > &results);
	void AddResults(Database &database, AligningParameters &parameters,
			std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database,
			std::vector<PresearchedResult> &added_results,
			std::vector<PresearchedResult> &results);
	void BuildResult(Query &query, Database &database,
			uint32_t database_chunk_id, uint32_t subject_id_in_chunk, int score,
			Coordinate &hit, Coordinate &start, Coordinate &end,
			std::vector<EditBlocks::EditOpType> &edits, Result &result);
};
#endif /* ALIGNER_H_ */
