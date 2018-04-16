/*
 * database_chunk.h
 *
 *  Created on: 2012/11/08
 *      Author: shu
 */

#ifndef DATABASE_CHUNK_H_
#define DATABASE_CHUNK_H_

#include <fstream>
#include <string>
#include <stdint.h>
#include <tr1/memory>
#include "alphabet_coder.h"
#include "seed_searcher.h"

class DatabaseChunk {
public:
	DatabaseChunk();
	DatabaseChunk(std::vector<Sequence> &sequences, AlphabetCoder &coder,
			AlphabetCoder::Code sequence_delimiter,
			SeedSearcherDatabaseParameters::BuildParameters &seed_search_parameters_build_parameters);
	DatabaseChunk(std::string filename_prefix);

	bool Build(std::vector<Sequence> &sequences, AlphabetCoder &coder,
			AlphabetCoder::Code sequence_delimiter,
			SeedSearcherDatabaseParameters::BuildParameters &seed_search_parameters_build_parameters);
	bool Clear();
	uint32_t GetNumberSequences() {
		return number_sequences_;
	}

	uint32_t GetConcatenatedSequenceLength() {
		return concatenated_sequences_length_;
	}

	std::string GetName(uint32_t id) {
#pragma omp critical(load_lock_for_names)
		{
			if (names_.size() == 0) {
				LoadNames(filename_prefix_);
			}
		}
		return names_[id];
	}

	uint32_t GetOffsets(uint32_t id) {
#pragma omp critical(load_lock_for_offsets)
		{
			if (offsets_.size() == 0) {
				LoadOffsets(filename_prefix_);
			}
		}
		return offsets_[id];
	}

	AlphabetCoder::Code *GetConcatenatedSequence() {
#pragma omp critical(load_lock_for_concatenated_sequence)
		{
			if (concatenated_sequence_.size() == 0) {
				LoadConcatenatedSequence(filename_prefix_);
			}
		}
		return &concatenated_sequence_[0];
	}

	SeedSearcherDatabaseParameters &GetSeedSearcherParameters() {
#pragma omp critical(load_lock_for_seed_searcher_parameters)
		{
			if (!setted_seed_searcher_parameters_) {
				LoadSeedSearcherParameters(filename_prefix_);
			}
		}
		return seed_searcher_parameters_;
	}

	uint32_t GetId(uint32_t position);

	bool Load(std::string filename_prefix);
	bool Save(std::string filename_prefix);
private:

	static std::string GetInformationFileName(std::string filename_prefix) {
		return filename_prefix + ".inf";
	}
	static std::string GetOffsetsFileName(std::string filename_prefix) {
		return filename_prefix + ".off";
	}
	static std::string GetNamesFileName(std::string filename_prefix) {
		return filename_prefix + ".nam";
	}
	static std::string GetSequencesFileName(std::string filename_prefix) {
		return filename_prefix + ".seq";
	}
	static std::string GetIndexFileName(std::string filename_prefix) {
		return filename_prefix + ".src";
	}

	void EncodeSequences(AlphabetCoder &coder, std::vector<Sequence> &sequences,
			AlphabetCoder::Code sequence_delimiter);

	bool LoadInfomation(std::string filename_prefix);
	bool LoadOffsets(std::string filename_prefix);
	bool LoadNames(std::string filename_prefix);
	bool LoadConcatenatedSequence(std::string filename_prefix);
	bool LoadSeedSearcherParameters(std::string filename_prefix);

	bool SaveInfomation(std::string filename_prefix);
	bool SaveOffsets(std::string filename_prefix);
	bool SaveNames(std::string filename_prefix);
	bool SaveConcatenatedSequence(std::string filename_prefix);
	bool SaveSeedSearcherParameters(std::string filename_prefix);

	bool building_;
	std::string filename_prefix_;
	uint32_t number_sequences_;
	uint32_t concatenated_sequences_length_;
	std::vector<std::string> names_;
	std::vector<uint32_t> offsets_;
	std::vector<AlphabetCoder::Code> concatenated_sequence_;
	bool setted_seed_searcher_parameters_;
	SeedSearcherDatabaseParameters seed_searcher_parameters_;
};

#endif /* DATABASE_CHUNK_H_ */
