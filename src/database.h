/*
 * database.h
 *
 *  Created on: 2012/10/31
 *      Author: shu
 */

#ifndef DATABASE_H_
#define DATABASE_H_

#include <vector>
#include <string>
#include <stdint.h>
#include <tr1/memory>
#include "alphabet_coder.h"
#include "alphabet_type.h"
#include "fasta_sequence_reader.h"
#include "database_chunk.h"
#include "sequence.h"
#include "seed_searcher.h"

class Database {
public:
	typedef int ChunkBuildOption;
	static const ChunkBuildOption ChunkBuildOptionNotBuild = -2;
	static const ChunkBuildOption ChunkBuildOptionAllBuild = -1;

	typedef struct {
		unsigned int chunk_size;
		ChunkBuildOption chunk_build_option;
		std::tr1::shared_ptr<AlphabetType> sequence_type_ptr;
		AlphabetCoder::Code sequence_delimiter;
		SeedSearcherDatabaseParameters::BuildParameters seed_search_parameters_build_parameters;
	} Parameters;

	Database();
	Database(std::string filename_prefix);
	Database(std::istream &is, Parameters &parameters);

	bool Load(std::string filename_prefix);
	bool Save(std::string filename_prefix);

	bool SetChunk(int id);
	void ResetChunk();
	bool NextChunk();

	uint32_t GetNumberChunks() {
		if (!setting_informations_) {
			SetInfomations();
		}
		return number_chunks_;
	}

	int GetChunkId() {
		return chunk_id_;
	}

	uint32_t GetNumberSequencesInChunk() {
		return chunk_.GetNumberSequences();
	}

	uint32_t GetConcatenatedSequenceLength() {
		return chunk_.GetConcatenatedSequenceLength();
	}

	std::string GetName(uint32_t id) {
		return chunk_.GetName(id);
	}

	uint32_t GetOffset(uint32_t id) {
		return chunk_.GetOffsets(id);
	}

	AlphabetCoder::Code *GetConcatenatedSequence() {
		return chunk_.GetConcatenatedSequence();
	}

	uint32_t GetId(uint32_t position) {
		return chunk_.GetId(position);
	}

	SeedSearcherDatabaseParameters &GetSeedSearcherParameters() {
		return chunk_.GetSeedSearcherParameters();
	}

	uint64_t GetDatabaseTotalLenght() {
		if (!setting_informations_) {
			SetInfomations();
		}
		return database_length_;
	}

	uint64_t GetNumberTotalSequences() {
		if (!setting_informations_) {
			SetInfomations();
		}
		return number_sequences_;
	}

	uint32_t GetMaxSequenceLength() {
		if (!setting_informations_) {
			SetInfomations();
		}
		return max_sequence_length_;
	}

	AlphabetCoder::Code GetSequenceDelimiter() {
		return sequence_delimiter_;
	}

private:
	static std::string GetInformationFileName(std::string filename_prefix) {
		return filename_prefix + ".inf";
	}

	unsigned int ReadSequences(std::vector<Sequence> &sequences);
	bool SetInfomations();
	bool BuildNextChunk();
	bool LoadInfomation(std::string filename_prefix);
	bool SaveInformation(std::string filename_prefix);

	bool saving_;
	bool setting_informations_;
	Parameters parameters_;
	FastaSequenceReader reader_;
	std::string filename_prefix_;
	int number_chunks_;
	uint32_t max_sequence_length_;
	uint64_t database_length_;
	uint64_t number_sequences_;
	AlphabetCoder::Code sequence_delimiter_;
	int chunk_id_;
	DatabaseChunk chunk_;
	Sequence next_sequence_;

};

#endif /* DATABASE_H_ */
