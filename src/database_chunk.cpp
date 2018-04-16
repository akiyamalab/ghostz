/*
 * database_chunk.cpp
 *
 *  Created on: 2012/11/08
 *      Author: shu
 */

#include <fstream>
#include <stdint.h>
#include <string>
#include <assert.h>
#include <limits.h>
#include "sequence.h"
#include "alphabet_coder.h"
#include "seed_searcher.h"
#include "database_chunk.h"
using namespace std;

DatabaseChunk::DatabaseChunk() :
    building_(false), filename_prefix_(""), number_sequences_(0), concatenated_sequences_length_(0), names_(
        0), offsets_(0), concatenated_sequence_(0), setted_seed_searcher_parameters_(false), seed_searcher_parameters_() {
}

DatabaseChunk::DatabaseChunk(std::vector<Sequence> &sequences, AlphabetCoder &coder,
    AlphabetCoder::Code sequence_delimiter,
    SeedSearcherDatabaseParameters::BuildParameters &seed_search_parameters_build_parameters) :
    building_(false), filename_prefix_(""), number_sequences_(0), concatenated_sequences_length_(0), names_(
        0), offsets_(0), concatenated_sequence_(0), setted_seed_searcher_parameters_(false), seed_searcher_parameters_()  {
  Build(sequences, coder, sequence_delimiter, seed_search_parameters_build_parameters);
}

DatabaseChunk::DatabaseChunk(string filename_prefix) :
    building_(false), filename_prefix_(""), number_sequences_(0), concatenated_sequences_length_(0), names_(
        0), offsets_(0), concatenated_sequence_(0), seed_searcher_parameters_()  {
  Load(filename_prefix);
}

bool DatabaseChunk::Build(std::vector<Sequence> &sequences, AlphabetCoder &coder,
    AlphabetCoder::Code sequence_delimiter,
    SeedSearcherDatabaseParameters::BuildParameters &seed_search_parameters_build_parameters) {
  number_sequences_ = sequences.size();
  names_.resize(sequences.size());
  offsets_.resize(sequences.size() + 1);
  concatenated_sequences_length_ = 1;
  for (unsigned int i = 0; i < sequences.size(); ++i) {
    concatenated_sequences_length_ += sequences[i].GetSequenceData().size() + 1;
    names_[i] = sequences[i].GetName();
  }
  EncodeSequences(coder, sequences, sequence_delimiter);
  seed_searcher_parameters_.Build(&concatenated_sequence_[0], concatenated_sequence_.size(),
      seed_search_parameters_build_parameters);
  setted_seed_searcher_parameters_ = true;
  building_ = true;
  return true;
}

bool DatabaseChunk::Clear() {
  number_sequences_ = 0;
  concatenated_sequences_length_ = 0;
  names_.clear();
  offsets_.clear();
  concatenated_sequence_.clear();
  setted_seed_searcher_parameters_ = false;
  seed_searcher_parameters_ = SeedSearcherDatabaseParameters();
  building_ = false;
  return true;
}

uint32_t DatabaseChunk::GetId(uint32_t position) {
#pragma omp critical(load_lock_for_GetId)
  {
    if (offsets_.size() == 0) {
      LoadOffsets(filename_prefix_);
    }
  }
  if (offsets_[number_sequences_ - 1] <= position && position < concatenated_sequences_length_) {
    return number_sequences_ - 1;
  }

  uint32_t left = 0;
  uint32_t right = number_sequences_ - 2;
  uint32_t mid;
  while (left <= right) {
    mid = (left + right) / 2;
    if (offsets_[mid] <= position && position < offsets_[mid + 1]) {
      return mid;
    } else if (offsets_[mid] < position) {
      left = mid + 1;
    } else {
      right = mid - 1;
    }
  }
  return UINT_MAX;
}

bool DatabaseChunk::Load(std::string filename_prefix) {
  filename_prefix_ = filename_prefix;
  names_.clear();
  offsets_.clear();
  concatenated_sequence_.clear();
  setted_seed_searcher_parameters_ = false;
  return LoadInfomation(filename_prefix);
}
bool DatabaseChunk::Save(std::string filename_prefix) {
  bool ret = false;
  if (building_ == true) {
    ret = true;
    ret &= SaveInfomation(filename_prefix);
    ret &= SaveNames(filename_prefix);
    ret &= SaveOffsets(filename_prefix);
    ret &= SaveConcatenatedSequence(filename_prefix);
    ret &= SaveSeedSearcherParameters(filename_prefix);
  }
  return ret;
}

void DatabaseChunk::EncodeSequences(AlphabetCoder &coder, vector<Sequence> &sequences,
    AlphabetCoder::Code sequence_delimiter) {
  concatenated_sequence_.resize(concatenated_sequences_length_);
  offsets_.resize(sequences.size() + 1);

  concatenated_sequence_[0] = sequence_delimiter;
  uint32_t offset = 1;
  for (unsigned int i = 0; i < sequences.size(); ++i) {
    offsets_[i] = offset;
    string sequence = sequences[i].GetSequenceData();
    coder.Encode(&sequence[0], sequence.size(), &concatenated_sequence_[offset]);
    offset += sequence.size();
    concatenated_sequence_[offset] = sequence_delimiter;
    ++offset;
  }
  offsets_[sequences.size()] = offset;
}

bool DatabaseChunk::LoadInfomation(std::string filename_prefix) {
  string filename = GetInformationFileName(filename_prefix);
  ifstream in;
  in.open(filename.c_str(), ios::binary);
  if (in) {
    in.read((char *) &number_sequences_, sizeof(number_sequences_));
    in.read((char *) &concatenated_sequences_length_, sizeof(concatenated_sequences_length_));

    in.close();
    return true;
  }
  return false;
}
bool DatabaseChunk::LoadOffsets(std::string filename_prefix) {
  string filename = GetOffsetsFileName(filename_prefix);
  ifstream in;
  in.open(filename.c_str(), ios::binary);
  if (in) {
    offsets_.resize(number_sequences_ + 1);
    in.read((char *) &offsets_[0], sizeof(offsets_[0]) * (number_sequences_ + 1));
    in.close();
    return true;
  }
  return false;
}

bool DatabaseChunk::LoadNames(std::string filename_prefix) {
  ifstream in;
  string line;
  string filename = GetNamesFileName(filename_prefix);
  in.open(filename.c_str());
  if (in) {
    names_.resize(number_sequences_);
    uint32_t i;
    for (i = 0; i < number_sequences_ && !in.eof(); ++i) {
      std::getline(in, line);
      names_[i] = line;
    }
    assert(i == number_sequences_);
    in.close();
    return true;
  }
  return false;
}
bool DatabaseChunk::LoadConcatenatedSequence(std::string filename_prefix) {
  string filename = GetSequencesFileName(filename_prefix);
  ifstream in;
  in.open(filename.c_str(), ios::binary);
  if (in) {
    concatenated_sequence_.resize(concatenated_sequences_length_);
    in.read((char *) &concatenated_sequence_[0],
        sizeof(concatenated_sequence_[0]) * concatenated_sequences_length_);
    in.close();
    return true;
  }
  return false;
}
bool DatabaseChunk::LoadSeedSearcherParameters(std::string filename_prefix) {
  string filename = GetIndexFileName(filename_prefix);
  ifstream in;
  in.open(filename.c_str(), ios::binary);
  if (in) {
    seed_searcher_parameters_.Build(GetConcatenatedSequence(), GetConcatenatedSequenceLength(), in);
    in.close();
    setted_seed_searcher_parameters_ = true;
    return true;
  }
  return false;
}

bool DatabaseChunk::SaveInfomation(std::string filename_prefix) {
  string filename = GetInformationFileName(filename_prefix);
  ofstream out;
  out.open(filename.c_str(), ios::binary);
  out.write((char *) &number_sequences_, sizeof(number_sequences_));
  out.write((char *) &concatenated_sequences_length_, sizeof(concatenated_sequences_length_));
  out.close();
  return true;
}

bool DatabaseChunk::SaveOffsets(std::string filename_prefix) {
  ofstream out;
  string filename = GetOffsetsFileName(filename_prefix);
  out.open(filename.c_str(), ios::binary);
  out.write((char *) &offsets_[0], sizeof(offsets_[0]) * offsets_.size());
  out.close();
  return true;
}
bool DatabaseChunk::SaveNames(std::string filename_prefix) {
  ofstream out;
  string filename = GetNamesFileName(filename_prefix);
  out.open(filename.c_str());
  for (vector<string>::iterator i = names_.begin(); i != names_.end(); ++i) {
    out << *i << endl;
  }
  out.close();
  return true;
}

bool DatabaseChunk::SaveConcatenatedSequence(std::string filename_prefix) {
  ofstream out;
  string filename = GetSequencesFileName(filename_prefix);
  out.open(filename.c_str(), ios::binary);
  out.write((char *) &concatenated_sequence_[0],
      sizeof(concatenated_sequence_[0]) * concatenated_sequence_.size());
  out.close();
  return true;
}

bool DatabaseChunk::SaveSeedSearcherParameters(std::string filename_prefix) {
  ofstream out;
  string filename = GetIndexFileName(filename_prefix);
  out.open(filename.c_str(), ios::binary);
  seed_searcher_parameters_.Save(out);
  out.close();
  return true;
}
