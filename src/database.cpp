/*
 * database.cpp
 *
 *  Created on: 2012/10/31
 *      Author: shu
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <string>
#include "database_chunk.h"
#include "alphabet_coder.h"
#include "reduced_alphabet_coder.h"
#include "fasta_sequence_reader.h"
#include "sequence.h"
#include "database.h"

using namespace std;

Database::Database(std::string filename_prefix) :
    saving_(false), setting_informations_(false), filename_prefix_(filename_prefix), number_chunks_(
        0),  max_sequence_length_(0), database_length_(0), number_sequences_(0), sequence_delimiter_(
        0), chunk_id_(-1), chunk_(), next_sequence_("", "") {
  Load(filename_prefix);
}
Database::Database(std::istream &is, Parameters &parameters) :
    saving_(false), setting_informations_(false), parameters_(parameters), reader_(is), filename_prefix_(
        ""), number_chunks_(0),max_sequence_length_(0), database_length_(0), number_sequences_(
        0), sequence_delimiter_(parameters.sequence_delimiter), chunk_id_(-1), chunk_(), next_sequence_(
        "", "") {
  BuildNextChunk();
}

bool Database::Load(std::string filename_prefix) {
  saving_ = true;
  setting_informations_ = true;
  LoadInfomation(filename_prefix);
  SetChunk(0);
  return true;
}

bool Database::Save(std::string filename_prefix) {
  filename_prefix_ = filename_prefix;
  int current_id = GetChunkId();
  ResetChunk();
  do {
    stringstream ss;
    ss << filename_prefix_ << "_" << GetChunkId();
    chunk_.Save(ss.str());
  } while (BuildNextChunk());
  SaveInformation(filename_prefix_);
  saving_ = true;
  setting_informations_ = true;
  SetChunk(current_id);
  return true;
}

bool Database::SetChunk(int id) {
  if (id == chunk_id_) {
    return true;
  }
  if (!saving_) {
    reader_.Seek(0);
    for (int i = 0; i < id; ++i) {
      BuildNextChunk();
    }
  } else {
    if (id >= number_chunks_) {
      return false;
    }

    stringstream prefix;
    prefix.str("");
    prefix << filename_prefix_ << "_" << id;
    chunk_.Load(prefix.str());
    chunk_id_ = id;
    return true;
  }
  return false;
}

void Database::ResetChunk() {
  SetChunk(0);
}

bool Database::NextChunk() {
  if (!saving_) {
    return BuildNextChunk();
  } else if (chunk_id_ + 1 < number_chunks_) {
    SetChunk(chunk_id_ + 1);
    return true;
  } else {
    chunk_id_ = number_chunks_;
    return false;
  }
}

unsigned int Database::ReadSequences(std::vector<Sequence> &sequences) {
  unsigned int sum_length = 0;
  string name;
  string sequence;
  sequences.clear();
  bool reader_ret = true;
  if (next_sequence_.GetName().size() == 0 && next_sequence_.GetSequenceData().size() == 0) {
    reader_ret = reader_.Read(name, sequence);
    next_sequence_ = Sequence(name, sequence);
  }

  while (reader_ret) {
    unsigned int new_sum_length = sum_length + next_sequence_.GetSequenceData().size();
    if (new_sum_length >= parameters_.chunk_size) {
      break;
    }
    sequences.push_back(next_sequence_);
    sum_length = new_sum_length;
    reader_ret = reader_.Read(name, sequence);
    if (reader_ret) {
      next_sequence_ = Sequence(name, sequence);
    } else {
      next_sequence_ = Sequence("", "");
      break;
    }
  }
  return sum_length;
}

bool Database::SetInfomations() {
  if (setting_informations_) {
    return true;
  } else if (reader_.IsRead()) {
    uint32_t id = GetChunkId();
    while (BuildNextChunk()) {
    }
    setting_informations_ = true;
    SetChunk(id);
    return true;
  } else {
    return false;
  }
}
bool Database::BuildNextChunk() {
  vector<Sequence> sequences;
  AlphabetCoder coder(*(parameters_.sequence_type_ptr));

  vector<AlphabetCoder::Code> encoded_sequences;
  vector<unsigned int> encoded_sequence_offsets;
  database_length_ += ReadSequences(sequences);
  if (sequences.empty()) {
    number_chunks_ = chunk_id_ + 1;
    return false;
  }
  for (size_t i = 0; i < sequences.size(); ++i) {
    max_sequence_length_ = max(max_sequence_length_,
        (uint32_t) sequences[i].GetSequenceData().size());
  }
  number_sequences_ += sequences.size();
  if (parameters_.chunk_build_option == Database::ChunkBuildOptionAllBuild
      || parameters_.chunk_build_option == (chunk_id_ + 1)) {
    chunk_.Build(sequences, coder, parameters_.sequence_delimiter,
        parameters_.seed_search_parameters_build_parameters);
  } else {
    chunk_.Clear();
  }
  ++chunk_id_;
  return true;
}
bool Database::LoadInfomation(std::string filename_prefix) {
  filename_prefix_ = filename_prefix;
  ifstream in;
  string filename = GetInformationFileName(filename_prefix_);
  in.open(filename.c_str(), ios::binary);
  if (in) {
    in.read((char *) &number_chunks_, sizeof(number_chunks_));
    in.read((char *) &max_sequence_length_, sizeof(max_sequence_length_));
    in.read((char *) &database_length_, sizeof(database_length_));
    in.read((char *) &number_sequences_, sizeof(number_sequences_));
    in.read((char *) &sequence_delimiter_, sizeof(sequence_delimiter_));
    in.close();
    return true;
  }
  return false;
}
bool Database::SaveInformation(std::string filename_prefix) {
  filename_prefix_ = filename_prefix;
  string filename = GetInformationFileName(filename_prefix);
  ofstream out;
  out.open(filename.c_str(), ios::binary);
  out.write((char *) &number_chunks_, sizeof(number_chunks_));
  out.write((char *) &max_sequence_length_, sizeof(max_sequence_length_));
  out.write((char *) &database_length_, sizeof(database_length_));
  out.write((char *) &number_sequences_, sizeof(number_sequences_));
  out.write((char *) &sequence_delimiter_, sizeof(sequence_delimiter_));
  out.close();
  return true;
}
