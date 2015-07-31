/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SEQREADER_HPP
#define SEQREADER_HPP

#include "kraken_headers.hpp"
#include "Queue.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <memory>
#include <thread>
#include <algorithm>

using namespace boost::filesystem;

namespace kraken {

typedef struct{
	std::vector<uint32_t> taxa;
	std::vector<uint8_t> ambig_list;
	std::unordered_map<uint32_t, uint32_t> hit_counts;

	uint64_t current_bin_key;
	int64_t current_min_pos = 1;
	int64_t current_max_pos = 0;

	uint32_t hits = 0;  // only maintained if in quick mode
	uint32_t taxon = 0;
} SeqClassifyInfo;

typedef struct {
	std::string id;
	std::string header_line;  // id + optional description
	std::string seq;
	std::string quals;

	// only used for BCL reader, otherwise null
	SeqClassifyInfo* readInfo;
} DNASequence;


class DNASequenceReader {
public:
	virtual DNASequence next_sequence() = 0;
	virtual bool is_valid() = 0;
	virtual ~DNASequenceReader() {}
};

class FastaReader : public DNASequenceReader {
public:
	FastaReader(std::string filename);
	DNASequence next_sequence();
	bool is_valid();

private:
	std::ifstream file;
	std::string linebuffer;
	bool valid;
};

class FastqReader : public DNASequenceReader {
public:
	FastqReader(std::string filename);
	DNASequence next_sequence();
	bool is_valid();

private:
	std::ifstream file;
	bool valid;
};

class BCLReader : public DNASequenceReader {
public:
	typedef std::vector<DNASequence> TDNABuffer;

	BCLReader(std::string filename);
	BCLReader(std::string filename, int length);
	//BCLReader(std::string filename, int length, std::vector<SeqClassifyInfo> *runInfoList);
	DNASequence next_sequence();
	bool is_valid();

private:
	std::vector<path> lanePaths;
	std::vector<std::vector<path> > cyclePaths;

	// a map of (tile number, bit) to show if to process a read
	std::unordered_map<int, std::vector<bool> > readsFinished;
	// a vector of flags to indicate if to process a tile file
	std::vector<bool> tileFinished;

	// A list of progress for the reads of one tile
	//std::vector<std::unique_ptr<SeqClassifyInfo> > runInfoList;
	std::vector<SeqClassifyInfo*> runInfoList;

	Queue<std::unique_ptr<TDNABuffer> > concurrentBufferQueue;
	std::unique_ptr<TDNABuffer> sequenceBuffer;
	path basecalls_path;

	int tile_num;
	int lane_num;
	bool valid, _valid;
	int read_length;

	bool fillSequenceBuffer();
	void addSequenceBuffer(int lane_num, int tile_num);
	void getCyclePaths();
	void init(std::string filename);
};
}

#endif
