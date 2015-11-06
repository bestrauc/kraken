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
#include "filemanager.hpp"
#include <boost/filesystem.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>



using namespace boost::filesystem;

namespace kraken {

typedef struct{
	std::vector<uint32_t> taxa;
	std::vector<uint8_t> ambig_list;
	std::unordered_map<uint32_t, uint32_t> hit_counts;

	/*uint64_t current_bin_key;
	int64_t current_min_pos = 1;
	int64_t current_max_pos = 0;*/

	uint32_t hits = 0;  // only maintained if in quick mode
	uint32_t taxon = 0;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version){
	    ar & taxa;
	    ar & ambig_list;
	    ar & hit_counts;
	    ar & hits;
	    ar & taxon;
	}
} SeqClassifyInfo;

typedef std::vector<std::shared_ptr<SeqClassifyInfo> > TRunInfoList;

struct RunInfoContainer{
	TRunInfoList runInfoList;

	RunInfoContainer(int count, int lane, int tile) : count(count), lane_num(lane), tile_num(tile) {
		processing_lock.lock(); // lock writing until all sequences have been read
	};

	std::atomic_uint count;
	std::mutex  processing_lock;

	int lane_num;
	int tile_num;

	void increment_count(){
		++count; // atomic increase

		// Release run information for writing when all information has been updated.
		if (count == runInfoList.size())
			processing_lock.unlock();
	}

};

struct DNASequence{
	std::string id;
	std::string header_line;  // id + optional description
	std::string seq;
	std::string quals;

	// only used for BCL reader, otherwise null
	std::shared_ptr<SeqClassifyInfo> readInfo;
	std::shared_ptr<RunInfoContainer> runContainer;

	/*DNASequence() = default;
	DNASequence(DNASequence&&) = default;
	DNASequence& operator=(DNASequence&&) = default;
	DNASequence(const DNASequence&) = delete;
	DNASequence& operator=(const DNASequence&) = delete;
    ~DNASequence() = default;*/
};


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
	typedef std::unordered_map<int, path> TTilePathMap;

	BCLReader(std::string filename);
	BCLReader(std::string filename, int length);
	//BCLReader(std::string filename, int length, std::vector<SeqClassifyInfo> *runInfoList);
	DNASequence next_sequence();
	bool is_valid();

private:
	// The path data structures for the lane, cycles and temporary files.
	path basecalls_path;
	std::vector<path> lanePaths;
	std::vector<std::vector<path> > cyclePaths;
	std::unordered_map<int, TTilePathMap> tmpPaths;

	BCLFileManager fileManager;


	// data structures indicating read status
	std::unordered_map<int, std::vector<bool> > readsFinished;
	std::vector<bool> tileFinished;

	// A list of progress for the reads of one tile
	//std::vector<std::unique_ptr<SeqClassifyInfo> > runInfoList;
	std::unique_ptr<TDNABuffer> sequenceBuffer;

	// Threading data
	// ------------------------
	std::unordered_map<int, std::unordered_map<int, std::mutex> > writeLocks;

	// Concurrent queues for the threads
	Queue<std::unique_ptr<TDNABuffer> > concurrentBufferQueue;
	Queue<std::shared_ptr<RunInfoContainer> > concurrentRunInfoQueue;

	// Information about the state of the reader
	//int tile_num;
	//int lane_num;
	int read_length;
	bool valid,_valid;

	// Private function definitions
	bool fillSequenceBuffer();
	void addSequenceBuffer(TileInfo tile);
	void saveRunInfo();
	void readInfo(std::shared_ptr<RunInfoContainer> runInfoContainer);
	void writeInfo(std::shared_ptr<RunInfoContainer> runInfoContainer);

	// Utility functions
	void getFilePaths();
	void init(std::string filename);
};
}

#endif
