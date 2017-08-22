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

//#include <boost/serialization/vector.hpp>
//#include <boost/archive/binary_oarchive.hpp>
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/serialization/shared_ptr.hpp>


using namespace boost::filesystem;

namespace kraken {

struct SeqClassifyInfo{
	std::vector<uint32_t> taxa;
	std::vector<uint8_t> ambig_list;
	std::unordered_map<uint32_t, uint32_t> hit_counts;

	uint16_t pos = 0;
	uint32_t hits = 0;  // only maintained if in quick mode
	uint32_t last_ambig = 0;
	uint64_t last_kmer = 0;
	uint16_t processed_len = 0;
	bool skip = false;
	bool first = true;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version){
	    ar & taxa;
	    ar & ambig_list;
	    ar & hit_counts;
	    ar & hits;
	}
};

typedef std::vector<std::shared_ptr<SeqClassifyInfo> > TSeqClassifyInfoList;

struct RunInfoContainer{
	//TSeqClassifyInfoList runInfoList;
	uint64_t runsize;

	RunInfoContainer() :
			count(0), lane_num(0), tile_num(0), processed_nt(0) {	};

	RunInfoContainer(int lane, int tile, int processed) :
		count(0), lane_num(lane), tile_num(tile), processed_nt(processed) {
		//processing_lock.lock(); // lock writing until all sequences have been read
	};

	std::atomic<std::uint64_t> count;
	std::mutex  processing_lock;

	int lane_num;
	int tile_num;
	int processed_nt = 0;


	void increment_count(uint64_t inc){
		count += inc; // atomic increase

		uint64_t runsize_compare = runsize;
		// Release run information for writing when all information has been updated.
		// use atomic load and swap to not double unlock
		if (count.compare_exchange_strong(runsize_compare, 0)){
			processing_lock.unlock();
		}
	}
};

struct DNASequence{
	std::string id;
	std::string header_line;  // id + optional description
	std::string seq;
	std::string quals;

	// only used for BCL reader, otherwise null
	std::shared_ptr<SeqClassifyInfo> readInfo = std::make_shared<SeqClassifyInfo>();
};

typedef struct {
	std::shared_ptr<RunInfoContainer> runContainer;
    std::vector<DNASequence> seqs;
} WorkUnit;


class DNASequenceReader {
public:
	virtual DNASequence next_sequence() = 0;
	virtual bool is_valid() = 0;

	// default implementation of the workunit generation using only next_sequence
	virtual size_t next_workunit(size_t work_nt_size, WorkUnit& work_unit) {
		size_t total_nt = 0;
		while (total_nt < work_nt_size){
			DNASequence dna = this->next_sequence();
			if (!this->is_valid()){
				break;
			}

			work_unit.seqs.push_back(dna);
			total_nt += dna.seq.size();
		}

		return total_nt;
	}

	virtual ~DNASequenceReader() {}
};

class FastaReader : public DNASequenceReader {
public:
	FastaReader(std::string filename);
	FastaReader(std::string filename, int length);
	DNASequence next_sequence();
	bool is_valid();

private:
	std::ifstream file;
	std::string linebuffer;
	bool valid;
	int classify_length;
};

class FastqReader : public DNASequenceReader {
public:
	FastqReader(std::string filename);
	FastqReader(std::string filename, int length);
	DNASequence next_sequence();
	bool is_valid();

private:
	std::ifstream file;
	bool valid;
	int classify_length;
};

class BCLReader : public DNASequenceReader {
public:
//	typedef std::vector<DNASequence> WorkUnit;
	typedef std::unordered_map<int, path> TTilePathMap;

	//BCLReader(std::string filename);
	BCLReader(std::string filename, int length, int start_cycle, int step, std::vector<int> target_tiles,
						std::vector<int> target_lanes);
	DNASequence next_sequence();
    virtual size_t next_workunit(size_t work_nt_size, WorkUnit& work_unit);
	bool is_valid();


private:
	// The path data structures for the lane, cycles and temporary files.
	path basecalls_path;
	std::vector<path> lanePaths;
	std::vector<std::vector<path> > cyclePaths;
	std::unordered_map<int, TTilePathMap> tmpPaths;

	BCLFileManager fileManager;

	// data structures indicating read status
	std::vector<bool> tileFinished;
	std::unordered_map<int, std::unordered_map<int, std::shared_ptr<RunInfoContainer> > > runInfoMap;
	std::unordered_map<int, std::unordered_map<int, TSeqClassifyInfoList> > runMap;
	//std::unordered_map<int, std::unordered_map<int, std::mutex> > runInfoLocks;

	// A list of progress for the reads of one tile
	std::unique_ptr<WorkUnit> sequenceBuffer;

	// Threading data
	// ------------------------
	std::unordered_map<int, std::unordered_map<int, std::mutex> > writeLocks;

	// Concurrent queues for the threads
	Queue<std::unique_ptr<WorkUnit> > concurrentBufferQueue;
	Queue<std::shared_ptr<RunInfoContainer> > concurrentRunInfoQueue;

	// Information about the state of the reader
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
