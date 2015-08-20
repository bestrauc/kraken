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

#include "kraken_headers.hpp"
#include "seqreader.hpp"

// for timing and debug
#include <sys/time.h>
#include <ctime>

#ifdef DEBUG
#   define LOG(x) cerr << x
#else
#   define LOG(x) do {} while (0)
#endif

uint64_t GetTimeMs64(){
	struct timeval tv;

	gettimeofday(&tv, NULL);

	uint64_t ret = tv.tv_usec;
	/* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
	ret /= 1000;

	/* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
	ret += (tv.tv_sec * 1000);

	return ret;
}

using namespace std;
namespace fs = boost::filesystem;

namespace kraken {
FastaReader::FastaReader(string filename) {
	file.open(filename.c_str());
	if (file.rdstate() & ifstream::failbit) {
		err(EX_NOINPUT, "can't open %s", filename.c_str());
	}
	valid = true;
}

DNASequence FastaReader::next_sequence() {
	DNASequence dna;

	if (! valid || ! file.good()) {
		valid = false;
		return dna;
	}
	string line;

	if (linebuffer.empty()) {
		getline(file, line);
	}
	else {
		line = linebuffer;
		linebuffer.clear();
	}

	if (line[0] != '>') {
		warnx("malformed fasta file - expected header char > not found");
		valid = false;
		return dna;
	}
	dna.header_line = line.substr(1);
	istringstream seq_id(dna.header_line);
	seq_id >> dna.id;

	ostringstream seq_ss;

	while (file.good()) {
		getline(file, line);
		if (line[0] == '>') {
			linebuffer = line;
			break;
		}
		else {
			seq_ss << line;
		}
	}
	dna.seq = seq_ss.str();

	if (dna.seq.empty()) {
		warnx("malformed fasta file - zero-length record (%s)", dna.id.c_str());
		valid = false;
		return dna;
	}

	return dna;
}

bool FastaReader::is_valid() {
	return valid;
}

FastqReader::FastqReader(string filename) {
	file.open(filename.c_str());
	if (file.rdstate() & ifstream::failbit) {
		err(EX_NOINPUT, "can't open %s", filename.c_str());
	}
	valid = true;
}

DNASequence FastqReader::next_sequence() {
	DNASequence dna;

	if (! valid || ! file.good()) {
		valid = false;
		return dna;
	}

	string line;
	getline(file, line);
	if (line.empty()) {
		valid = false;  // Sometimes FASTQ files have empty last lines
		return dna;
	}
	if (line[0] != '@') {
		if (line[0] != '\r')
			warnx("malformed fastq file - sequence header (%s)", line.c_str());
		valid = false;
		return dna;
	}
	dna.header_line = line.substr(1);
	istringstream line_ss(dna.header_line);

	line_ss >> dna.id;
	getline(file, dna.seq);

	getline(file, line);
	if (line.empty() || line[0] != '+') {
		if (line[0] != '\r')
			warnx("malformed fastq file - quality header (%s)", line.c_str());
		valid = false;
		return dna;
	}
	getline(file, dna.quals);

	return dna;
}

bool FastqReader::is_valid() {
	return valid;
}

//-----------------------------------------------------------
//              BCL READER helper functions
//-----------------------------------------------------------


// BCLReader helper functions
char numToDNA(int base){
	switch (base)
	{
	case 0:
		return 'A';
	case 1:
		return 'C';
	case 2:
		return 'G';
	case 3:
		return 'T';
	default:
		return 'N';
	}
}

bool scanFilter(const fs::path &filter_path, std::vector<bool> &tile_index){
	std::ifstream in_file;
	in_file.open(filter_path.c_str(), std::ios::in | std::ios::binary);
	in_file.seekg(8);

	uint32_t N;
	in_file.read((char*)&N, 4);

	tile_index.resize(N);

	tile_index.assign(std::istreambuf_iterator<char>(in_file), std::istreambuf_iterator<char>());

	return true;
}

// Scan the tile given in the tile_path and save sequences into the buffer
bool scanTile(int tile_num, const fs::path &tile_path, std::vector<bool> &tile_filter,
		std::shared_ptr<RunInfoContainer> &runInfo, std::unique_ptr<BCLReader::TDNABuffer> &buffer){

	std::ifstream in_file;
	in_file.open(tile_path.c_str(), std::ios::in | std::ios::binary);

	if (!in_file.is_open()){
		std::cout << "Unable to open file: " << tile_path << "\n";
		return false;
	}

	// Process file
	uint32_t N; // number of clusters (sequences)
	uint32_t pos = 0;
	//static char raw_buffer[8192]; //or 4096?
	static char raw_buffer[16384]; //or 4096?

	in_file.read((char*)&N, 4);

	// This should hold if all files are correct and were correctly read
	assert(N == tile_filter.size());

	buffer->resize(N);

	// If runInfo is not yet initialized, resize it and insert .
	if (runInfo->runInfoList.size() != N){
		runInfo->runInfoList.resize(N);
//		std::generate_n(std::back_inserter(*runInfo), N, std::make_shared<SeqClassifyInfo>);
	}

	// Variables for timing measurement. (On one line to comment it out easily.)
	uint64_t t1=0, t2=0, read_time=0, process_time=0;

	while(!in_file.eof()){
		t1 = GetTimeMs64();

		in_file.read(raw_buffer, sizeof(raw_buffer));

		read_time += (GetTimeMs64() - t1);
		t2 = GetTimeMs64();

		#pragma omp parallel for num_threads(4) reduction(+:process_time)
		for (unsigned i=0; i < in_file.gcount(); ++i){
			uint8_t base, qual;
			uint32_t index = pos+i;
			uint32_t rev_index = N - index - 1;

			// Skip those sequences that did not pass the quality filter.
			if (tile_filter[index] == 0)
				continue;

			if (raw_buffer[i] == 0){ // no call if all 0 bits
				base = 4; // will be converted to 'N'
				qual = 0; // no quality
			}
			else{
				base = (raw_buffer[i] & 3);  // mask the lower 2 bits
				qual = ((raw_buffer[i] >> 2) & 63); // get bits 3-8
			}

			// Convert values to base and PHRED score.
			char baseChar = numToDNA(base);  // convert 0->A,..,3->T
			char qualChar = (char)(qual+33); // PHRED score to ASCII

			// If tile scanned first time, set id and runInfo.
			if (buffer->at(rev_index).id.size()==0){
				std::stringstream tmp;
				tmp << tile_num << "_" << (pos+i);
				buffer->at(rev_index).id = tmp.str();
				buffer->at(rev_index).runContainer = runInfo;
				runInfo->runInfoList.at(rev_index) = std::make_shared<SeqClassifyInfo>();
				buffer->at(rev_index).readInfo 	= runInfo->runInfoList.at(rev_index);
				//runInfo->at(rev_index) = buffer->at(rev_index).readInfo;

				// Create new classification info if empty, read old one otherwise
				/*if (runInfo->empty())
					buffer->at(rev_index).readInfo 	= std::unique_ptr<SeqClassifyInfo>(new SeqClassifyInfo());
				else
					;// read old info here*/
			}

			// Save the qualities into the read buffer.
			buffer->at(rev_index).seq += baseChar;
			buffer->at(rev_index).quals += qualChar;

			process_time += (GetTimeMs64() - t2);
		}

		pos += in_file.gcount();
	}

	//LOG("Tile processing time " << process_time + read_time << "\n");
	//std::cout << "Tile processing time " << process_time << ". Read time " << read_time << "\n";
	return true;
}

// Comparator function for sorting the cycles (C1.1,..,C201.1,..)
bool cyclePathCompare(const fs::path &p1, const fs::path &p2){
	std::string s1 = p1.filename().string();
	std::string s2 = p2.filename().string();

	int i1 = std::stoi(s1.substr(1, s1.find(".")-1));
	int i2 = std::stoi(s2.substr(1, s2.find(".")-1));

	return i1 < i2;
}

// Given a tile, advance it to the next tile number.
// (I.e. 1216 -> 1301 or 1316 -> 2101)
int getNextTile(int tile){
	// Get the components of the tile number
	int i1 = tile / 1000;               // Flowcell side
	int i2 = (tile - 1000*i1) / 100;    // Swath of the cell on this side
	int i3 = tile - 1000*i1 - 100*i2;   // One of 16 positions on this swath

	int overflow = 0;

	i3 +=1;
	if (i3 > 16){
		i3 = 1;
		overflow = 1;
	}

	i2 += overflow;
	overflow = 0;
	if (i2 > 3){
		i2 = 1;
		overflow = 1;
	}

	i1 += overflow;
	if (i1 > 2){
		return 0;
	}

	return i1*1000+i2*100+i3;
}

//-----------------------------------------------------------
//              BCL READER class
//-----------------------------------------------------------
void BCLReader::init(string filename){
	basecalls_path = fs::path(filename);
	if (!fs::exists(basecalls_path)){
		err(EX_NOINPUT, "%s not found.", filename.c_str());
	}

	if (fs::is_regular_file(basecalls_path)){
		valid = false;
		err(EX_NOINPUT, "File given. 'BaseCalls' directory expected.");
	}

	// Scan target directory and get all paths for each of the lane folders.
	getFilePaths();

	// We have two values for false. "_valid" is used internally and
	// is false if no new buffer can be filled. "valid" is used for the
	// external interface and is false if the last buffer is exhausted.
	// _valid will always be set to false before valid
	valid = true;
	_valid = true;
}


BCLReader::BCLReader(string filename)
: tile_num(1101), lane_num(1), read_length(-1) {
	this->init(filename);
}

BCLReader::BCLReader(string filename, int length)
: tile_num(1101), lane_num(1), read_length(length) {
	this->init(filename);
}

/*BCLReader::BCLReader(string filename, int length, std::vector<SeqClassifyInfo> *runInfoList)
: runInfoList(runInfoList), tile_num(1101), lane_num(1), read_length(length) {
	this->init(filename);
}*/

DNASequence BCLReader::next_sequence() {
	DNASequence dna;

	// We can't read anything anymore and have noting buffered either -> end.
	if (!_valid && concurrentBufferQueue.empty() && sequenceBuffer->empty()){
		valid = false;
		return dna;
	}

	// We execute this for the first time, start a sequence reader.
	if (sequenceBuffer == nullptr){
		LOG("Spawning new sequence reader thread. (" << tile_num << ")\n");
		fillSequenceBuffer();
	}

	// If the sequence buffer is empty, we get the next buffer, which
	// should have been filled by a concurrent thread. If the next buffer
	// is not ready yet, the function blocks there and waits for it.
	if (sequenceBuffer == nullptr || sequenceBuffer->empty()){

		// spawn a writer thread for temp information, if such information exists
		// TODO: Ensure that writing only is done when nothing is being read (to improve IO performance)
		// TODO: Refactor start of writer thread into own function.
		if (sequenceBuffer != nullptr && sequenceBuffer->empty()){
			saveRunInfo();
		}

#ifdef DEBUG
		if (concurrentBufferQueue.empty())
			LOG("\nWaiting for buffer to be filled...\n");
#endif

		// Get buffered reads.
		sequenceBuffer = std::move(concurrentBufferQueue.pop()); // this is blocking

		LOG(std::cout << "Found buffer. Size: " << sequenceBuffer->size() << "\n");

		// After consuming a buffer from the queue, spawn a new sequence reader.
		if (_valid){
			LOG("Spawning new sequence reader thread. (" << tile_num << ")\n");
			fillSequenceBuffer();
		}
	}

	//std::cout << sequenceBuffer->size() << "\n";

	// If the sequence is empty, it did not pass
	// the quality filter and will be skipped.
	do{
		dna = std::move(sequenceBuffer->back());
		sequenceBuffer->pop_back();
	} while (!sequenceBuffer->empty() && dna.seq.empty());

	// If we have exhausted the remaining buffer without
	// finding an unfiltered sequence, recursively call
	// next_sequence to create a new buffer and search again.
	if (sequenceBuffer->empty() && dna.seq.empty())
		return this->next_sequence();

	return dna;
}

bool BCLReader::is_valid() {
	return valid;
}


// Fill a buffer with sequences of one tile.
// Each call of this function will advance the tile number. If all
// If all tiles in a lane are processed, will continue with tile 1 of next lane.
bool BCLReader::fillSequenceBuffer(){
	// Set directory pointer to the next lane folder (L00..) if we
	// have read all tiles in the current lane (i.e. reached tile 2316).
	// 0 is a pseudo-tile name indicating that we have read all tiles.
	if (tile_num == 0){
		tile_num = 1101;

		// If we don't have another lane directory, end scan.
		if ((size_t)(lane_num-1 + 1) >= cyclePaths.size())
			_valid = false;
		else
			lane_num++;
	}

	// If current lane directory cannot be found.
	if (!_valid)
		return false;

	// Start thread that reads the BCL tile file into a buffer.
	std::thread sequenceReader(&BCLReader::addSequenceBuffer, this, lane_num, tile_num);

	// Tile processing thread runs independently, we will block later to wait.
	sequenceReader.detach();


	// Get the next tile number (processed in the next call of this function.)
	tile_num = getNextTile(tile_num);

	// Check if the next tile exists. If not, end reading.
	// XXX: Makes it abort when less than 96 tiles are found. Always intended?
	std::string tile_str(cyclePaths[lane_num-1][0].string() + "/s_" + std::to_string(lane_num) + "_" + std::to_string(tile_num) + ".bcl");
	if (tile_num != 0 && !fs::exists(fs::path(tile_str)))
		_valid = false;
	return false;

	return true;
}

// Creates a new sequence buffer, fills it with reads from the tile numbered
// "tile_num" and adds it to the Queue of buffers which Kraken consumes.
void BCLReader::addSequenceBuffer(int lane_num, int tile_num){
	LOG("Entered addSequenceBuffer " << tile_num << "\n");
	// Create new buffer to hold the reads.
	std::unique_ptr<TDNABuffer> buffer(new TDNABuffer());
	// Create new buffer for the classification state.
	std::shared_ptr<RunInfoContainer> runInfo(new RunInfoContainer(0, lane_num, tile_num));

	// create tile string for path construction
	std::string tile_str("/s_" + std::to_string(lane_num) + "_" + std::to_string(tile_num));

	LOG(tile_str << "\n");
	std::cout << tile_str << "\n";

	// Make index of .filter file for current tile.
	std::string s(lanePaths[lane_num-1].string() + tile_str + ".filter");


	// If temporary progress file exists, read it
	// otherwise the SeqClassifyInfo is initialized empty later
	if (fs::exists(tmpPaths[lane_num][tile_num])){
		// wait for writing to finish, just in case
		// (don't have to unlock again, since mutex is destroyed after)
		//writeLocks[lane_num][tile_num].lock();
		// put reading here


	}

	std::vector<bool> tile_filter;
	scanFilter(fs::path(s), tile_filter);

	// Process current tile in each cycle directory.
	for (unsigned i=0; i < cyclePaths[lane_num-1].size(); ++i){
		std::string s(cyclePaths[lane_num-1][i].string() + tile_str + ".bcl");

		// Add bases of tile in the current cycle to the buffered reads.
		if (scanTile(tile_num, fs::path(s), tile_filter, runInfo, buffer) == false){
			_valid = false;
		}
	}

	// Add the read buffer to the Queue which holds the precomputed buffers.
	// (It is a concurrent queue, so multiple threads can access it.)
	concurrentBufferQueue.push(std::move(buffer));
	concurrentRunInfoQueue.push(std::move(runInfo));
	LOG("Ended addSequenceBuffer" << tile_num << "\n");
}

void BCLReader::saveRunInfo(){
	std::shared_ptr<RunInfoContainer> tmp = std::move(concurrentRunInfoQueue.pop());

	// Lock the temporary file until the thread finished writing.
	writeLocks[tmp->lane_num] = std::unordered_map<int, std::mutex>();
	writeLocks[tmp->lane_num][tmp->tile_num].lock();
	std::thread RunInfoWriter(&BCLReader::writeInfo, this, std::move(tmp));
	RunInfoWriter.detach();
}

void BCLReader::writeInfo(std::shared_ptr<RunInfoContainer> runInfoContainer){
	std::ofstream ofs(tmpPaths[runInfoContainer->lane_num][runInfoContainer->tile_num].string());

		// wait for the runInfoList
		std::cout << std::this_thread::get_id() << " - Waiting for the release of the write lock...\n";
		runInfoContainer->processing_lock.lock();
		std::cout << std::this_thread::get_id() << " - Obtained the write lock...\n";
	    {
	        boost::archive::binary_oarchive oa(ofs);
	        for (TRunInfoList::iterator it = runInfoContainer->runInfoList.begin();
	        		it != runInfoContainer->runInfoList.end(); ++it){
	        	oa << *(*it);
	        }
	    }

	    // Note: We don't have to release processing_lock again since writeRunInfo is destroyed.
	    // Release the writeLock to indicate that writing is finished.
	    writeLocks[runInfoContainer->lane_num][runInfoContainer->tile_num].unlock();
}

// Utility functions

void BCLReader::getFilePaths(){
	copy_if(fs::directory_iterator(basecalls_path), fs::directory_iterator(), back_inserter(lanePaths),
			[&](const fs::path& p){
		return (p.filename().string()[0] == 'L');
	}
	);

	sort(lanePaths.begin(), lanePaths.end());

	for (fs::path &path : lanePaths){
		cyclePaths.push_back(vector<fs::path>());
		copy_if(fs::directory_iterator(path), fs::directory_iterator(), back_inserter(cyclePaths.back()),
				[&](const fs::path& p){
					return (fs::is_directory(p) && p.filename().string()[0] == 'C');
				}
		);

		sort(cyclePaths.back().begin(), cyclePaths.back().end(), cyclePathCompare);
		copy(cyclePaths.back().begin(), cyclePaths.back().end(), ostream_iterator<fs::path>(cout, "\n"));
	}

	std::cout << "Found " << lanePaths.size() << " lane directories.\n";
	copy(lanePaths.begin(), lanePaths.end(), ostream_iterator<fs::path>(cout, "\n"));

	std::cout << "Setting up paths for temporary directories.\n";
	for (size_t i=0; i < lanePaths.size(); ++i){
		tmpPaths.insert(std::make_pair<int, TTilePathMap>(i+1, TTilePathMap()));

		size_t j=1101;
		do{
			std::string s(lanePaths[i].string() + "/" + std::to_string(j) + ".filter");
			tmpPaths[i+1][j] = fs::path(s);
			std::cout << tmpPaths[i+1][j].string() << "\n";
		} while( (j = getNextTile(j)) );
	}
}

} // namespace
