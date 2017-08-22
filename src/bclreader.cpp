//
// Created by benjamin on 18.08.17.
//

// for timing and debug
#include "seqreader.hpp"

using namespace std;
namespace fs = boost::filesystem;

namespace kraken {
  //-----------------------------------------------------------
  //              BCL READER helper functions
  //-----------------------------------------------------------
  // BCLReader helper functions
  char numToDNA(int base) {
    switch (base) {
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

  bool scanFilter(const path &filter_path, vector<bool> &tile_index) {
    ifstream in_file;
    in_file.open(filter_path.c_str(), ios_base::in | ios_base::binary);
    in_file.seekg(8);

    uint32_t N;
    in_file.read((char *) &N, 4);

    tile_index.resize(N);

    tile_index.assign(istreambuf_iterator<char>(in_file), istreambuf_iterator<char>());

    return true;
  }

  // Scan the tile given in the tile_path and save sequences into the buffer
  bool scanTile(int tile_num, const path &tile_path, TSeqClassifyInfoList &seqClassifyInfo,
                shared_ptr<RunInfoContainer> &runInfo, unique_ptr<WorkUnit> &buffer) {
    ifstream in_file;
    in_file.open(tile_path.c_str(), ios_base::in | ios_base::binary);

    if (!in_file.is_open()) {
      cout << "Unable to open file: " << tile_path << "\n";
      return false;
    }

    //TSeqClassifyInfoList &seqClassifyInfo = runInfo->seqClassifyInfo;

    // Process file
    uint32_t N; // number of clusters (sequences)
    uint32_t pos = 0;
    static char raw_buffer[16384]; //or 4096?

    in_file.read((char *) &N, 4);

    // This should hold if all files are correct and were correctly read
    //assert(N == tile_filter.size());

    buffer->seqs.resize(N);

    bool old_classinfo = true;

    // If the classification metadata for the sequences is not yet initialized, resize it and insert .
    if (seqClassifyInfo.size() != N) {
      seqClassifyInfo.resize(N);
      old_classinfo = false;
      runInfo->runsize = N;
    }

    buffer->runContainer = runInfo;

    // Variables for timing measurement. (On one line to comment it out easily.)
    uint64_t t1 = 0, t2 = 0, read_time = 0, process_time = 0;

    while (!in_file.eof()) {
      //t1 = GetTimeMs64();

      in_file.read(raw_buffer, sizeof(raw_buffer));

      //read_time += (GetTimeMs64() - t1);
      //t2 = GetTimeMs64();

#pragma omp parallel for num_threads(4) reduction(+:process_time)
      for (unsigned i = 0; i < in_file.gcount(); ++i) {
        uint8_t base, qual;
        uint32_t index = pos + i;
        uint32_t rev_index = N - index - 1;

        // Skip those sequences that did not pass the quality filter.
        //if (tile_filter[index] == 0)
        //	continue;

        if (raw_buffer[i] == 0) { // no call if all 0 bits
          base = 4; // will be converted to 'N'
          qual = 0; // no quality
        } else {
          base = (raw_buffer[i] & 3);  // mask the lower 2 bits
          qual = ((raw_buffer[i] >> 2) & 63); // get bits 3-8
        }

        // Convert values to base and PHRED score.
        char baseChar = numToDNA(base);  // convert 0->A,..,3->T
        char qualChar = (char) (qual + 33); // PHRED score to ASCII

        // If tile scanned first time, set id and runInfo.
        if (buffer->seqs.at(rev_index).id.size() == 0) {
          stringstream tmp;
          tmp << tile_num << "_" << index;
          buffer->seqs.at(rev_index).id = tmp.str();

          if (!old_classinfo)
            seqClassifyInfo.at(rev_index) = make_shared<SeqClassifyInfo>();

          buffer->seqs.at(rev_index).readInfo = seqClassifyInfo.at(rev_index);
          seqClassifyInfo.at(rev_index)->pos = runInfo->processed_nt;
        }

        // Save the qualities into the read buffer.
        buffer->seqs.at(rev_index).seq += baseChar;
        buffer->seqs.at(rev_index).quals += qualChar;
        buffer->seqs.at(rev_index).readInfo->processed_len++;

        //process_time += (GetTimeMs64() - t2);
      }

      pos += in_file.gcount();
    }

    //LOG("Tile processing time " << process_time + read_time << "\n");
    //std::cout << "Tile processing time " << process_time << ". Read time " << read_time << "\n";
    return true;
  }

  //-----------------------------------------------------------
  //              BCL READER class
  //-----------------------------------------------------------
  void BCLReader::init(string filename) {

    // We have two values for false. "_valid" is used internally and
    // is false if no new buffer can be filled. "valid" is used for the
    // external interface and is false if the last buffer is exhausted.
    // _valid will always be set to false before valid
    _valid = fileManager.is_valid();
    valid = _valid;

    // start thread that waits to write temporary sequence info
    //std::thread RunInfoWrite//r(&BCLReader::saveRunInfo, this);
    //RunInfoWriter.detach();

    //fileManager.cyclePaths.size()
  }

  BCLReader::BCLReader(std::string file_name, int length, int start_cycle, int step,
                       std::vector<int> target_tiles, std::vector<int> target_lanes)
      : fileManager(file_name, length, start_cycle, step, target_tiles, target_lanes) {
    //, writerThread(&BCLReader::saveRunInfo, this) {
    this->init(file_name);
  }

  size_t BCLReader::next_workunit(size_t work_nt_size, WorkUnit &work_unit) {
    if (!_valid && concurrentBufferQueue.empty() && sequenceBuffer->seqs.empty()) {
      valid = false;
      return 0;
    }

    // fill a sequence buffer if none exists
    if (sequenceBuffer == nullptr) {
      fillSequenceBuffer();
    }

    // get new reads from buffer and start a new reader afterwards
    if (sequenceBuffer == nullptr || sequenceBuffer->seqs.empty()) {
      // Get buffered reads.
      sequenceBuffer = move(concurrentBufferQueue.pop()); // this is blocking

      // After consuming a buffer from the queue, spawn a new sequence reader (unless there is nothing left).
      if (_valid) {
        fillSequenceBuffer();
      }
    }

    size_t total_nt = 0;

    while (total_nt < work_nt_size && !sequenceBuffer->seqs.empty()) {
      work_unit.runContainer = sequenceBuffer->runContainer;
      work_unit.seqs.emplace_back(sequenceBuffer->seqs.back());
      sequenceBuffer->seqs.pop_back();
      total_nt += work_unit.seqs.back().seq.size();
    }

    return total_nt;
  }

// TODO: Not implemented for BCLReader. Use next_workunit to load whole tiles for use in classify.cpp.
  DNASequence BCLReader::next_sequence() {
    DNASequence dna;

    return dna;
  }

  bool BCLReader::is_valid() {
    return valid;
  }

  // Fill a buffer with sequences of one tile.
  // Each call of this function will advance the tile number. If all
  // If all tiles in a lane are processed, will continue with tile 1 of next lane.
  bool BCLReader::fillSequenceBuffer() {
    // get the next tile number
    TileInfo tile = fileManager.getTile();

    // If we are at the end.
    if (!fileManager.is_valid()) {
      cout << "At end\n";
      _valid = false;
      return false;
    }

    // Start thread that reads the BCL tile file into a buffer.
    thread sequenceReader(&BCLReader::addSequenceBuffer, this, tile);
    //	addSequenceBuffer(tile);

    // Tile processing thread runs independently, we will block later to wait.
    sequenceReader.detach();

    return true;
  }

  // Creates a new sequence buffer, fills it with reads from the tile numbered
  // "tile_num" and adds it to the Queue of buffers which Kraken consumes.
  void BCLReader::addSequenceBuffer(TileInfo tile) {
    // Create new buffer to hold the reads.
    unique_ptr<WorkUnit> buffer(new WorkUnit());


    if (runInfoMap[tile.lane_num][tile.tile_num] == nullptr)
      runInfoMap[tile.lane_num][tile.tile_num] = make_shared<RunInfoContainer>(tile.lane_num, tile.tile_num,
                                                                               tile.last_cycle);

    if (!runInfoMap[tile.lane_num][tile.tile_num]->processing_lock.try_lock()) {
      cout << "Waiting for last tile \n";
      runInfoMap[tile.lane_num][tile.tile_num]->processing_lock.lock();
    }

    //runInfoMap[tile.lane_num][tile.tile_num]->processing_lock.lock();
    runInfoMap[tile.lane_num][tile.tile_num]->lane_num = tile.lane_num;
    runInfoMap[tile.lane_num][tile.tile_num]->tile_num = tile.tile_num;
    runInfoMap[tile.lane_num][tile.tile_num]->processed_nt = tile.last_cycle;
    runInfoMap[tile.lane_num][tile.tile_num]->count = 0;

    // create tile string for path construction
    string tile_str("/s_" + to_string(tile.lane_num) + "_" + to_string(tile.tile_num));

    // Make index of .filter file for current tile.
    //std::string s(lanePaths[lane_num-1].string() + tile_str + ".filter");
    //std::vector<bool> tile_filter;
    //scanFilter(fs::path(s), tile_filter);

    cout << "Filling sequence buffer: " << tile.first_cycle << " .. " << tile.last_cycle << "\n";

    // Process current tile in each cycle directory.
    for (int i = tile.first_cycle; i < tile.last_cycle; ++i) {
      string s(fileManager.cyclePaths[tile.lane_num - 1][i].string() + tile_str + ".bcl");
//      std::cout << s << " - \n";

      // Add bases of tile in the current cycle to the buffered reads.
      if (scanTile(tile.tile_num, path(s), runMap[tile.lane_num][tile.tile_num],
                   runInfoMap[tile.lane_num][tile.tile_num], buffer) == false) {
        _valid = false;
      }
    }

    // Add the read buffer to the Queue which holds the precomputed buffers.
    // (It is a concurrent queue, so multiple threads can access it.)
    concurrentBufferQueue.push(move(buffer));
    //concurrentRunInfoQueue.push(runInfoMap[tile.lane_num][tile.tile_num]);
  }

  void BCLReader::saveRunInfo() {
    //std::cout << "Waiting to save temp files\n";
    while (valid) {
      shared_ptr<RunInfoContainer> tmp = move(concurrentRunInfoQueue.pop());
      cout << "Received temp run info\n";

      // Lock the temporary file until the thread finished writing.
      writeLocks[tmp->lane_num] = unordered_map<int, mutex>();
      writeLocks[tmp->lane_num][tmp->tile_num].lock();
      cout << "Locked writer for temp files\n";
      writeInfo(move(tmp));

      //std::thread RunInfoWriter(&BCLReader::writeInfo, this, std::move(tmp));
      //RunInfoWriter.detach();
    }
  }

  void BCLReader::readInfo(shared_ptr<RunInfoContainer> runInfoContainer) {
    ifstream ifs(fileManager.tmpPaths[runInfoContainer->lane_num][runInfoContainer->tile_num].string());

    {
      // TODO: make a compression archive
      // e.g.: http://stackoverflow.com/questions/4961155/boostiostream-zlib-compressing-multiple-files-into-one-archive
//      boost::archive::binary_iarchive ia(ifs);
      // ia >> runInfoContainer->runInfoList;
    }

    // Release the writeLock to indicate that writing is finished.
    cout << "Unlocking temp file reader\n";
    writeLocks[runInfoContainer->lane_num][runInfoContainer->tile_num].unlock();
  }

  void BCLReader::writeInfo(shared_ptr<RunInfoContainer> runInfoContainer) {
    ofstream ofs(fileManager.tmpPaths[runInfoContainer->lane_num][runInfoContainer->tile_num].string());

    // wait for the runInfoList
    cout << this_thread::get_id() << " - Waiting for the release of the write lock...\n";
    runInfoContainer->processing_lock.lock();
    cout << this_thread::get_id() << " - Obtained the write lock...\n";

    cout << runInfoContainer->processed_nt << " " << runInfoContainer->lane_num << " " << runInfoContainer->tile_num
         << " " << fileManager.tmpPaths[runInfoContainer->lane_num][runInfoContainer->tile_num].string() << "\n";
    {
      // TODO: make a compression archive
      // e.g.: http://stackoverflow.com/questions/4961155/boostiostream-zlib-compressing-multiple-files-into-one-archive
//      boost::archive::binary_oarchive oa(ofs);
      // oa << runInfoContainer->runInfoList;
    }

    // Note: We don't have to release processing_lock again since writeRunInfo is destroyed.
    // Release the writeLock to indicate that writing is finished.
    writeLocks[runInfoContainer->lane_num][runInfoContainer->tile_num].unlock();
  }
}