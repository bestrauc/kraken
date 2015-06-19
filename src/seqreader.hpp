/*
 * Copyright 2013-2014, Derrick Wood <dwood@cs.umd.edu>
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
  typedef struct {
    std::string id;
    std::string seq;
    std::string quals;
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
    DNASequence next_sequence();
    bool is_valid();

    private:
    Queue<std::unique_ptr<TDNABuffer> > concurrentBufferQueue;
    std::unique_ptr<TDNABuffer> sequenceBuffer;
    std::vector<path> cyclePaths;
    boost::filesystem::directory_iterator lane_dir_iter;
    path basecalls_path;

    int tile_num;
    int lane_num;
    bool valid, _valid;
    int read_length;

    bool fillSequenceBuffer();
    void addSequenceBuffer(int lane_num, int tile_num, std::vector<path> &cyclePaths);
    void getCyclePaths();
    void init(std::string filename);
  };
}

#endif
