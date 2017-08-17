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
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include <unordered_set>
//#include "seqreader.hpp"

const size_t DEF_WORK_UNIT_SIZE = 500000;

using namespace std;
using namespace kraken;

void parse_command_line(int argc, char **argv);

void usage(int exit_code = EX_USAGE);

void process_file(char *filename);

void classify_sequence(DNASequence &dna, ostringstream &koss, ostringstream &coss, ostringstream &uoss,
                       uint64_t &seq_count);

void classify_partial_sequence(DNASequence &dna, uint32_t &call);

void classify_finalize(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss, uint32_t call);

string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig);

set<uint32_t> get_ancestry(uint32_t taxon);

void report_stats(struct timeval time1, struct timeval time2);

int Num_threads = 1;
string DB_filename, Index_filename, Nodes_filename;
bool Quick_mode = false;

// Kraken modifications here
enum Input_mode {
    FASTQ, FASTA, BCL
};
Input_mode File_input = FASTA;

size_t length = -1;
int max_tile = 2316;
// ------------------------

bool Print_classified = false;
bool Print_unclassified = false;
bool Print_kraken = true;
bool Populate_memory = false;
bool Only_classified_kraken_output = false;
uint32_t Minimum_hit_count = 1;
unordered_map<uint32_t, uint32_t> Parent_map;
unordered_map<uint32_t, string> taxLevel_map;
KrakenDB Database;
string Classified_output_file, Unclassified_output_file, Kraken_output_file;
ostream *Classified_output;
ostream *Unclassified_output;
ostream *Kraken_output;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;

uint64_t total_classified = 0;
uint64_t total_sequences = 0;
uint64_t total_bases = 0;

uint64_t GetTimeMs64() {
#ifdef _WIN32
    /* Windows */
    FILETIME ft;
    LARGE_INTEGER li;

    /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
     * to a LARGE_INTEGER structure. */
    GetSystemTimeAsFileTime(&ft);
    li.LowPart = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;

    uint64 ret = li.QuadPart;
    ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
    ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

    return ret;
#else
    /* Linux */
    struct timeval tv;

    gettimeofday(&tv, NULL);

    uint64_t ret = tv.tv_usec;
    /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
    ret /= 1000;

    /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
    ret += (tv.tv_sec * 1000);

    return ret;
#endif
}


int main(int argc, char **argv) {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif

    parse_command_line(argc, argv);
    if (!Nodes_filename.empty()) {
        Parent_map = build_parent_map(Nodes_filename, taxLevel_map);
    }

    if (Populate_memory)
        cerr << "Loading database... ";

    QuickFile db_file;
    db_file.open_file(DB_filename);
    if (Populate_memory)
        db_file.load_file();
    Database = KrakenDB(db_file.ptr());
    KmerScanner::set_k(Database.get_k());

    QuickFile idx_file;
    idx_file.open_file(Index_filename);
    if (Populate_memory)
        idx_file.load_file();
    KrakenDBIndex db_index(idx_file.ptr());
    Database.set_index(&db_index);

    if (Populate_memory)
        cerr << "complete." << endl;

    if (Print_classified) {
        if (Classified_output_file == "-")
            Classified_output = &cout;
        else
            Classified_output = new std::ofstream(Classified_output_file.c_str());
    }

    if (Print_unclassified) {
        if (Unclassified_output_file == "-")
            Unclassified_output = &cout;
        else
            Unclassified_output = new std::ofstream(Unclassified_output_file.c_str());
    }

    if (!Kraken_output_file.empty()) {
        if (Kraken_output_file == "-")
            Print_kraken = false;
        else
            Kraken_output = new std::ofstream(Kraken_output_file.c_str());
    } else
        Kraken_output = &cout;

    struct timeval tv1, tv2;
    gettimeofday(&tv1, NULL);
    for (int i = optind; i < argc; i++)
        process_file(argv[i]);
    gettimeofday(&tv2, NULL);

    Kraken_output->flush();

    report_stats(tv1, tv2);

    return 0;
}

void report_stats(struct timeval time1, struct timeval time2) {
    time2.tv_usec -= time1.tv_usec;
    time2.tv_sec -= time1.tv_sec;
    if (time2.tv_usec < 0) {
        time2.tv_sec--;
        time2.tv_usec += 1000000;
    }
    double seconds = time2.tv_usec;
    seconds /= 1e6;
    seconds += time2.tv_sec;

    cerr << "\r";
    fprintf(stderr,
            "%llu sequences (%.2f Mbp) processed in %.3fs (%.1f Kseq/m, %.2f Mbp/m).\n",
            (unsigned long long) total_sequences, total_bases / 1.0e6, seconds,
            total_sequences / 1.0e3 / (seconds / 60),
            total_bases / 1.0e6 / (seconds / 60));
    fprintf(stderr, "  %llu sequences classified (%.2f%%)\n",
            (unsigned long long) total_classified, total_classified * 100.0 / total_sequences);
    fprintf(stderr, "  %llu sequences unclassified (%.2f%%)\n",
            (unsigned long long) (total_sequences - total_classified),
            (total_sequences - total_classified) * 100.0 / total_sequences);
}

void check_workunit(WorkUnit &unit){
    std::unordered_set<int> lengths;
    std::unordered_set<int> poses;
    std::unordered_set<std::shared_ptr<RunInfoContainer> > pointers;

    for (auto& dna: unit){
        lengths.insert(dna.seq.size());
        poses.emplace(dna.readInfo->pos);
        pointers.insert(dna.runContainer);
    }

    std::cout << "Found sizes: \n";
    for (const auto &elem : lengths)
        std::cout << elem << " ";

    std::cout << "\nFound poses: \n";
    for (const auto &elem : poses)
        std::cout << elem << " ";

    std::cout << "\nFound pointers: \n";
    for (const auto &elem : pointers)
        std::cout << elem << " ";
    std::cout << "\n";
}

void process_file(char *filename) {
    string file_str(filename);
    DNASequenceReader *reader;
    DNASequence dna;

    if (File_input == FASTQ)
        reader = new FastqReader(file_str, length);
    else if (File_input == FASTA)
        reader = new FastaReader(file_str, length);
    else if (File_input == BCL)
        reader = new BCLReader(file_str, length, max_tile);
    else {
        cerr << "File type not recognized\n";
    }

    uint64_t tax_level_time = 0, read_time = 0, process_time = 0;

    // save a log of the classifications at each iteration
    // (iteration, call) -> count
    std::map<int, std::map<int, int> > classified_count;

    //uint64_t t = GetTimeMs64();
#pragma omp parallel reduction(+:read_time, process_time, tax_level_time)
    {
        WorkUnit work_unit;
        ostringstream kraken_output_ss, classified_output_ss, unclassified_output_ss;

        while (reader->is_valid()) {
            work_unit.clear();
            size_t total_nt = 0;
#pragma omp critical(get_input)
            {
                uint64_t t1 = GetTimeMs64();
                total_nt = reader->next_workunit(Work_unit_size, work_unit);
                read_time += (GetTimeMs64() - t1);
            }

            if (total_nt == 0)
                break;

            kraken_output_ss.str("");
            classified_output_ss.str("");
            unclassified_output_ss.str("");

            // sequence count for this batch of reads
            // for BCL files, count individually
            uint64_t seq_count = 0; //= (File_input != BCL) ? work_unit.size() : 0;

            std::cout << "Work unit size: " << work_unit.size() << "\n";
            std::cout << "Sequence length: " << work_unit[0].seq.size() << " to " << work_unit.back().seq.size() << "\n";

            check_workunit(work_unit);

            uint64_t count_pos = 0;

            uint64_t t1 = GetTimeMs64();
            for (size_t j = 0; j < work_unit.size(); j++) {
                bool full_length = work_unit[j].readInfo->pos >= length;
                int last_size = work_unit[j].seq.size();
                count_pos += full_length;

                classify_sequence(work_unit[j], kraken_output_ss, classified_output_ss,
                                  unclassified_output_ss, seq_count);
            }

            std::cout << count_pos << " counted.\n";

//            std::cout << seq_count << " of " << work_unit.size() << " at iteration " <<
//                      work_unit[0].readInfo->pos << " - " << work_unit.back().readInfo->pos << "\n";

            process_time += (GetTimeMs64() - t1);

#pragma omp critical(write_output)
            {
                if (Print_kraken)
                    (*Kraken_output) << kraken_output_ss.str();
                if (Print_classified)
                    (*Classified_output) << classified_output_ss.str();
                if (Print_unclassified)
                    (*Unclassified_output) << unclassified_output_ss.str();
                total_sequences += count_pos; //seq_count;
                total_bases += total_nt;
                cerr << "\rProcessed " << total_sequences << " sequences (" << total_bases << " bp) ...\n";
            }
        }
    }  // end parallel section

    std::cout << "Read: " << read_time << ", Process: " << process_time << ", Tax Level: " << tax_level_time << "\n";

    delete reader;
}

void classify_sequence(DNASequence &dna, ostringstream &koss, ostringstream &coss, ostringstream &uoss,
                       uint64_t &seq_count) {
    uint32_t call = 0;

    // classify the sequence we're given
    classify_partial_sequence(dna, call);

    // we can finalize if: Full fasta or fastq sequence given and called
    // or if the BCL file was sufficiently accurately classified
    bool tax_call = check_tax_level(call, "genus", Parent_map, taxLevel_map);

    if (!dna.runContainer || tax_call || dna.readInfo->pos >= length) {
//        std::cout << "Entering because " << !dna.runContainer << " " << tax_call << " " <<
//                  dna.readInfo->pos << "\n";
//        seq_count += (File_input == BCL);
        seq_count++;
        dna.readInfo->skip = true;
        classify_finalize(dna, koss, coss, uoss, call);
    }
}

void classify_partial_sequence(DNASequence &dna, uint32_t &call) {
    vector<uint32_t> taxa;
    vector<uint8_t> ambig_list;

    uint64_t *kmer_ptr;
    uint32_t taxon = 0;

    uint64_t current_bin_key;
    int64_t current_min_pos = 1;
    int64_t current_max_pos = 0;

    // check !first because k-mer size requirement may not hold for incremental steps
    if (!dna.readInfo->first || dna.seq.size() >= Database.get_k()) {
        KmerScanner scanner(dna.seq, 0, ~0, !dna.readInfo->first, dna.readInfo->last_ambig,
                            dna.readInfo->last_kmer);

        while ((kmer_ptr = scanner.next_kmer()) != NULL) {
            taxon = 0;
            if (scanner.ambig_kmer()) {
                dna.readInfo->ambig_list.push_back(1);
            } else {
                dna.readInfo->ambig_list.push_back(0);
                uint32_t *val_ptr = Database.kmer_query(
                        Database.canonical_representation(*kmer_ptr),
                        &current_bin_key,
                        &current_min_pos, &current_max_pos
                );

                taxon = val_ptr ? *val_ptr : 0;

                if (taxon) {
                    dna.readInfo->hit_counts[taxon]++;
                    if (Quick_mode && ++dna.readInfo->hits >= Minimum_hit_count)
                        break;
                }
            }

            dna.readInfo->last_kmer = *kmer_ptr;
            dna.readInfo->last_ambig = scanner.get_ambig();
            dna.readInfo->taxa.push_back(taxon);
        }

        if (Quick_mode)
            call = dna.readInfo->hits >= Minimum_hit_count ? taxon : 0;
        else {
            call = resolve_tree2(dna.readInfo->hit_counts, Parent_map);
        }

        dna.readInfo->first = false;
    }
}

void classify_finalize(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss, uint32_t call) {

    if (call)
#pragma omp atomic
        total_classified++;

    if (Print_unclassified || Print_classified) {
        ostringstream *oss_ptr = call ? &coss : &uoss;
        bool print = call ? Print_classified : Print_unclassified;
        if (print) {
            if (File_input == FASTQ) {
                (*oss_ptr) << "@" << dna.header_line << endl
                           << dna.seq << endl
                           << "+" << endl
                           << dna.quals << endl;
            } else {
                (*oss_ptr) << ">" << dna.header_line << endl
                           << dna.seq << endl;
            }
        }
    }

    if (!Print_kraken)
        return;

    if (call) {
        koss << "C\t";
    } else {
        if (Only_classified_kraken_output)
            return;

        koss << "U\t";
    }

    koss << dna.id << "\t" << call << "\t" << dna.readInfo->processed_len << "\t";

    if (Quick_mode) {
        koss << "Q:" << dna.readInfo->hits;
    } else {
        if (dna.readInfo->taxa.empty())
            koss << "0:0";
        else
            koss << hitlist_string(dna.readInfo->taxa, dna.readInfo->ambig_list);
    }

    koss << endl;
}


string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig) {
    int64_t last_code;
    int code_count = 1;
    ostringstream hitlist;

    if (ambig[0]) { last_code = -1; }
    else { last_code = taxa[0]; }

    for (size_t i = 1; i < taxa.size(); i++) {
        int64_t code;
        if (ambig[i]) { code = -1; }
        else { code = taxa[i]; }

        if (code == last_code) {
            code_count++;
        } else {
            if (last_code >= 0) {
                hitlist << last_code << ":" << code_count << " ";
            } else {
                hitlist << "A:" << code_count << " ";
            }
            code_count = 1;
            last_code = code;
        }
    }
    if (last_code >= 0) {
        hitlist << last_code << ":" << code_count;
    } else {
        hitlist << "A:" << code_count;
    }
    return hitlist.str();
}

set<uint32_t> get_ancestry(uint32_t taxon) {
    set<uint32_t> path;

    while (taxon > 0) {
        path.insert(taxon);
        taxon = Parent_map[taxon];
    }
    return path;
}

void parse_command_line(int argc, char **argv) {
    int opt;
    long long sig;

    if (argc > 1 && strcmp(argv[1], "-h") == 0)
        usage(0);
    while ((opt = getopt(argc, argv, "d:i:t:u:n:m:o:qfbC:U:Ml:x:")) != -1) {
        switch (opt) {
            case 'd' :
                DB_filename = optarg;
                break;
            case 'i' :
                Index_filename = optarg;
                break;
            case 't' :
                sig = atoll(optarg);
                if (sig <= 0)
                    errx(EX_USAGE, "can't use nonpositive thread count");
#ifdef _OPENMP
                if (sig > omp_get_num_procs())
                    errx(EX_USAGE, "thread count exceeds number of processors");
                Num_threads = sig;
                omp_set_num_threads(Num_threads);
#endif
                break;
            case 'n' :
                Nodes_filename = optarg;
                break;
            case 'q' :
                Quick_mode = true;
                break;
            case 'm' :
                sig = atoll(optarg);
                if (sig <= 0)
                    errx(EX_USAGE, "can't use nonpositive minimum hit count");
                Minimum_hit_count = sig;
                break;
            case 'f' :;
                File_input = FASTQ;
                break;
            case 'b' :
                File_input = BCL;
                break;
            case 'c' :
                Only_classified_kraken_output = true;
                break;
            case 'C' :
                Print_classified = true;
                Classified_output_file = optarg;
                break;
            case 'U' :
                Print_unclassified = true;
                Unclassified_output_file = optarg;
                break;
            case 'o' :
                Kraken_output_file = optarg;
                break;
            case 'u' :
                sig = atoll(optarg);
                if (sig <= 0)
                    errx(EX_USAGE, "can't use nonpositive work unit size");
                Work_unit_size = sig;
                break;
            case 'M' :
                Populate_memory = true;
                break;
            case 'l' :
                length = atoi(optarg);
                break;
            case 'x' :
                max_tile = atoi(optarg);
                break;
            default:
                usage();
                break;
        }
    }

    if (DB_filename.empty()) {
        cerr << "Missing mandatory option -d" << endl;
        usage();
    }
    if (Index_filename.empty()) {
        cerr << "Missing mandatory option -i" << endl;
        usage();
    }
    if (Nodes_filename.empty() && !Quick_mode) {
        cerr << "Must specify one of -q or -n" << endl;
        usage();
    }
    if (optind == argc) {
        cerr << "No sequence data files specified" << endl;
    }
}

void usage(int exit_code) {
    cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
         << endl
         << "Options: (*mandatory)" << endl
         << "* -d filename      Kraken DB filename" << endl
         << "* -i filename      Kraken DB index filename" << endl
         << "  -n filename      NCBI Taxonomy nodes file" << endl
         << "  -o filename      Output file for Kraken output" << endl
         << "  -t #             Number of threads" << endl
         << "  -u #             Thread work unit size (in bp)" << endl
         << "  -q               Quick operation" << endl
         << "  -m #             Minimum hit count (ignored w/o -q)" << endl
         << "  -C filename      Print classified sequences" << endl
         << "  -U filename      Print unclassified sequences" << endl
         << "  -l length        Length of reads (BCL mode)" << endl
         << "  -x max_tile      Maximum tile number (BCL mode, e.g. 1116)" << endl
         << "  -f               Input is in FASTQ format" << endl
         << "  -b               Input is in BCL format" << endl
         << "  -c               Only include classified reads in output" << endl
         << "  -M               Preload database files" << endl
         << "  -h               Print this message" << endl
         << endl
         << "At least one FASTA or FASTQ file must be specified." << endl
         << "Kraken output is to standard output by default." << endl;
    exit(exit_code);
}
