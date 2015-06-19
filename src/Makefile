CXX = g++
#CXX = clang++
CXXFLAGS = -Wall -fopenmp -std=c++0x
PROGS = db_sort set_lcas classify make_file_to_taxon_map db_shrink

.PHONY: all install clean

all: CXX += -DDEBUG -O3
all: $(PROGS)

debug: CXX += -DDEBUG -g
debug: $(PROGS)

install: $(PROGS)
	cp $(PROGS) $(KRAKEN_DIR)

clean:
	rm -f $(PROGS) *.o

db_shrink: krakendb.o quickfile.o

db_sort: krakendb.o quickfile.o

set_lcas: krakendb.o quickfile.o krakenutil.o seqreader.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) set_lcas.cpp krakendb.o quickfile.o krakenutil.o seqreader.o -o set_lcas $(LOADLIBES) $(LDLIBS) -lboost_system -lboost_filesystem

classify: krakendb.o quickfile.o krakenutil.o seqreader.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) classify.cpp krakendb.o quickfile.o krakenutil.o seqreader.o -o classify $(LOADLIBES) $(LDLIBS) -lboost_system -lboost_filesystem

make_file_to_taxon_map: quickfile.o

krakenutil.o: krakenutil.cpp krakenutil.hpp
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp

seqreader.o: seqreader.cpp seqreader.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c seqreader.cpp

quickfile.o: quickfile.cpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c quickfile.cpp
