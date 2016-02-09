CXX = g++
#CXX = clang++
#LDFLAGS=-L/home/RKID1/strauchb/library/lib
override CXXFLAGS += -Wall -fopenmp -std=c++11
PROGS = db_sort set_lcas classify make_seqid_to_taxid_map db_shrink

.PHONY: all install clean

all: CXX += -O3
all: $(PROGS)

debug: CXX += -DDEBUG -g
debug: $(PROGS)

install: $(PROGS)
	cp $(PROGS) $(KRAKEN_DIR)

clean:
	rm -f $(PROGS) *.o

db_shrink: krakendb.o quickfile.o

db_sort: krakendb.o quickfile.o

set_lcas: krakendb.o quickfile.o krakenutil.o seqreader.o filemanager.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) set_lcas.cpp krakendb.o quickfile.o krakenutil.o seqreader.o filemanager.o -o set_lcas $(LOADLIBES) $(LDLIBS) -lboost_serialization -lboost_system -lboost_filesystem

classify: classify.cpp krakendb.o quickfile.o krakenutil.o seqreader.o filemanager.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) classify.cpp krakendb.o quickfile.o krakenutil.o seqreader.o filemanager.o -o classify $(LOADLIBES) $(LDLIBS) -lboost_serialization -lboost_system -lboost_filesystem

make_seqid_to_taxid_map: quickfile.o

krakenutil.o: krakenutil.cpp krakenutil.hpp
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp

seqreader.o: seqreader.cpp seqreader.hpp quickfile.hpp filemanager.hpp
	$(CXX) $(CXXFLAGS) -c seqreader.cpp 

quickfile.o: quickfile.cpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c quickfile.cpp
	
filemanager.o: filemanager.cpp filemanager.hpp
	$(CXX) $(CXXFLAGS) -c filemanager.cpp