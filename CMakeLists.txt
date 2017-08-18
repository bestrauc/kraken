cmake_minimum_required(VERSION 3.7)
project(kraken)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp" )


set(SOURCE_FILES
        src/classify.cpp
        src/filemanager.cpp
        src/filemanager.hpp
        src/krakendb.cpp
        src/krakendb.hpp
        src/kraken_headers.hpp
        src/krakenutil.cpp
        src/krakenutil.hpp
        src/Queue.h
        src/quickfile.cpp
        src/quickfile.hpp
        src/seqreader.cpp
        src/bclreader.cpp
        src/seqreader.hpp)


find_package(Boost REQUIRED COMPONENTS serialization system filesystem)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(classify ${SOURCE_FILES})
target_link_libraries(classify ${Boost_LIBRARIES})
