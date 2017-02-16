cmake_minimum_required(VERSION 3.6)

add_definitions(-Wall -O3 -std=gnu++1z -mavx2 -mlzcnt -pthread -lnuma)
add_subdirectory(integerlist-utils)
add_subdirectory(pcm)
add_executable(numapr numa-PageRank.C ligra.h polymer.h polymer-wgh.h graph.h utils.h IO.h parallel.h gettime.h quickSort.h ./pcm/cpucounters.cpp ./pcm/msr.cpp ./pcm/pci.cpp ./pcm/client_bw.cpp)
target_link_libraries(numapr cmp_multi)
