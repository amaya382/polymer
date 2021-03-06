ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif
#CILK = 1
# # no compare and swap!
# ifdef OPENMP
# PCC = g++
# PCFLAGS = -fopenmp -mcx16 -O3 -DOPENMP $(INTT) $(INTE)

ifdef CILK
PCC = g++
#-cilk
PCFLAGS = -fcilkplus -lcilkrts -O2 -DCILK $(INTT) $(INTE)
PLFLAGS = -fcilkplus -lcilkrts

else ifdef MKLROOT
PCC = icpc
PCFLAGS = -O3 -DCILKP $(INTT) $(INTE)

else
PCC = g++
PCFLAGS = -O3 $(INTT) $(INTE)
endif

#PCFLAGS = -fcilkplus -lcilkrts -O2 -DCILK $(INTT) $(INTE)
#PLFLAGS = -fcilkplus -lcilkrts

COMMON= ligra.h polymer.h polymer-wgh.h graph.h utils.h IO.h parallel.h gettime.h quickSort.h

ALL= DegreeCount ConvertToBinary #PartitionGraphToEdgeList
#MYAPPS= numa-BP numa-PageRank numa-PageRank-bin numa-PageRank-pull numa-PageRank-write numa-PageRankDelta numa-Components numa-BFS numa-BFS-async-pipe numa-SPMV numa-BellmanFord ConvertToJSON ConvertTmp
MYAPPS = numa-PageRank
MYHEADER= ligra-rewrite.h ligra-numa.h
LIBS_I_NEED= -pthread -lnuma

#all: $(ALL) $(MYAPPS)
all: $(MYAPPS)

debug: PCFLAGS = -fcilkplus -lcilkrts -O0 -g -DCILK $(INTT) $(INTE)
debug: all

% : %.C $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $< $(LIBS_I_NEED) ./pcm/cpucounters.cpp ./pcm/msr.cpp ./pcm/pci.cpp ./pcm/client_bw.cpp -I ./pcm -std=c++14 -mavx2 -mlzcnt -DTYPE=6

.PHONY : clean

clean :
	rm -f *.o $(ALL) $(MYAPPS)

