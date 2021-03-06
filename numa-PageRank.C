﻿/* 
 * This code is part of the project "NUMA-aware Graph-structured Analytics"
 * 
 *
 * Copyright (C) 2014 Institute of Parallel And Distributed Systems (IPADS), Shanghai Jiao Tong University
 *     All rights reserved 
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 * 
 * For more about this software, visit:
 *
 *     http://ipads.se.sjtu.edu.cn/projects/polymer.html
 *
 */

#include "polymer.h"
#include "gettime.h"
#include "math.h"

#include <pthread.h>
#include <sys/mman.h>
#include <numa.h>
#include <sys/syscall.h>

#include <sched.h>
#include <string>

#include "pcm/cpucounters.h"

//#include <papi.h>
#define NUM_EVENTS 3

using namespace std;

#define PAGE_SIZE (4096)
int N_NODES = 0;
int N_CORES_PER_NODE = 0;
int N_USE_NODES = 0;
int N_USE_CORES_PER_NODE = 0;

//int CORES_PER_NODE = 6;
int NODE_USED = -1;

volatile int shouldStart = 0;

double *p_curr_global = NULL;
double *p_next_global = NULL;

double *p_ans = NULL;
int vPerNode = 0;
//int numOfNode = 0;

bool needResult = false;

pthread_barrier_t barr;
pthread_barrier_t global_barr;
pthread_mutex_t mut;

volatile int global_counter = 0;
volatile int global_toggle = 0;

vertices *Frontier;

struct PR_F {
    double *p_curr, *p_next;
    //vertex *V;
    int rangeLow;
    int rangeHi;
    intE *out_degrees;

    PR_F(double *_p_curr, double *_p_next, intE *_out_degrees, int _rangeLow, int _rangeHi) :
            p_curr(_p_curr), p_next(_p_next), out_degrees(_out_degrees), rangeLow(_rangeLow), rangeHi(_rangeHi) {}

    inline void *nextPrefetchAddr(intT index) {
        return &p_curr[index];
    }

    inline bool update(intT s, intT d) { //update function applies PageRank equation
        p_next[d] += p_curr[s] / out_degrees[s];
        return 1;
    }

    inline double getCurrVal(intT i) {
        return p_curr[i];
    }

    inline bool updateValVer(intT s, double val, intT d) {
        writeAdd(&p_next[d], val / out_degrees[s]);
        return true;
    }

    inline bool updateAtomic(intT s, intT d) { //atomic Update
        writeAdd(&p_next[d], p_curr[s] / out_degrees[s]);
        //return (p_curr[s] / V[s].getOutDegree()) >= 0;
        /*
        if (d == 110101) {
            cout << "Update from " << s << "\t" << std::scientific << std::setprecision(9) << p_curr[s]/V[s].getOutDegree() << " -- " << p_next[d] << "\n";
        }
        */
        return 1;
    }

    inline void initFunc(void *dataPtr, intT d) {
        *(double *) dataPtr = 0.0;
    }

    inline bool reduceFunc(void *dataPtr, intT s, bool print_info = false) {
        *(double *) dataPtr += p_curr[s] / (double) out_degrees[s];
        if (print_info) {
            //cout << "reduce: " << s << " " << std::scientific << std::setprecision(9) << p_curr[s] / (double)V[s].getOutDegree() << " " << p_curr[s] << " " << V[s].getOutDegree() << *(double *)dataPtr << "\n";
        }
        return true;
    }

    inline bool combineFunc(void *dataPtr, intT d) {
        double val = *(double *) dataPtr;
        writeAdd((double *) &p_next[d], val);
        /*
        if (d == 77) {
            cout << "combine result: " << std::scientific << std::setprecision(9) << val << " " << p_next[d] << "\n";
        }
        */
        return true;
    }

    inline bool cond(intT d) { return true; } //does nothing
};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F {
    double damping;
    double addedConstant;
    double *p_curr;
    double *p_next;

    PR_Vertex_F(double *_p_curr, double *_p_next, double _damping, intT n) :
            p_curr(_p_curr), p_next(_p_next),
            damping(_damping), addedConstant((1 - _damping) * (1 / (double) n)) {}

    inline bool operator()(intT i) {
        p_next[i] = damping * p_next[i] + addedConstant;
        return 1;
    }
};

//resets p
struct PR_Vertex_Reset {
    double *p_curr;

    PR_Vertex_Reset(double *_p_curr) :
            p_curr(_p_curr) {}

    inline bool operator()(intT i) {
        p_curr[i] = 0.0;
        return 1;
    }
};

struct PR_worker_arg {
    void *GA;
    int maxIter;
    int tid;
    int numOfNode;
    int rangeLow;
    int rangeHi;
};

struct PR_subworker_arg {
    void *GA;
    int maxIter;
    int tid;
    int subTid;
    int startPos;
    int endPos;
    int rangeLow;
    int rangeHi;
    double **p_curr_ptr;
    double **p_next_ptr;
    double damping;
    pthread_barrier_t *node_barr;
    LocalFrontier *localFrontier;
    volatile int *barr_counter;
    volatile int *toggle;
};

template<class F>
bool *edgeMapDenseForwardOTHER(graph0 GA, vertices *frontier, F f, LocalFrontier *next,
                               intE *nghs, bool part = false, int start = 0, int end = 0) {
    intT numVertices = GA.n;

    int currNodeNum = 0;
    bool *currBitVector = frontier->getArr(currNodeNum);
    int nextSwitchPoint = frontier->getSize(0);
    int currOffset = 0;
    int counter = 0;

    intT outEdgesCount = 0;
    bool *nextB = next->b;

    int startPos = 0;
    int endPos = numVertices;
    if (part) {
        startPos = start;
        endPos = end;
        currNodeNum = frontier->getNodeNumOfIndex(startPos);
        //printf("nodeNum: %d %d\n", currNodeNum, endPos);
        currBitVector = frontier->getArr(currNodeNum);
        nextSwitchPoint = frontier->getOffset(currNodeNum + 1);
        currOffset = frontier->getOffset(currNodeNum);
    }

    GA.map_out_nghs([&f, startPos, endPos](uintT ngh, uintT curr, uintT idx) {
      if(startPos <= curr && curr < endPos) {
        f.updateValVer(curr, f.getCurrVal(curr), ngh);
      }
    });

  /*
    for (long i = startPos; i < endPos; i++) {
        if (i == nextSwitchPoint) {
            currOffset += frontier->getSize(currNodeNum);
            nextSwitchPoint += frontier->getSize(currNodeNum + 1);
            currNodeNum++;
            currBitVector = frontier->getArr(currNodeNum);
            //printf("OK\n");
        }
        if (currBitVector[i - currOffset]) {
            auto val = f.getCurrVal(i);
#if TYPE < 5
            G[i].traverseOutNgh([&f, &next, i, val](uintT ngh) {
//                if(f.cond(ngh) && f.updateValVer(i, val, ngh)){
//                    next->setBit(ngh, true);
//                }
                f.updateValVer(i, val, ngh);
            });
#else
            GA.traverseOutNgh(i, [&f, i, val](uintT ngh) {
                f.updateValVer(i, val, ngh);
            });
#endif
        }
    }
    */
    return NULL;
}

void *PageRankSubWorker(void *arg) {
    PR_subworker_arg *my_arg = (PR_subworker_arg *) arg;
    graph0 &GA = *(graph0 *)my_arg->GA;
    const intT n = GA.n;
    int maxIter = my_arg->maxIter;
    int tid = my_arg->tid;
    int subTid = my_arg->subTid;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    int P_CORES_PER_NODE = N_CORES_PER_NODE / 2;
    int offset = subTid < P_CORES_PER_NODE ? 0 : (N_NODES - 1) * P_CORES_PER_NODE;
    int core = tid * P_CORES_PER_NODE + subTid + offset;
    CPU_SET(core, &cpuset);
    sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set_t), &cpuset);
    //cerr << "On " + to_string(sched_getcpu()) + "\n";

    intE *nghs = (intE *)numa_alloc_local(sizeof(intE)*2000);

    pthread_barrier_t *local_barr = my_arg->node_barr;
    LocalFrontier *output = my_arg->localFrontier;

    double *p_curr = *(my_arg->p_curr_ptr);
    double *p_next = *(my_arg->p_next_ptr);

    double damping = my_arg->damping;
    int currIter = 0;
    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    int start = my_arg->startPos;
    int end = my_arg->endPos;

    Custom_barrier globalCustom(&global_counter, &global_toggle, Frontier->numOfNodes);
    Custom_barrier localCustom(my_arg->barr_counter, my_arg->toggle, N_USE_CORES_PER_NODE);

    Subworker_Partitioner subworker(N_USE_CORES_PER_NODE);
    subworker.tid = tid;
    subworker.subTid = subTid;
    subworker.dense_start = start;
    subworker.dense_end = end;
    subworker.global_barr = &global_barr;
    subworker.local_custom = localCustom;
    subworker.subMaster_custom = globalCustom;

    if (subTid == 0) {
        Frontier->getFrontier(tid)->m = rangeHi - rangeLow;
    }

    pthread_barrier_wait(local_barr);
    pthread_barrier_wait(&global_barr);

    while (1) {
        if (maxIter > 0 && currIter >= maxIter)
            break;
        currIter++;
        if (subTid == 0)
            Frontier->calculateNumOfNonZero(tid);
        if (subTid == 0) {
            //{parallel_for(long i=output->startID;i<output->endID;i++) output->setBit(i, false);}
        }

        pthread_barrier_wait(&global_barr);
        //pthread_barrier_wait(local_barr);

        //edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output,0,DENSE_FORWARD, false, true, subworker);
        clearLocalFrontier(output, subworker.tid, subworker.subTid, subworker.numOfSub);
        output->sparseCounter = 0;
        subworker.globalWait();
        //edgeMapDenseReduce(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output,false,subworker);

        subworker.localWait();
        struct timeval startT, midT, endT;
        struct timezone tz = {0, 0};
        gettimeofday(&startT, &tz);
        //edgeMapDenseForward(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output, true, subworker.dense_start, subworker.dense_end);
        edgeMapDenseForwardOTHER(GA, Frontier, PR_F(p_curr, p_next, GA.out_degrees, rangeLow, rangeHi), output, nghs, true,
                                 subworker.dense_start, subworker.dense_end);
        //edgeMapDenseForwardDynamic(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output, subworker);
        gettimeofday(&midT, &tz);
        subworker.localWait();
        gettimeofday(&endT, &tz);
        if (subworker.isSubMaster()) {
            double time1 = ((double) startT.tv_sec) + ((double) startT.tv_usec) / 1000000.0;
            double time2 = ((double) midT.tv_sec) + ((double) midT.tv_usec) / 1000000.0;
            double time3 = ((double) endT.tv_sec) + ((double) endT.tv_usec) / 1000000.0;
            double duration0 = time2 - time1;
            double duration1 = time3 - time2;
//            printf("time of %d: %lf-%lf on %d, %d, %d, %d, %d\n", subworker.tid * CORES_PER_NODE + subworker.subTid,
//                   duration0, duration1, sched_getcpu(), rangeLow, rangeHi, subworker.dense_start, subworker.dense_end);
        }

        output->isDense = true;

        pthread_barrier_wait(&global_barr);
        //pthread_barrier_wait(local_barr);
        if (subTid == 0) {
            //printf("next active: %d\n", output->m);
        }

        vertexMap(Frontier, PR_Vertex_F(p_curr, p_next, damping, n),
            tid, subTid, N_USE_CORES_PER_NODE);
        //vertexCounter(GA, output, tid, subTid, CORES_PER_NODE);
        output->m = 1;

        pthread_barrier_wait(&global_barr);
        //pthread_barrier_wait(local_barr);

        vertexMap(Frontier, PR_Vertex_Reset(p_curr), tid, subTid, N_USE_CORES_PER_NODE);
        pthread_barrier_wait(&global_barr);
        //pthread_barrier_wait(local_barr);
        swap(p_curr, p_next);
        if (subworker.isSubMaster()) {
            pthread_barrier_wait(&global_barr);
            switchFrontier(tid, Frontier, output);
        } else {
            output = Frontier->getFrontier(tid);
            pthread_barrier_wait(&global_barr);
        }
        //pthread_barrier_wait(local_barr);
    }

    if (subworker.isMaster()) {
        p_ans = p_curr;
    }
    pthread_barrier_wait(local_barr);

    numa_free(nghs, sizeof(intE)*2000);
    return NULL;
}

pthread_barrier_t timerBarr;

void *PageRankThread(void *arg) {
    PR_worker_arg *my_arg = (PR_worker_arg *) arg;
    graph0 &GA = *(graph0 *) my_arg->GA;
    int maxIter = my_arg->maxIter;
    int tid = my_arg->tid;

    char nodeString[10];
    sprintf(nodeString, "%d", tid);
    struct bitmask *nodemask = numa_parse_nodestring(nodeString);
    numa_bind(nodemask);

    int rangeLow = my_arg->rangeLow;
    int rangeHi = my_arg->rangeHi;

    /*
    if (tid == 0) {
        //printf("average is: %lf\n", GA.m / (float) (my_arg->numOfNode));
    }
    pthread_barrier_wait(&barr);
    intT degreeSum = 0;
    for (intT i = rangeLow; i < rangeHi; i++) {
        degreeSum += GA.in_degrees[i];
    }
    cerr << to_string(tid) + " : degree count: " + to_string(degreeSum) + "\n";
    */

    graph0 localGraph = graphFilter2Direction0(GA, rangeLow, rangeHi);

    pthread_barrier_wait(&barr);
    if (tid == 0)
        GA.del();
    pthread_barrier_wait(&barr);

    int sizeOfShards[N_USE_CORES_PER_NODE];

    subPartitionByDegree(localGraph, N_USE_CORES_PER_NODE, sizeOfShards, sizeof(double), true, true);
    //intT localDegrees = (intT *)malloc(sizeof(intT) * localGraph.n);

    for (int i = 0; i < N_USE_CORES_PER_NODE; i++) {
        //printf("subPartition: %d %d: %d\n", tid, i, sizeOfShards[i]);
    }

    while (shouldStart == 0);
    pthread_barrier_wait(&timerBarr);
    cerr << "over filtering\n";
    /*
    if (0 != __cilkrts_set_param("nworkers","1")) {
    printf("set failed: %d\n", tid);
    }
    */

    const intT n = GA.n;
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    int numOfT = my_arg->numOfNode;

    int blockSize = rangeHi - rangeLow;

    //printf("blockSizeof %d: %d low: %d high: %d\n", tid, blockSize, rangeLow, rangeHi);

    double one_over_n = 1 / (double) n;


    double *p_curr = p_curr_global;
    double *p_next = p_next_global;
    bool *frontier = (bool *) numa_alloc_local(sizeof(bool) * blockSize);

    /*
    double* p_curr = (double *)malloc(sizeof(double) * blockSize);
    double* p_next = (double *)malloc(sizeof(double) * blockSize);
    bool* frontier = (bool *)malloc(sizeof(bool) * blockSize);
    */

    /*
    if (tid == 0)
    startTime();
    */
    double mapTime = 0.0;
    struct timeval start, end;
    struct timezone tz = {0, 0};

    for (intT i = rangeLow; i < rangeHi; i++) p_curr[i] = one_over_n;
    for (intT i = rangeLow; i < rangeHi; i++) p_next[i] = 0; //0 if unchanged
    for (intT i = 0; i < blockSize; i++) frontier[i] = true;
    if (tid == 0)
        Frontier = new vertices(numOfT);

    //printf("register %d: %p\n", tid, frontier);

    LocalFrontier *current = new LocalFrontier(frontier, rangeLow, rangeHi);

    bool *next = (bool *) numa_alloc_local(sizeof(bool) * blockSize);
    for (intT i = 0; i < blockSize; i++) next[i] = false;
    LocalFrontier *output = new LocalFrontier(next, rangeLow, rangeHi);

    pthread_barrier_wait(&barr);

    Frontier->registerFrontier(tid, current);

    pthread_barrier_wait(&barr);

    if (tid == 0)
        Frontier->calculateOffsets();

    pthread_barrier_t localBarr;
    pthread_barrier_init(&localBarr, NULL, N_USE_CORES_PER_NODE + 1);

    int startPos = 0;

    pthread_t subTids[N_USE_CORES_PER_NODE];

    volatile int local_custom_counter = 0;
    volatile int local_toggle = 0;

    for (int i = 0; i < N_USE_CORES_PER_NODE; i++) {
        PR_subworker_arg *arg = (PR_subworker_arg *) malloc(sizeof(PR_subworker_arg));
        arg->GA = (void *) (&localGraph);
        arg->maxIter = maxIter;
        arg->tid = tid;
        arg->subTid = i;
        arg->rangeLow = rangeLow;
        arg->rangeHi = rangeHi;
        arg->p_curr_ptr = &p_curr;
        arg->p_next_ptr = &p_next;
        arg->damping = damping;
        arg->node_barr = &localBarr;
        arg->localFrontier = output;

        arg->barr_counter = &local_custom_counter;
        arg->toggle = &local_toggle;

        arg->startPos = startPos;
        arg->endPos = startPos + sizeOfShards[i];
        startPos = arg->endPos;

        pthread_create(&subTids[i], NULL, PageRankSubWorker, (void *) arg);
    }

    pthread_barrier_wait(&barr);

    pthread_barrier_wait(&localBarr);

    pthread_barrier_wait(&localBarr);

    pthread_barrier_wait(&barr);
    intT round = 0;
    /*
    while(1){
    if (maxIter > 0 && round >= maxIter)
        break;
    round++;

    pthread_barrier_wait(&localBarr);

    //edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V,rangeLow,rangeHi),output,GA.m/20,DENSE_FORWARD);
    pthread_barrier_wait(&localBarr);
    
    //vertexMap(Frontier, PR_Vertex_F(p_curr, p_next, damping, n), tid);      

    pthread_barrier_wait(&barr);
    pthread_barrier_wait(&localBarr);
    //vertexMap(Frontier,PR_Vertex_Reset(p_curr), tid);
    pthread_barrier_wait(&localBarr);

    swap(p_curr,p_next);
    if (tid == 0) {
        p_ans = p_curr;
    }

    switchFrontier(tid, Frontier, output);

    pthread_barrier_wait(&localBarr);    
    pthread_barrier_wait(&barr);
    }
    */
    return NULL;
}

struct PR_Hash_F {
    int shardNum;
    int vertPerShard;
    int n;

    PR_Hash_F(int _n, int _shardNum) : n(_n), shardNum(_shardNum), vertPerShard(_n / _shardNum) {}

    inline int hashFunc(int index) {
        if (index >= shardNum * vertPerShard) {
            return index;
        }
        int idxOfShard = index % shardNum;
        int idxInShard = index / shardNum;
        return (idxOfShard * vertPerShard + idxInShard);
    }

    inline int hashBackFunc(int index) {
        if (index >= shardNum * vertPerShard) {
            return index;
        }
        int idxOfShard = index / vertPerShard;
        int idxInShard = index % vertPerShard;
        return (idxOfShard + idxInShard * shardNum);
    }
};

tuple<SystemCounterState, SystemCounterState> PageRank(graph0 &GA, int maxIter) {
    N_NODES = numa_num_configured_nodes();
    N_CORES_PER_NODE = numa_num_configured_cpus() / N_NODES;

//    numOfNode = numa_num_configured_nodes();
    vPerNode = GA.n / N_USE_NODES;
//    CORES_PER_NODE = numa_num_configured_cpus() / numOfNode;
    if (NODE_USED != -1)
        N_USE_NODES = NODE_USED;
    pthread_barrier_init(&barr, NULL, N_USE_NODES);
    pthread_barrier_init(&timerBarr, NULL, N_USE_NODES + 1);
    pthread_barrier_init(&global_barr, NULL, N_USE_CORES_PER_NODE * N_USE_NODES);
    pthread_mutex_init(&mut, NULL);
    int sizeArr[N_USE_NODES];
    PR_Hash_F hasher(GA.n, N_USE_NODES);
    //graphHasher(GA, hasher);
    //graphAllEdgeHasher(GA, hasher);
    partitionByDegree0(GA, N_USE_NODES, sizeArr, sizeof(double));
    /*
    intT vertPerPage = PAGESIZE / sizeof(double);
    intT subShardSize = ((GA.n / numOfNode) / vertPerPage) * vertPerPage;
    for (int i = 0; i < numOfNode - 1; i++) {
    sizeArr[i] = subShardSize;
    }
    sizeArr[numOfNode - 1] = GA.n - subShardSize * (numOfNode - 1);
    */
    /*
    int accum = 0;
    for (int i = 0; i < N_USE_NODES; i++) {
        intT degreeSum = 0;
        for (intT j = accum; j < accum + sizeArr[i]; j++) {
            degreeSum += GA.in_degrees[j];
        }
        cerr << to_string(i) + ": degree sum: " + to_string(degreeSum) + "\n";
        accum += sizeArr[i];
    }
    */
    //return;

    p_curr_global = (double *) mapDataArray(N_USE_NODES, sizeArr, sizeof(double));
    p_next_global = (double *) mapDataArray(N_USE_NODES, sizeArr, sizeof(double));

    cerr << "start create " + to_string(N_USE_NODES) << "threads\n";
    pthread_t tids[N_USE_NODES];
    int prev = 0;
    for (int i = 0; i < N_USE_NODES; i++) {
        PR_worker_arg *arg = (PR_worker_arg *) malloc(sizeof(PR_worker_arg));
        arg->GA = (void *) (&GA);
        arg->maxIter = maxIter;
        arg->tid = i;
        arg->numOfNode = N_USE_NODES;
        arg->rangeLow = prev;
        arg->rangeHi = prev + sizeArr[i];
        prev = prev + sizeArr[i];

        pthread_create(&tids[i], NULL, PageRankThread, (void *) arg);
    }
    shouldStart = 1;

    pthread_barrier_wait(&timerBarr);
    //nextTime("Graph Partition");
    nextTime("partition over");
    auto mid = getSystemCounterState();
    cerr << "all created\n";
    for (int i = 0; i < N_USE_NODES; i++) {
        pthread_join(tids[i], NULL);
    }
    nextTime("PageRank");
    auto end = getSystemCounterState();

    if (needResult) {
        for (intT i = 0; i < GA.n; i++) {
            cerr << i << "\t" << std::scientific << std::setprecision(9) << p_ans[hasher.hashFunc(i)] << "\n";
            //cout << i << "\t" << std::scientific << std::setprecision(9) << p_ans[i] << "\n";
        }
    }

    return make_tuple(mid, end);
}

int parallel_main(int argc, char *argv[]) {
    char *iFile;
    bool binary = false;
    bool symmetric = false;
    int maxIter = 20;
    needResult = false;

    iFile = argv[1];
    maxIter = atoi(argv[2]);
    N_USE_NODES = atoi(argv[3]);
    N_USE_CORES_PER_NODE = atoi(argv[4]);

//    if (argc > 3) NODE_USED = atoi(argv[3]);
//    if (argc > 4) if ((string) argv[4] == (string) "-result") needResult = true;
//    if (argc > 5) if ((string) argv[5] == (string) "-s") symmetric = true;
//    if (argc > 6) if ((string) argv[6] == (string) "-b") binary = true;
    numa_set_interleave_mask(numa_all_nodes_ptr);

    PCM *m = PCM::getInstance();
    m->program(PCM::DEFAULT_EVENTS, NULL);
    auto start = getSystemCounterState();

    startTime();
    if (symmetric) {
//        graph<symmetricVertex> G =
//                readGraph<symmetricVertex>(iFile, symmetric, binary);
//        PageRank(G, maxIter);
        //G.del();
    } else {
        //graph<asymmetricVertex> G =
                //readGraph<asymmetricVertex>(iFile, symmetric, binary);
        graph0 G = readGraphFromFile0(iFile, symmetric);

        auto tup = PageRank(G, maxIter);
        auto mid = get<0>(tup);
        auto end = get<1>(tup);

        cout << getBytesReadFromMC(start, mid)
            << "\t" << getBytesReadFromMC(mid, end) << endl;
        cout << getBytesWrittenToMC(start, mid)
            << "\t" << getBytesWrittenToMC(mid, end) << endl;
        cout << getIORequestBytesFromMC(start, mid)
            << "\t" << getIORequestBytesFromMC(mid, end) << endl;
        cout << getL2CacheMisses(start, mid)
            << "\t" << getL2CacheMisses(mid, end) << endl;
        cout << getL3CacheMisses(start, mid)
            << "\t" << getL3CacheMisses(mid, end) << endl;
        cout << to_string(getL2CacheHitRatio(start, mid))
            << "\t" << to_string(getL2CacheHitRatio(mid, end)) << endl;
        cout << to_string(getL3CacheHitRatio(start, mid))
            << "\t" << to_string(getL3CacheHitRatio(mid, end)) << endl;

        //G.del();
    }
    return 0;
}
