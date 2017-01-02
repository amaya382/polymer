#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "parallel.h"

using namespace std;

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

struct symmetricVertex {
    intE *neighbors;
    intT degree;
    intT fakeDegree;

    void del() { free(neighbors); }

    symmetricVertex(intE *n, intT d) : neighbors(n), degree(d) {}

    uintE getInNeighbor(intT j) { return neighbors[j]; }

    uintE getOutNeighbor(intT j) { return neighbors[j]; }

    intE *getInNeighborPtr() { return neighbors; }

    intE *getOutNeighborPtr() { return neighbors; }

    intT getInDegree() { return degree; }

    intT getOutDegree() { return degree; }

    intT getFakeInDegree() { return fakeDegree; }

    intT getFakeDegree() { return fakeDegree; }

    void setInNeighbors(intE *_i) { neighbors = _i; }

    void setOutNeighbors(intE *_i) { neighbors = _i; }

    void setInDegree(intT _d) { degree = _d; }

    void setOutDegree(intT _d) { degree = _d; }

    void setFakeDegree(intT _d) { fakeDegree = _d; }

    void setFakeInDegree(intT _d) { fakeDegree = _d; }

    void flipEdges() {}
};


#include <vector>
#include <cstring>
#include <bitset>
template<typename uint_t>
uint64_t encode(vector <uint_t> &in, vector <uint8_t> &out) {
    // out[0..n_chunks-1] are reserved for flags
    // out[n_chunks..] are used to store compressed data
    uint64_t size = in.size();
    uint64_t n_chunks = (size + 3) / 4;
    uint64_t used = n_chunks;

    for (uint64_t i = 0; i < n_chunks; i++) {
        bitset<8> flags;
        uint64_t block = i * 4;
        for (uint8_t j = 0; j < 4 && block + j < size; j++) {
            if (in[block + j] <= 0xFF) {
                memcpy(&out[used], &in[block + j], 1);
                used++;
            } else if (in[block + j] <= 0xFFFF) {
                memcpy(&out[used], &in[block + j], 2);
                used += 2;
                flags = 1 << (3 - j) * 2;
            } else if (in[block + j] <= 0xFFFFFF) {
                memcpy(&out[used], &in[block + j], 4);
                used += 4;
                flags = 2 << (3 - j) * 2;
            } else {
                memcpy(&out[used], &in[block + j], 8);
                used += 8;
                flags = 3 << (3 - j) * 2;
            }
        }
        out[i] = static_cast<uint8_t>(flags.to_ulong());
    }

    return used; // byte
}

template<typename uint_t>
uint64_t encode0(uint_t *in, uint64_t size, uint8_t *out) {
    // out[0..n_chunks-1] are reserved for flags
    // out[n_chunks..] are used to store compressed data
    uint64_t n_chunks = (size + 3) / 4;
    uint64_t used = n_chunks;

    for (uint64_t i = 0; i < n_chunks; i++) {
        bitset<8> flags;
        uint64_t block = i * 4;
        for (uint8_t j = 0; j < 4 && block + j < size; j++) {
            if (in[block + j] <= 0xFF) {
                memcpy(&out[used], &in[block + j], 1);
                used++;
            } else if (in[block + j] <= 0xFFFF) {
                memcpy(&out[used], &in[block + j], 2);
                used += 2;
                flags = 1 << (3 - j) * 2;
            } else if (in[block + j] <= 0xFFFFFF) {
                memcpy(&out[used], &in[block + j], 4);
                used += 4;
                flags = 2 << (3 - j) * 2;
            } else {
                memcpy(&out[used], &in[block + j], 8);
                used += 8;
                flags = 3 << (3 - j) * 2;
            }
        }
        out[i] = static_cast<uint8_t>(flags.to_ulong());
    }

    return used; // byte
}

template<typename uint_t>
void decode(vector <uint8_t> &in, uint64_t size, vector <uint_t> &out) {
    uint64_t n_chunks = (size + 3) / 4;
    uint64_t used = n_chunks;

    for (uint64_t i = 0; i < n_chunks; i++) {
        uint64_t block = i * 4;
        for (uint8_t j = 0; j < 4 && block + j < size; j++) {
            switch (in[i] >> (3 - j) * 2 & 0b00000011) {
                case 0:
                    memcpy(&out[block + j], &in[used], 1);
                    used++;
                    break;
                case 1:
                    memcpy(&out[block + j], &in[used], 2);
                    used += 2;
                    break;
                case 2:
                    memcpy(&out[block + j], &in[used], 4);
                    used += 4;
                    break;
                case 3:
                    memcpy(&out[block + j], &in[used], 8);
                    used += 8;
                    break;
            }
        }
    }
}

template<typename uint_t>
void decode0(uint8_t *in, uint64_t size, uint_t *out) {
    uint64_t n_chunks = (size + 3) / 4;
    uint64_t used = n_chunks;

    for (uint64_t i = 0; i < n_chunks; i++) {
        uint64_t block = i * 4;
        for (uint8_t j = 0; j < 4 && block + j < size; j++) {
            switch (in[i] >> (3 - j) * 2 & 0b00000011) {
                case 0:
                    memcpy(&out[block + j], &in[used], 1);
                    used++;
                    break;
                case 1:
                    memcpy(&out[block + j], &in[used], 2);
                    used += 2;
                    break;
                case 2:
                    memcpy(&out[block + j], &in[used], 4);
                    used += 4;
                    break;
                case 3:
                    memcpy(&out[block + j], &in[used], 8);
                    used += 8;
                    break;
            }
        }
    }
}

template<typename uint_t>
void decode1(uint8_t *in, uint64_t size, vector<uint_t> &out) {
    uint64_t n_chunks = (size + 3) / 4;
    uint64_t used = n_chunks;

    for (uint64_t i = 0; i < n_chunks; i++) {
        uint64_t block = i * 4;
        for (uint8_t j = 0; j < 4 && block + j < size; j++) {
            switch (in[i] >> (3 - j) * 2 & 0b00000011) {
                case 0:
                    memcpy(&out[block + j], &in[used], 1);
                    used++;
                    break;
                case 1:
                    memcpy(&out[block + j], &in[used], 2);
                    used += 2;
                    break;
                case 2:
                    memcpy(&out[block + j], &in[used], 4);
                    used += 4;
                    break;
                case 3:
                    memcpy(&out[block + j], &in[used], 8);
                    used += 8;
                    break;
            }
        }
    }
}

struct asymmetricVertex {
    intE *inNeighbors;
    intE *outNeighbors;
    intT outDegree;
    intT fakeInDegree;
    intT fakeOutDegree;
    intT inDegree;

    uint8_t *in;
    uint8_t *out;

    void del() {}

    asymmetricVertex(intE *iN, intE *oN, intT id, intT od) : inNeighbors(iN), outNeighbors(oN), inDegree(id),
                                                             outDegree(od) {}

    uintE getInNeighbor(intT j) { return inNeighbors[j]; }

    uintE getOutNeighbor(intT j) { return outNeighbors[j]; }

    intE *getInNeighborPtr() { return inNeighbors; }

    intE *getOutNeighborPtr() { return outNeighbors; }

    intT getInDegree() { return inDegree; }

    intT getOutDegree() { return outDegree; }

    intT getFakeInDegree() { return fakeInDegree; }

    intT getFakeDegree() { return fakeOutDegree; }

    void setInNeighbors(intE *_i) {
        inNeighbors = _i;
        in = (uint8_t *) _i;
    }

    void setOutNeighbors(intE *_i) {
        outNeighbors = _i;
        out = (uint8_t *) _i;
    }

    void setInDegree(intT _d) { inDegree = _d; }

    void setOutDegree(intT _d) { outDegree = _d; }

    void setFakeDegree(intT _d) { fakeOutDegree = _d; }

    void setFakeInDegree(intT _d) { fakeInDegree = _d; }

    void flipEdges() {
        swap(inNeighbors, outNeighbors);
        swap(inDegree, outDegree);
    }

//    vector<uint8_t> in;
//    vector<uint8_t> out;

    void getOutNgh(vector<uintE> &data,uint64_t size){ decode1<uintE>(out,size,data);}
//    void setOutNgh(vector<uintE> &data){encode<uintE>(data,out);}
//    uint64_t setOutNgh(uintE *data, uint64_t size){return decode0<uintE>(out, size, data);}
//    void getOutNgh(uintE *data, uint64_t size){encode0<uintE>(data, size, out);}
};

template<class vertex>
struct graph {
    vertex *V;
    intT n;
    uintT m;
    intE *allocatedInplace;
    intE *inEdges;
    intT *flags;

    graph(vertex *VV, intT nn, uintT mm)
            : V(VV), n(nn), m(mm), allocatedInplace(NULL), flags(NULL) {}

    graph(vertex *VV, intT nn, uintT mm, intE *ai, intE *_inEdges = NULL)
            : V(VV), n(nn), m(mm), allocatedInplace(ai), inEdges(_inEdges), flags(NULL) {}

    void del() {
        if (flags != NULL) free(flags);
    }

    void transpose() {
        if (sizeof(vertex) == sizeof(asymmetricVertex)) {
            parallel_for (intT i = 0; i < n; i++) {
                V[i].flipEdges();
            }
        }
    }
};

struct symmetricWghVertex {
    //weights are stored in the entry after the neighbor ID
    //so size of neighbor list is twice the degree
    intE *neighbors;
    intT degree;
    intT fakeDegree;

    void del() { free(neighbors); }

    symmetricWghVertex(intE *n, intT d) : neighbors(n), degree(d) {}

    intE getInNeighbor(intT j) { return neighbors[2 * j]; }

    intE getOutNeighbor(intT j) { return neighbors[2 * j]; }

    intE *getInNeighborPtr() { return neighbors; }

    intE *getOutNeighborPtr() { return neighbors; }

    intE getInWeight(intT j) { return neighbors[2 * j + 1]; }

    intE getOutWeight(intT j) { return neighbors[2 * j + 1]; }

    intT getInDegree() { return degree; }

    intT getOutDegree() { return degree; }

    intT getFakeDegree() { return fakeDegree; }

    intT getFakeInDegree() { return fakeDegree; }

    void setInNeighbors(intE *_i) { neighbors = _i; }

    void setOutNeighbors(intE *_i) { neighbors = _i; }

    void setInDegree(intT _d) { degree = _d; }

    void setOutDegree(intT _d) { degree = _d; }

    void setFakeDegree(intT _d) { fakeDegree = _d; }

    void setFakeInDegree(intT _d) { fakeDegree = _d; }
};

struct asymmetricWghVertex {
    //weights are stored in the entry after the neighbor ID
    //so size of neighbor list is twice the degree
    intE *inNeighbors;
    intE *outNeighbors;
    intT inDegree;
    intT outDegree;
    intT fakeOutDegree;
    intT fakeInDegree;

    void del() {
        free(inNeighbors);
        free(outNeighbors);
    }

    asymmetricWghVertex(intE *iN, intE *oN, intT id, intT od) : inNeighbors(iN), outNeighbors(oN), inDegree(id),
                                                                outDegree(od) {}

    intE getInNeighbor(intT j) { return inNeighbors[2 * j]; }

    intE getOutNeighbor(intT j) { return outNeighbors[2 * j]; }

    intE *getInNeighborPtr() { return inNeighbors; }

    intE *getOutNeighborPtr() { return outNeighbors; }

    intE getInWeight(intT j) { return inNeighbors[2 * j + 1]; }

    intE getOutWeight(intT j) { return outNeighbors[2 * j + 1]; }

    intT getInDegree() { return inDegree; }

    intT getOutDegree() { return outDegree; }

    intT getFakeDegree() { return fakeOutDegree; }

    intT getFakeInDegree() { return fakeInDegree; }

    void setInNeighbors(intE *_i) { inNeighbors = _i; }

    void setOutNeighbors(intE *_i) { outNeighbors = _i; }

    void setInDegree(intT _d) { inDegree = _d; }

    void setOutDegree(intT _d) { outDegree = _d; }

    void setFakeDegree(intT _d) { fakeOutDegree = _d; }

    void setFakeInDegree(intT _d) { fakeInDegree = _d; }
};

template<class vertex>
struct wghGraph {
    vertex *V;
    intT n;
    uintT m;
    intE *allocatedInplace;
    intE *inEdges;
    intT *flags;

    wghGraph(vertex *VV, intT nn, uintT mm)
            : V(VV), n(nn), m(mm), allocatedInplace(NULL), flags(NULL) {}

    wghGraph(vertex *VV, intT nn, uintT mm, intE *ai, intE *_inEdges = NULL)
            : V(VV), n(nn), m(mm), allocatedInplace(ai), inEdges(_inEdges), flags(NULL) {}

    void del() {
        if (flags != NULL) free(flags);
        if (allocatedInplace == NULL)
            for (intT i = 0; i < n; i++) V[i].del();
        else { free(allocatedInplace); }
        free(V);
        if (inEdges != NULL) free(inEdges);
    }
};
