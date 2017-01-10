#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "parallel.h"

#include <vector>
#include <cstring>
#include <functional>
#include <immintrin.h>

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

template<typename uint_t>
inline void decode(uint8_t *in, uint64_t size, uint_t *out) {
    uint64_t n_chunks = (size + 3) / 4;
    uint64_t used = n_chunks;

    for (uint64_t i = 0; i < n_chunks; i++) {
        uint64_t block = i * 4;
        for (uint8_t j = 0; j < 4 && block + j < size; j++) {
            uint8_t n_bytes = 0b00000001 << (in[i] >> (3 - j) * 2 & 0b00000011);
            out[block + j] = 0;
            memcpy(&out[block + j], &in[used], n_bytes);
            used += n_bytes;
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

    void setInNeighbors(intE *_i) { inNeighbors = _i; }

    void setOutNeighbors(intE *_i) { outNeighbors = _i; }

    void setInDegree(intT _d) { inDegree = _d; }

    void setOutDegree(intT _d) { outDegree = _d; }

    void setFakeDegree(intT _d) { fakeOutDegree = _d; }

    void setFakeInDegree(intT _d) { fakeInDegree = _d; }

    void flipEdges() {
        swap(inNeighbors, outNeighbors);
        swap(inDegree, outDegree);
    }

#if TYPE == 0
    // 4ALL
    template<typename Func>
    void traverseOutNgh(Func f){
        uint64_t n_chunks = (fakeOutDegree + 3) / 4;
        uint64_t used = n_chunks;

        uintT ngh = 0;
        for (uint64_t i = 0; i < n_chunks; i++) {
            uint64_t block = i * 4;
            for (uint8_t j = 0; j < 4 && block + j < fakeOutDegree; j++) {
                uint8_t n_bytes = 0b00000001 << (out[i] >> (3 - j) * 2 & 0b00000011);
                uint64_t mask = 0xFFFFFFFFFFFFFFFFull >> (8 - n_bytes) * 8;
                ngh += (reinterpret_cast<uintT *>(&out[used]))[0] & mask;
                f(ngh);
                used += n_bytes;
            }
        }
    }
#elif TYPE == 1
    // 2ALL
    template<typename Func>
    void traverseOutNgh(Func f){
        uint64_t n_chunks = (fakeOutDegree + 7) / 8;
        uint64_t used = n_chunks;

        uintT ngh = 0;
        for (uint64_t i = 0; i < n_chunks; i++) {
            uint64_t block = i * 8;
            for (uint8_t j = 0; j < 8 && block + j < fakeOutDegree; j++) {
                bool flag = out[i] >> (7 - j) & 0b00000001;
                if (flag) {
                    ngh += (reinterpret_cast<uintT *>(&out[used]))[0]
                        & 0xFFFFFFFFFFFFFFFFull;
                    f(ngh);
                    used += 8;
                } else {
                    ngh += (reinterpret_cast<uintT *>(&out[used]))[0]
                        & 0xFFFF;
                    f(ngh);
                    used += 2;
                }
            }
        }
    }
#elif TYPE == 2
    // HEAD4ALL
    template<typename Func>
    void traverseOutNgh(Func f){
        uint64_t n_chunks = (fakeOutDegree + 2) / 4;
        uint64_t used = n_chunks;

        uintT ngh = 0;
        if (fakeOutDegree > 0) {
            ngh += (reinterpret_cast<uintT *>(out))[0];
            f(ngh);
            used += 8;
        } else {
            return;
        }

        for (uint64_t i = 0; i < n_chunks; i++) {
            uint64_t block = i * 4;
            for (uint8_t j = 0; j < 4 && block + j < fakeOutDegree - 1; j++) {
                uint8_t n_bytes = 0b00000001 << (out[i + 8] >> (3 - j) * 2 & 0b00000011);
                uint64_t mask = 0xFFFFFFFFFFFFFFFFull >> (8 - n_bytes) * 8;
                ngh += (reinterpret_cast<uintT *>(&out[used]))[0] & mask;
                f(ngh);
                used += n_bytes;
            }
        }
    }
#elif TYPE == 3
    // HEAD2ALL
    template<typename Func>
    void traverseOutNgh(Func f){
        uint64_t n_chunks = (fakeOutDegree + 6) / 8;
        uint64_t used = n_chunks;

        uintT ngh = 0;
        if (fakeOutDegree > 0) {
            ngh += (reinterpret_cast<uintT *>(out))[0];
            f(ngh);
            used += 8;
        } else {
            return;
        }

        for (uint64_t i = 0; i < n_chunks; i++) {
            uint64_t block = i * 8;
            for (uint8_t j = 0; j < 8 && block + j < fakeOutDegree - 1; j++) {
                bool flag = out[i + 8] >> (7 - j) & 0b00000001;
                if (flag) {
                    ngh += (reinterpret_cast<uintT *>(&out[used]))[0]
                        & 0xFFFFFFFFFFFFFFFFull;
                    f(ngh);
                    used += 8;
                } else {
                    ngh += (reinterpret_cast<uintT *>(&out[used]))[0]
                        & 0xFFFF;
                    f(ngh);
                    used += 2;
                }
            }
        }
    }
#elif TYPE == 4
    // No optimization
    template<typename Func>
    void traverseOutNgh(Func f){
        uintT ngh = 0;
        auto p = reinterpret_cast<uintT *>(out);
        for (uint64_t i = 0; i < fakeOutDegree; i++){
            ngh += p[i];
            f(ngh);
        }
    }
#elif TYPE == 5
    // TODO: probably too slow
    template<typename uint_t = uint32_t>
    uint_t *unpack(uint8_t *packed, uint8_t n_bits, uint_t *xs, uint8_t size = 8) {
        constexpr auto BIT_SIZE_OF_T = sizeof(uint_t) * 8;
        auto *_packed = reinterpret_cast<uint_t *>(packed);
        auto mask = 0xFFFFFFFF;
        uint8_t acc = n_bits;
        for (auto i = 0; i < size; i++, acc += n_bits) {
            auto block = acc / BIT_SIZE_OF_T;
            auto surplus = acc % BIT_SIZE_OF_T < n_bits ? acc % BIT_SIZE_OF_T % n_bits : 0;
            if (surplus) {
                xs[i] = _packed[block - 1] >> (BIT_SIZE_OF_T - (n_bits - surplus));
                xs[i] |= (_packed[block] & (mask >> (BIT_SIZE_OF_T - surplus))) << (n_bits - surplus);
            } else {
                if (!(acc % BIT_SIZE_OF_T)) {
                    xs[i] = _packed[block - 1] >> (BIT_SIZE_OF_T - n_bits);
                } else {
                    xs[i] = (_packed[block] >> ((acc % BIT_SIZE_OF_T) - n_bits)) >> (BIT_SIZE_OF_T - n_bits);
                }
            }
        }
        return xs;
    }

    void getOutNghs(uint32_t *ref) {
        if (fakeOutDegree > 0) {
            auto ref_offset = 0;
            auto out_offset = 0;
            auto head = reinterpret_cast<uint32_t *>(out)[0];
            uint32_t prev_scalar[1] = { head };
            ref[ref_offset++] = head;
            out_offset += sizeof(uint32_t);
            if (fakeOutDegree > 1) {
                auto n_blocks = (fakeOutDegree - 1) / 8;

                if (n_blocks) {
                    out_offset += (n_blocks + 1) / 2;

                    auto prev = _mm256_broadcastd_epi32(_mm_load_si128(reinterpret_cast<__m128i *>(prev_scalar)));
                    alignas(256) uint32_t xs[8];
                    for (auto i = 0; i < n_blocks; i++) {
                        auto s = ((out[sizeof(uint32_t) + (i / 2)] >> (i % 2) * 4) & 0b00001111) * 2;
                        auto curr = _mm256_load_si256(reinterpret_cast<__m256i *>(unpack(out + out_offset, s, xs)));
                        _mm256_storeu_si256(reinterpret_cast<__m256i *>(ref + ref_offset),
                            _mm256_add_epi32(prev, curr));
                        out_offset += s;
                        ref_offset += 8;
                        *prev_scalar = ref[ref_offset - 1];
                        prev = _mm256_broadcastd_epi32(_mm_load_si128(reinterpret_cast<__m128i *>(prev_scalar)));
                    }
                }

                if (fakeOutDegree - ref_offset > 0) {
                    auto flag_idx = out_offset;
                    out_offset++;
                    for (; ref_offset < fakeOutDegree; ref_offset++) {
                        if (out[flag_idx] >> (7 - (fakeOutDegree - ref_offset)) & 0b00000001) {
                            *prev_scalar += reinterpret_cast<uint32_t *>(out + out_offset)[0];
                            out_offset += 4;
                        } else {
                            *prev_scalar += reinterpret_cast<uint16_t *>(out + out_offset)[0];
                            out_offset += 2;
                        }
                        ref[ref_offset] = *prev_scalar;
                    }
                }
            }
        } else {
            return;
        }
    }
    template<typename Func>
    void traverseOutNgh(Func f) {
        if (fakeOutDegree > 0) {
            auto ref_offset = 0;
            auto out_offset = 0;
            auto head = reinterpret_cast<uint32_t *>(out)[0];
            uint32_t prev_scalar[1] = { head };
            f(head);
            ref_offset++;
            out_offset += sizeof(uint32_t);
            if (fakeOutDegree > 1) {
                auto n_blocks = (fakeOutDegree - 1) / 8;

                if (n_blocks) {
                    out_offset += (n_blocks + 1) / 2;

                    auto prev = _mm256_broadcastd_epi32(_mm_load_si128(reinterpret_cast<__m128i *>(prev_scalar)));
                    alignas(256) uint32_t xs[8];
                    for (auto i = 0; i < n_blocks; i++) {
                        auto s = ((out[sizeof(uint32_t) + (i / 2)] >> (i % 2) * 4) & 0b00001111) * 2;
                        auto curr = _mm256_load_si256(reinterpret_cast<__m256i *>(unpack(out + out_offset, s, xs)));
                        _mm256_storeu_si256(reinterpret_cast<__m256i *>(xs),
                            _mm256_add_epi32(prev, curr));
                        for (auto j = 0; j < 8; j++) {
                            f(xs[i]);
                        }
                        ref_offset += 8;
                        out_offset += s;
                        *prev_scalar = xs[7];
                        prev = _mm256_broadcastd_epi32(_mm_load_si128(reinterpret_cast<__m128i *>(prev_scalar)));
                    }
                }

                if (fakeOutDegree - ref_offset > 0) {
                    auto flag_idx = out_offset;
                    out_offset++;
                    for (; ref_offset < fakeOutDegree; ref_offset++) {
                        if (out[flag_idx] >> (7 - (fakeOutDegree - ref_offset)) & 0b00000001) {
                            *prev_scalar += reinterpret_cast<uint32_t *>(out + out_offset)[0];
                            out_offset += 4;
                        }
                        else {
                            *prev_scalar += reinterpret_cast<uint16_t *>(out + out_offset)[0];
                            out_offset += 2;
                        }
                        f(*prev_scalar);
                    }
                }
            }
        }
        else {
            return;
        }
    }
#endif

    void getOutNgh(uintE *data, uint64_t size){ decode<uintE>(out,size,data);}
    void setOutNeighbors(uint8_t *_i) { out = _i; }
    void setInNeighbors(uint8_t *_i) { in = _i; }
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
