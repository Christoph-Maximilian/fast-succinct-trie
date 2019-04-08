/**
Copyright 2013 Carnegie Mellon University

Authors: Dong Zhou, David G. Andersen and Michale Kaminsky

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Modified by Huanchen Zhang
*/

#include <assert.h>
#include <stdlib.h>
#include <algorithm>

#include "bitmap-rank.h"
#include "popcount.h"
#include "shared.h"

#include <iostream>

// used for the LOWER level (sparse) of FST
BitmapRankPoppy::BitmapRankPoppy(uint64_t *bits, uint64_t nbits) {
    bits_ = bits;
    auto bit_length_64 = nbits % 64 == 0 ? nbits / 64 : nbits / 64 + 1;
    bitsVector.resize(bit_length_64);
    std::copy(bits, bits + bit_length_64, bitsVector.begin());
    nbits_ = nbits;
    basicBlockCount_ = nbits_ / kBasicBlockSize;
    rankLUT_ = new uint64_t[basicBlockCount_];
    //rankLUTVector.resize(basicBlockCount_);

    assert(posix_memalign((void **) &rankLUT_, kCacheLineSize, basicBlockCount_ * sizeof(uint32_t)) >= 0);

    uint64_t rankCum = 0;
    for (uint64_t i = 0; i < basicBlockCount_; i++) {
        rankLUT_[i] = rankCum;
        //rankLUTVector[i] = rankCum;
        rankCum += popcountLinear(bits_,
                                  i * kWordCountPerBasicBlock,
                                  kBasicBlockSize);
    }
    //TODO - Christoph changed this part
    //rankLUT_[basicBlockCount_ - 1] = rankCum;
    //rankLUTVector[basicBlockCount_ - 1] = rankCum;

    pCount_ = rankCum;
    mem_ = nbits / 8 + basicBlockCount_ * sizeof(uint32_t);
}

uint64_t BitmapRankPoppy::rank(uint64_t pos) {
    assert(pos <= nbits_);
    uint64_t blockId = pos >> kBasicBlockBits;// == pos / 2^9 == pos / 512
    // which superblock? -> blockid << 3 == blockid * 8
    // which position within superblock? -> pos & 511 = 00000001|11111111_2
    return rankLUT_[blockId] + popcountLinear(bits_, (blockId << 3), (pos & 511));
    // 511_10 = 111111111_2
}

uint64_t *BitmapRankPoppy::getBits() {
    return bits_;
}

uint64_t BitmapRankPoppy::getNbits() {
    return nbits_;
}

uint64_t BitmapRankPoppy::getMem() {
    return mem_;
}
