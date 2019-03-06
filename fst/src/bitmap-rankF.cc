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

#include "bitmap-rankF.h"
#include "popcount.h"
#include "shared.h"

#include <iostream>

BitmapRankFPoppy::BitmapRankFPoppy(uint64_t *bits, uint32_t nbits)
{
    bits_ = bits;
    nbits_ = nbits;
    basicBlockCount_ = nbits_ / kBasicBlockSize;

    assert(posix_memalign((void **) &rankLUT_, kCacheLineSize, basicBlockCount_ * sizeof(uint32_t)) >= 0);

    uint32_t rankCum = 0;
    for (uint32_t i = 0; i < basicBlockCount_; i++) {
	rankLUT_[i] = rankCum;
	rankCum += popcountLinear(bits_, 
				  i * kWordCountPerBasicBlock, 
				  kBasicBlockSize);
    }
    rankLUT_[basicBlockCount_-1] = rankCum;

    pCount_ = rankCum;
    mem_ = nbits / 8 + basicBlockCount_ * sizeof(uint32_t);
}

uint32_t BitmapRankFPoppy::rank(uint32_t pos)
{
    assert(pos <= nbits_);
    uint32_t blockId = pos >> kBasicBlockBits;
    uint32_t offset = pos & (uint32_t)63;
    if (offset)
	return rankLUT_[blockId] + popcount(bits_[blockId] >> (64 - offset));
    else
	return rankLUT_[blockId];
}

uint64_t* BitmapRankFPoppy::getBits() {
    return bits_;
}

uint32_t BitmapRankFPoppy::getNbits() {
    return nbits_;
}

uint32_t BitmapRankFPoppy::getMem() {
    return mem_;
}
