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

#ifndef _BITMAPRANKF_H_
#define _BITMAPRANKF_H_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <stdint.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "shared.h"

class BitmapRankF {
public:
    const int kWordSize = 64;
    const int kBasicBlockSize = 64;
    const int kBasicBlockBits = 6;
    const int kBasicBlockMask = kBasicBlockSize - 1;
    const int kWordCountPerBasicBlock = kBasicBlockSize / kWordSize;

    BitmapRankF() { pCount_ = 0; }
    virtual uint64_t rank(uint64_t pos) = 0;
    uint64_t pCount() { return pCount_; }
    
protected:
    uint64_t pCount_;
};

class BitmapRankFPoppy: public BitmapRankF {
public:
    BitmapRankFPoppy(uint64_t* bits, uint64_t nbits);
    ~BitmapRankFPoppy() {}

    uint64_t rank(uint64_t pos);

    uint64_t* getBits();
    uint64_t getNbits();
    uint64_t getMem();

    friend class FST;
    friend class FSTIter;
    
private:
    uint64_t* bits_;
    uint64_t  nbits_;
    uint64_t  mem_;

    uint32_t* rankLUT_;
    uint64_t  basicBlockCount_;
};

#endif /* _BITMAPRANKF_H_ */
