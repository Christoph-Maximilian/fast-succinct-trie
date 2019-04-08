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

#ifndef _BITMAPRANK_H_
#define _BITMAPRANK_H_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <stdint.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <vector>
#include "shared.h"

class BitmapRank {
public:
    const int kWordSize = 64;
    const int kBasicBlockSize = 512;
    const int kBasicBlockBits = 9;
    const int kBasicBlockMask = kBasicBlockSize - 1;
    const int kWordCountPerBasicBlock = kBasicBlockSize / kWordSize;

    BitmapRank() { pCount_ = 0; }    
    virtual uint64_t rank(uint64_t pos) = 0;
    uint64_t pCount() { return pCount_; }
    
protected:
    uint64_t pCount_;
};

class BitmapRankPoppy: public BitmapRank {
public:
    BitmapRankPoppy(uint64_t* bits, uint64_t nbits);
    ~BitmapRankPoppy(){
        delete[] bits_;
        delete[] rankLUT_;
    };

    uint64_t rank(uint64_t pos);


    uint64_t* getBits();
    uint64_t getNbits();
    uint64_t getMem();

    friend class FST;
    friend class FSTIter;
    
private:
    uint64_t* bits_;
    std::vector<uint64_t > bitsVector;
    uint64_t  nbits_;
    uint64_t  mem_;

    uint64_t* rankLUT_;
    std::vector<uint64_t > rankLUTVector;
    uint64_t  basicBlockCount_;
};

#endif /* _BITMAPRANK_H_ */
