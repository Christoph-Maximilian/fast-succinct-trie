#include <emmintrin.h>
#include <assert.h>
#include <string.h>

#include <vector>
#include <iostream>
#include <string>

#include <common.h>
#include "popcount.h"

using namespace std;

class FSTIter;
class FST;

//******************************************************
// Constants for FST
//******************************************************
const uint8_t TERM = 36; //$
const int CUTOFF_RATIO = 64;

//******************************************************
//Constants for rank and select
//******************************************************
const int kCacheLineSize = 64;
const int kWordSize = 64;

const int kBasicBlockSizeU = 64;
const int kBasicBlockBitsU = 6;
const int kBasicBlockMaskU = kBasicBlockSizeU - 1;
const int kWordCountPerBasicBlockU = kBasicBlockSizeU / kWordSize;

const int kBasicBlockSize = 512;
const int kBasicBlockBits = 9;
const int kBasicBlockMask = kBasicBlockSize - 1;
const int kWordCountPerBasicBlock = kBasicBlockSize / kWordSize;

const int kWordBits = 6;
const uint32_t kWordMask = ((uint32_t)1 << kWordBits) - 1;

const int skip = 64;
const int kSkipBits = 6;
const uint32_t kSkipMask = ((uint32_t)1 << kSkipBits) - 1;


//******************************************************
// Initilization functions for FST
//******************************************************
FST* load(vector<string> &keys, vector<uint64_t> &values, int longestKeyLen);
FST* load(vector<uint64_t> &keys, vector<uint64_t> &values);

//helpers
inline bool insertChar_cond(const uint8_t ch, vector<uint8_t> &c, vector<uint64_t> &t, vector<uint64_t> &s, int &pos, int &nc);
inline bool insertChar(const uint8_t ch, bool isTerm, vector<uint8_t> &c, vector<uint64_t> &t, vector<uint64_t> &s, int &pos, int &nc);


//******************************************************
// Fast Succinct Trie Class
//******************************************************
class FST {
public:
    FST();
    FST(int16_t cl, uint16_t th, int8_t fvp, int32_t lvp, uint32_t ncu, uint32_t ccu,
	uint32_t cUmem, uint32_t cUbbc, uint32_t tUmem, uint32_t tUbbc,
	uint32_t oUmem, uint32_t oUbbc, uint32_t vUm, uint32_t cmem,
	uint32_t tmem, uint32_t tbbc, uint32_t smem, uint32_t sbbc, uint32_t vm);
    virtual ~FST();

    friend inline bool insertChar_cond(const uint8_t ch, vector<uint8_t> &c, vector<uint64_t> &t, vector<uint64_t> &s, int &pos, int &nc);
    friend inline bool insertChar(const uint8_t ch, bool isTerm, vector<uint8_t> &c, vector<uint64_t> &t, vector<uint64_t> &s, int &pos, int &nc);

    friend FST* load(vector<string> &keys, vector<uint64_t> &values, int longestKeyLen);
    friend FST* load(vector<uint64_t> &keys, vector<uint64_t> &values);

    //point query
    bool lookup(const uint8_t* key, const int keylen, uint64_t &value);
    bool lookup(const uint64_t key, uint64_t &value);

    //range query
    bool lowerBound(const uint8_t* key, const int keylen, FSTIter &iter);
    bool lowerBound(const uint64_t key, FSTIter &iter);
    bool upperBound(const uint8_t* key, const int keylen, FSTIter &iter);
    bool upperBound(const uint64_t key, FSTIter &iter);


    //mem stats
    //-------------------------------------------------
    uint32_t cMemU();
    uint32_t cRankMemU();
    uint32_t tMemU();
    uint32_t tRankMemU();
    uint32_t oMemU();
    uint32_t oRankMemU();
    uint32_t valueMemU();

    uint64_t cMem();
    uint32_t tMem();
    uint32_t tRankMem();
    uint32_t sMem();
    uint32_t sSelectMem();
    uint64_t valueMem();

    uint64_t mem();

    //debug
    void printU();
    void print();

private:
    //bit/byte vector accessors
    //-------------------------------------------------
    inline uint64_t* cUbits_();
    inline uint32_t* cUrankLUT_();
    inline uint64_t* tUbits_();
    inline uint32_t* tUrankLUT_();
    inline uint64_t* oUbits_();
    inline uint32_t* oUrankLUT_();
    inline uint8_t* cbytes_();
    inline uint64_t* tbits_();
    inline uint32_t* trankLUT_();
    inline uint64_t* sbits_();
    inline uint32_t* sselectLUT_();
    inline uint64_t* valuesU_();
    inline uint64_t* values_();

    //rank and select init
    //-------------------------------------------------
    inline void cUinit(uint64_t* bits, uint32_t nbits);
    inline void tUinit(uint64_t* bits, uint32_t nbits);
    inline void oUinit(uint64_t* bits, uint32_t nbits);
    inline void tinit(uint64_t* bits, uint32_t nbits);
    inline void sinit(uint64_t* bits, uint32_t nbits);

    //rank and select support
    //-------------------------------------------------
    inline uint32_t cUrank(uint32_t pos);
    inline uint32_t tUrank(uint32_t pos);
    inline uint32_t oUrank(uint32_t pos);
    inline uint32_t trank(uint32_t pos);
    inline uint32_t sselect(uint32_t pos);

    //helpers
    //-------------------------------------------------
    inline bool isCbitSetU(uint64_t nodeNum, uint8_t kc);
    inline bool isTbitSetU(uint64_t nodeNum, uint8_t kc);
    inline bool isObitSetU(uint64_t nodeNum);
    inline bool isSbitSet(uint64_t pos);
    inline bool isTbitSet(uint64_t pos);
    inline uint64_t valuePosU(uint64_t nodeNum, uint64_t pos);
    inline uint64_t valuePos(uint64_t pos);

    inline uint64_t childNodeNumU(uint64_t pos);
    inline uint64_t childNodeNum(uint64_t pos);
    inline uint64_t childpos(uint64_t nodeNum);

    inline int nodeSize(uint64_t pos);
    inline bool simdSearch(uint64_t &pos, uint64_t size, uint8_t target);
    inline bool binarySearch(uint64_t &pos, uint64_t size, uint8_t target);
    inline bool linearSearch(uint64_t &pos, uint64_t size, uint8_t target);
    inline bool nodeSearch(uint64_t &pos, int size, uint8_t target);
    inline bool nodeSearch_lowerBound(uint64_t &pos, int size, uint8_t target);
    inline bool nodeSearch_upperBound(uint64_t &pos, int size, uint8_t target);

    inline bool binarySearch_lowerBound(uint64_t &pos, uint64_t size, uint8_t target);
    inline bool binarySearch_upperBound(uint64_t &pos, uint64_t size, uint8_t target);
    inline bool linearSearch_lowerBound(uint64_t &pos, uint64_t size, uint8_t target);
    inline bool linearSearch_upperBound(uint64_t &pos, uint64_t size, uint8_t target);

    inline bool nextItemU(uint64_t nodeNum, uint8_t kc, uint8_t &cc);
    inline bool prevItemU(uint64_t nodeNum, uint8_t kc, uint8_t &cc);

    inline bool nextLeftU(int keypos, uint64_t pos, FSTIter* iter);
    inline bool nextLeft(int keypos, uint64_t pos, FSTIter* iter);
    inline bool nextRightU(int keypos, uint64_t pos, FSTIter* iter);
    inline bool nextRight(int keypos, uint64_t pos, FSTIter* iter);

    inline bool nextNodeU(int keypos, uint64_t nodeNum, FSTIter* iter);
    inline bool nextNode(int keypos, uint64_t pos, FSTIter* iter);
    inline bool prevNodeU(int keypos, uint64_t nodeNum, FSTIter* iter);
    inline bool prevNode(int keypos, uint64_t pos, FSTIter* iter);

    //members
    //-------------------------------------------------
    int16_t cutoff_level_;
    uint16_t tree_height_;
    int8_t first_value_pos_; // negative means in valuesU_
    int32_t last_value_pos_; // negative means in valuesU_
    uint32_t nodeCountU_;
    uint32_t childCountU_;

    //D-Labels
    uint32_t cUnbits_;
    uint32_t cUpCount_;
    uint32_t cUmem_;
    uint32_t cUbasicBlockCount_;

    //D-HasChild
    uint32_t tUnbits_;
    uint32_t tUpCount_;
    uint32_t tUmem_;
    uint32_t tUbasicBlockCount_;

    //D-IsPrefixKey
    uint32_t oUnbits_;
    uint32_t oUpCount_;
    uint32_t oUmem_;
    uint32_t oUbasicBlockCount_;

    //D-values
    uint32_t valUmem_;

    //S-Labels
    uint32_t cmem_;

    //S-HasChild
    uint32_t tnbits_;
    uint32_t tpCount_;
    uint32_t tmem_;
    uint32_t tbasicBlockCount_;

    //S-LOUDS
    uint32_t  snbits_;
    uint32_t  spCount_;
    uint32_t  smem_;
    uint32_t  sselectLUTCount_;

    //S-values
    uint64_t valmem_;

    //data
    char data_[0];

    friend class FSTIter;
};

typedef struct {
    int32_t keyPos;
    int32_t valPos;
    bool isO;
} Cursor;


class FSTIter {
public:
    FSTIter();
    FSTIter(FST* idx);

    void clear ();

    inline void setVU (int keypos, uint64_t nodeNum, uint64_t pos);
    inline void setVU_R (int keypos, uint64_t nodeNum, uint64_t pos);
    inline void setKVU (int keypos, uint64_t nodeNum, uint64_t pos, bool o);
    inline void setKVU_R (int keypos, uint64_t nodeNum, uint64_t pos, bool o);
    inline void setV (int keypos, uint64_t pos);
    inline void setV_R (int keypos, uint64_t pos);
    inline void setKV (int keypos, uint64_t pos);
    inline void setKV_R (int keypos, uint64_t pos);

    string key ();
    uint64_t value ();
    bool operator ++ (int);
    bool operator -- (int);

private:
    FST* index;
    vector<Cursor> positions;

    uint32_t len;
    bool isBegin;
    bool isEnd;

    uint32_t cBoundU;
    uint64_t cBound;
    int cutoff_level;
    uint32_t tree_height;
    int8_t first_value_pos;
    uint32_t last_value_pos;

    friend class FST;
};

