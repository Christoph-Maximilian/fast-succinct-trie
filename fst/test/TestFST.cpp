//************************************************
// FST unit tests
//************************************************
#include "gtest/gtest.h"
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <string>
#include <bitset>
#include "FST.hpp"
#include "s2/s2cell_id.h"
#include <stdio.h>

#define TEST_SIZE 234369
#define RANGE_SIZE 10
#define CELLS 128
//#define PRINT_KEY_STRING

using namespace std;

const string testFilePath = "test/bulkload_sort";
const string testFileExamplePath = "test/bulkload_example";
const string testPolygonIdsPath = "test/polygon_ids";
const string testPointsIdsPath = "test/simple_point_test";

class UnitTest : public ::testing::Test {
public:
    virtual void SetUp() {}

    virtual void TearDown() {}
};

inline void printStatFST(FST *index) {
    cout << "mem = " << index->mem() << "\n";

    cout << "cMemU = " << index->cMemU() << "\n";
    cout << "tMemU = " << index->tMemU() << "\n";
    cout << "oMemU = " << index->oMemU() << "\n";
    cout << "keyMemU = " << index->keyMemU() << "\n";
    cout << "valueMemU = " << index->valueMemU() << "\n";

    cout << "cMem = " << index->cMem() << "\n";
    cout << "tMem = " << index->tMem() << "\n";
    cout << "sMem = " << index->sMem() << "\n";
    cout << "keyMem = " << index->keyMem() << "\n";
    cout << "valueMem = " << index->valueMem() << "\n";
}

inline int loadFile(string filePath, vector<string> &keys, vector<uint64_t> &values) {
    ifstream infile(filePath);
    string op;
    string key;
    uint64_t count = 0;
    int longestKeyLen = 0;
    while (count < TEST_SIZE && infile.good()) {
        auto tmp = key;
        infile >> key; //subject to change
        if (key == tmp) { break; }
        keys.push_back(key);
        values.push_back(count);
        if (key.length() > longestKeyLen) { longestKeyLen = key.length(); }
        count++;
    }
    return longestKeyLen;
}

std::string cell_id_to_string(const S2CellId cell_id, bool polygon = false) {
    // clear the trailing 1
    auto level = cell_id.level();
    uint64_t id = cell_id.id();
    // IMPORTANT: We do not kick the last 1, since it is important later
    // id &= ~(1UL << (64 - level * 2 - 4));
    std::bitset<64> cell_bits(id);
    //remove the face bits
    id = id << 3;
    std::bitset<64> cell_bits_without_face(id);

    auto used_bytes = static_cast<uint64_t >((2 * level + 7) / 8);
    char key[8];
    reinterpret_cast<uint64_t *>(key)[0] = __builtin_bswap64(id);

#ifdef PRINT_KEY_STRING
    std::cout << cell_bits << std::endl;
    std::cout << cell_bits_without_face << std::endl;
    std::cout << "Cell ID: " << cell_id.id() << "\nKEY: ";
    for (auto i = 0; i < used_bytes; i++) {
        std::cout << +key[i] << ":";
    }
    std::cout << std::endl;
#endif

    return std::string(key, used_bytes);
}

int loadPolygonIdsFile(vector<string> &keys, vector<uint64_t> &values, std::string file) {
    ifstream infile(file);
    std::string op;
    std::string cell_id;
    uint64_t count = 0;
    int longestKeyLen = 0;
    while (count < CELLS && infile.good()) {
        std::string tmp = cell_id;
        infile >> cell_id; //subject to change
        if (cell_id == tmp) { break; }
        std::bitset<64> baz(cell_id);
        S2CellId cell(baz.to_ulong());
        auto keystring = cell_id_to_string(cell);
        keys.push_back(keystring);
        values.push_back(count);
        if (cell_id.length() > longestKeyLen) { longestKeyLen = keystring.length(); }
        count++;
    }
    return longestKeyLen;
}

inline int loadFile_ptrValue(string filePath, vector<string> &keys, vector<uint64_t> &values) {
    ifstream infile(filePath);
    string op;
    string key;
    uint64_t count = 0;
    int longestKeyLen = 0;
    while (count < TEST_SIZE && infile.good()) {
        infile >> key; //subject to change
        keys.push_back(key);
        //values.push_back(count);
        values.push_back((uint64_t) (const_cast<char *>(keys[count].c_str())));
        if (key.length() > longestKeyLen)
            longestKeyLen = key.length();
        count++;
    }

    return longestKeyLen;
}

inline int loadMonoInt(vector<uint64_t> &keys) {
    for (uint64_t i = 0; i < TEST_SIZE; i++)
        keys.push_back(i);
    return sizeof(uint64_t);
}

inline int loadRandInt(vector<uint64_t> &keys) {
    srand(0);
    for (uint64_t i = 0; i < TEST_SIZE; i++) {
        uint64_t r = rand();
        keys.push_back(r);
    }
    sort(keys.begin(), keys.end());
    return sizeof(uint64_t);
}


//*****************************************************************
// FST TESTS
//*****************************************************************

TEST_F(UnitTest, ScanTest) {
    vector<string> keys;
    vector<uint64_t> values;
    int longestKeyLen = loadFile(testFilePath, keys, values);

    FST *index = new FST(2);
    index->load(keys, values, longestKeyLen);

    printStatFST(index);

    FSTIter iter(index);
    for (int i = 0; i < TEST_SIZE - 1; i++) {
        if (i > 0 && keys[i].compare(keys[i - 1]) == 0)
            continue;
        ASSERT_TRUE(index->lowerBound((uint8_t *) keys[i].c_str(), keys[i].length(), iter));
        ASSERT_EQ(values[i], iter.value());

        for (int j = 0; j < RANGE_SIZE; j++) {
            if (i + j + 1 < TEST_SIZE) {
                ASSERT_TRUE(iter++);
                ASSERT_EQ(values[i + j + 1], iter.value());
            } else {
                ASSERT_FALSE(iter++);
                ASSERT_EQ(values[TEST_SIZE - 1], iter.value());
            }
        }
    }
}

TEST_F(UnitTest, ScanMonoIntTest) {
    vector<uint64_t> keys;
    int longestKeyLen = loadMonoInt(keys);

    FST *index = new FST();
    index->load(keys, keys);

    printStatFST(index);

    FSTIter iter(index);
    for (int i = 0; i < TEST_SIZE - 1; i++) {
        uint64_t fetchedValue;
        index->lookup(keys[i], fetchedValue);
        ASSERT_TRUE(index->lowerBound(keys[i], iter));
        ASSERT_EQ(keys[i], iter.value());

        for (int j = 0; j < RANGE_SIZE; j++) {
            if (i + j + 1 < TEST_SIZE) {
                ASSERT_TRUE(iter++);
                ASSERT_EQ(keys[i + j + 1], iter.value());
            } else {
                ASSERT_FALSE(iter++);
                ASSERT_EQ(keys[TEST_SIZE - 1], iter.value());
            }
        }
    }
}

TEST_F(UnitTest, LookupTestExample) {
    vector<string> keys;
    vector<uint64_t> values;
    int longestKeyLen = loadFile(testFileExamplePath, keys, values);

    FST *index = new FST();
    index->load(keys, values, longestKeyLen);
    index->print_csv();
    printStatFST(index);
    std::string key = "grape";
    uint64_t value;
    index->lookup(reinterpret_cast<uint8_t *>(&(key[0])), key.size(), value);
    //print values saved in fst

    uint64_t fetchedValue;


}

TEST_F(UnitTest, LookupS2Prefixes) {
    vector<string> keys;
    vector<uint64_t> values;
    int longestKeyLen = loadPolygonIdsFile(keys, values, testPolygonIdsPath);

    FST *index = new FST();
    index->load(keys, values, longestKeyLen);
    //index->print();
    //printStatFST(index);

    //print values saved in fst
    uint64_t fetchedValue;
    /*
    for (int i = 0; i < keys.size(); i++) {
        ASSERT_TRUE(index->lookup(reinterpret_cast<uint8_t *>(&(keys[i][0])), keys[i].size(), fetchedValue));
        ASSERT_EQ(values[i], fetchedValue);
    }

    std::bitset<64> point_cell_id ("1000100111000010010000101100000011000000001100000000000000000000");
    S2CellId cell(point_cell_id.to_ulong());
    for (int i = 0; i < 128; i++) {
        auto keystring = keys[i];
        ASSERT_TRUE(index->lookup(reinterpret_cast<uint8_t *>(&keystring[0]), keystring.size(), fetchedValue));
        ASSERT_EQ(i, fetchedValue);
    }
     */

    vector<string> point_keys;
    vector<uint64_t> _;
    loadPolygonIdsFile(point_keys, _, testPointsIdsPath);
    for (int i = 0; i < 128; i++) {
        auto keystring = point_keys[i];
        //std::cout << keystring << std::endl;
        ASSERT_TRUE(index->lookup(reinterpret_cast<uint8_t *>(&keystring[0]), keystring.size(), fetchedValue));
        ASSERT_EQ(i, fetchedValue);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
