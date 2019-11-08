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
#include <iterator>     // std::back_inserter
#include <algorithm>

#define PRINT_KEY_STRING
#define TEST_SIZE 234369
#define RANGE_SIZE 10
#define CELLS 128
//#define PRINT_KEY_STRING

using namespace std;

const string testFilePath = "test/bulkload_sort";
const string testFileExamplePath = "test/bulkload_example";
const string testPolygonIdsPath = "test/polygon_ids_simplified";
const string testPointsIdsPath = "test/simple_point_test";
const string testSimplifiedPointsIdsPath = "test/point_ids_simplified";

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

std::vector<uint8_t> cell_id_to_array(const S2CellId cell_id, bool polygon = false) {
    // clear the trailing 1
    auto level = cell_id.level();
    uint64_t id = cell_id.id();
    // IMPORTANT: We do not kick the last 1, since it is important later
    // id &= ~(1UL << (64 - level * 2 - 4));
    std::bitset<64> cell_bits(id);
    //remove the face bits
    id = id << 3u;
    std::bitset<64> cell_bits_without_face(id);

    auto used_bytes = static_cast<uint64_t >((2 * level + 7) / 8);
    std::vector<uint8_t> key(used_bytes + 1);
    uint8_t shift = 64 - level * 2;
    // clear the last set bit
    id = (id >> shift) << shift;
    std::bitset<64> cell_bits_without_face_trailing(id);
    key[0] = level;
    reinterpret_cast<uint64_t *>(key.data() + 1)[0] = __builtin_bswap64(id);


#ifdef PRINT_KEY_STRING
    std::cout << cell_bits << std::endl;
    std::cout << cell_bits_without_face << std::endl;
    std::cout << cell_bits_without_face_trailing << std::endl;

    std::cout << "Cell ID: " << cell_id.id() << "\nKEY: ";
    for (auto i = 0; i < used_bytes + 1; i++) {
        std::cout << +key[i] << ":";
    }
    std::cout << std::endl;
#endif

    return key;
}

int loadPolygonIdsFile(vector<uint8_t> &keys, vector<uint64_t> &values, const std::string &file) {
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
        auto key_vector = cell_id_to_array(cell);
        std::copy(key_vector.begin(), key_vector.end(), std::back_inserter(keys));
        values.push_back(count);
        if (key_vector[0] > longestKeyLen) { longestKeyLen = key_vector[0]; }
        count++;
    }
    return longestKeyLen;
}


//*****************************************************************
// FST TESTS
//*****************************************************************

TEST_F(UnitTest, ScanTest) {
    vector<uint8_t> keys;
    vector<uint64_t> values;
    int longestKeyLen = loadPolygonIdsFile(keys, values, testPolygonIdsPath);

    vector<uint8_t> point_keys;
    vector<uint64_t> point_values;
    loadPolygonIdsFile(point_keys, point_values, testSimplifiedPointsIdsPath);

    FST *index = new FST(10);
    index->load(keys, values, longestKeyLen);

    uint64_t expected_value = 0;
    for (int i = 0; i < point_keys.size();) {
        uint64_t value;
        int keylen = point_keys[i];
        index->lookup(point_keys.data() + i + 1, keylen, value);
        assert(value == expected_value);
        expected_value++;
        i += (keylen * 2 + 7) / 8 + 1;
    }

    printStatFST(index);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
